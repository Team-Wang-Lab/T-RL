# T-RL

Tree-based reinforcement learning for estimating optimal dynamic treatment regimes.

# Download 

```{r}
library(devtools)
install_github("Team-Wang-Lab/T-RL")
```

# Example

A simulation example for 2 stage and 3 treatment per stage.

```{r}
library(rpart)
library(randomForest)

########################################################
# simulation Stage 2 Treatment 3

############################################
# Functions for simulation.
# function to summarize simulation results
summary2<-function(x){#x is a numerical vector or matrix
  s1<-summary(x)
  sd2<-function(y){
   return(sd(y,na.rm=T))
  }
  if(is.matrix(x)){
    SD<-apply(x,2,sd2)
  } else{
    SD<-sd2(x)
  }
  s2<-list(s1,SD)
  names(s2)<-c("summary","SD")#return summary and sd to each column
  return(s2)
}


# function to sample treatment A
# input matrix.pi as a matrix of sampling probabilities, which could be non-normalized
A.sim<-function(matrix.pi){
  N<-nrow(matrix.pi) # sample size
  K<-ncol(matrix.pi) # treatment options
  if(N<=1 | K<=1) stop("Sample size or treatment options are insufficient!")
  if(min(matrix.pi)<0) stop("Treatment probabilities should not be negative!")

  # normalize probabilities to add up to 1 and simulate treatment A for each row
  probs<-t(apply(matrix.pi,1,function(x){x/sum(x,na.rm = TRUE)}))
  A<-apply(probs,1,function(x) sample(0:(K-1),1,prob = x))
  return(A)
}
############## Simulation start ####################################
N<-500 # sample size of training data
N2<-1000 # sample size of test data
iter<-5 # replication

select1<-select2<-selects<-rep(NA,iter) # percent of optimality
EYs<-rep(NA,iter) # estimated mean counterfactual outcome

for(i in 1:iter){
  set.seed(i)
  x1<-rnorm(N)
  x2<-rnorm(N)
  x3<-rnorm(N)
  x4<-rnorm(N)
  x5<-rnorm(N)
  X0<-cbind(x1,x2,x3,x4,x5)

############### stage 1 data simulation ##############
# simulate A1, stage 1 treatment with K1=3
pi10<-rep(1,N); pi11<-exp(0.5*x4+0.5*x1); pi12<-exp(0.5*x5-0.5*x1)
# weights matrix
matrix.pi1<-cbind(pi10,pi11,pi12)

A1<-A.sim(matrix.pi1)
class.A1<-sort(unique(A1))

# propensity stage 1
pis1.hat<-M.propen(A1,cbind(x1,x4,x5))
#   pis1.hat<-M.propen(A1,rep(1,N))

# simulate stage 1 optimal g1.opt
g1.opt<-(x1>-1)*((x2>-0.5)+(x2>0.5))
# stage 1 outcome
R1 <- exp(1.5+0.3*x4-abs(1.5*x1-2)*(A1-g1.opt)^2) + rnorm(N,0,1)

############### stage 2 data simulation ##############
# A2, stage 2 treatment with K2=3
pi20<-rep(1,N); pi21<-exp(0.2*R1-0.5); pi22<-exp(0.5*x2)
matrix.pi2<-cbind(pi20,pi21,pi22)
A2<-A.sim(matrix.pi2)
class.A2<-sort(unique(A2))

# propensity stage 2
pis2.hat<-M.propen(A2,cbind(R1,x2))
#   pis2.hat<-M.propen(A2,rep(1,N))

# optimal g2.opt
g2.opt<-(x3>-1)*((R1>0)+(R1>2))
# stage 2 outcome R2
R2<-exp(1.18+0.2*x2-abs(1.5*x3+2)*(A2-g2.opt)^2)+rnorm(N,0,1)

############### stage 2 Estimation ###############################
# Backward induction
###########################################
  
  # conditional mean model using linear regression 
  REG2<-Reg.mu(Y=R2,As=cbind(A1,A2),H=cbind(X0,R1))
  mus2.reg<-REG2$mus.reg
  
  #################
  # DTRtree
  # input: outcome Y, treatment A, covariate history H, propensity pis.hat
  # lambda.pct is minimum purity improvement as a percent of the estimated counterfactaul mean at root node without splitting
  # minsplit is minimum node size
  
  tree2<-DTRtree(R2,A2,H=cbind(X0,A1,R1),pis.hat=pis2.hat,mus.reg=mus2.reg,lambda.pct=0.02,minsplit=max(0.05*N,20))

############### stage 1 Estimation ################################
# calculate pseudo outcome (PO)

# expected optimal stage 2 outcome 
E.R2.tree<-rep(NA,N)

# estimated optimal regime
g2.tree<-predict_DTR(tree2,newdata=data.frame(X0,A1,R1))

## use observed R2 + E(loss), modified Q learning as in Huang et al.2015

# random forest for the estimated mean
RF2<-randomForest(R2~., data=data.frame(A2,X0,A1,R1))
mus2.RF<-matrix(NA,N,length(class.A2))
for(d in 1L:length(class.A2)) mus2.RF[,d]<-predict(RF2,newdata=data.frame(A2=rep(class.A2[d],N),X0,A1,R1))

for(m in 1:N){
  E.R2.tree[m]<-R2[m] + mus2.RF[m,g2.tree[m]+1]-mus2.RF[m,A2[m]+1]
}

# pseudo outcomes
PO.tree<-R1+E.R2.tree

############################################
# DTRtree
tree1<-DTRtree(PO.tree,A1,H=X0,pis.hat=pis1.hat,lambda.pct=0.02,minsplit=max(0.05*N,20))

############################################
# prediction using new data
############################################
  set.seed(i+10000)
  x1<-rnorm(N2)
  x2<-rnorm(N2)
  x3<-rnorm(N2)
  x4<-rnorm(N2)
  x5<-rnorm(N2)
  X0<-cbind(x1,x2,x3,x4,x5)
# true optimal for regime at stage 1
  g1.opt<-(x1>-1)*((x2>-0.5)+(x2>0.5))

  R1<-exp(1.5+0.3*x4)+rnorm(N2,0,1)

  R2<-exp(1.18+0.2*x2)+rnorm(N2,0,1)#has noting to do with a1
  
####### stage 1 prediction #######

  # predict selection %

  g1.tree<-predict_DTR(tree1,newdata=data.frame(X0))
  
  select1[i]<-mean(g1.tree==g1.opt)
  
  R1.tree<-exp(1.5+0.3*x4-abs(1.5*x1-2)*(g1.tree-g1.opt)^2)+rnorm(N2,0,1)

####### stage 2 prediction #######

  g2.tree<-predict_DTR(tree2,newdata=data.frame(X0,A1=g1.tree,R1=R1.tree))
 # true optimal for regime at stage 2
  g2.opt.tree<-(x3>-1)*((R1.tree>0)+(R1.tree>2))

  select2[i]<-mean(g2.tree==g2.opt.tree)

  selects[i]<-mean(g1.tree==g1.opt & g2.tree==g2.opt.tree)

# predict R2
  R2.tree<-exp(1.18+0.2*x2-abs(1.5*x3+2)*(g2.tree-g2.opt.tree)^2)+rnorm(N2,0,1)
 
  EYs[i]<-mean(R1.tree+R2.tree)
  print(tree1);print(tree2)
  if(i%%50==0) print(i)
}

summary2(select1)
summary2(select2)
summary2(selects)

summary2(EYs)


```
