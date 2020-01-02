#' Tree-based Reinforcement Learning for estimating optimal DTR.
#' 
#' a tree-based reinforcement learning (T-RL) method to directly
#' estimate optimal DTRs in a multi-stage multi-treatment setting. At
#' each stage, T-RL builds an unsupervised decision tree that directly handles
#' the problem of optimization with multiple treatment comparisons, through a
#' purity measure constructed with augmented inverse probability weighted estimators.
#' 
#' @param Y A vector of outcome of interest.
#' @param A A vector of observed treatment options. 
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param pis.hat Estimated propensity score matrix. 
#' @param m.method Method for calculating estimated conditional mean.
#' @param mus.reg Regression-based conditional mean outcome.
#' @param depth Maximum tree depth.
#' @param lambda.pct Minimal percent change in purity measure for split.
#' @param minsplit Minimal node size.
#' @param lookahead Whether or not to look into a further step of splitting to find the best split.
#' @export
DTRtree<-function(Y,A,H,pis.hat=NULL,m.method=c("AIPW","randomForest"),
                  mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F){
  # initialization
  # indicator for subset data
  n<-length(Y)#number of people
  I.node<-rep(1,n)#indicator of nodes
  class.A<-sort(unique(A))
  output<-matrix(NA,1,5)
  colnames(output)<-c("node","X","cutoff","mEy","trt")
  
  # estimate mus.hat if not given
  if(m.method[1]=="AIPW"){
    # estimate propenstiy matrix if not given, using all data
    # same propensity for all subset data
    if(is.null(pis.hat)) pis.hat<-M.propen(A=A,Xs=H)
    if(is.null(mus.reg)) mus.reg<-Reg.mu(Y=Y,As=A,H=H)$mus.reg
    mus.hat<-mus.AIPW(Y=Y,A=A,pis.hat=pis.hat,mus.reg=mus.reg)
  } else if(m.method[1]=="randomForest"){
    require(randomForest)
    RF<-randomForest(Y~., data=data.frame(A,H))
    mus.hat<-matrix(NA,n,length(class.A))
    for(i in 1L:length(class.A)) mus.hat[,i]<-predict(RF,newdata=data.frame(A=rep(class.A[i],n),H))
  } else{
    stop("The method for estimating conditional means is not available!")
  }
  
  # expected outcome at root
  root<-Opt.A(A,mus.hat)
  Ey0<-root$Ey.opt1
  
  # split if improved at least lambda, as a percent of Ey0
  lambda<-abs(Ey0)*lambda.pct
  
  for(k in 1L:depth){#depth is the most number of split to reach one terminal node
    output<-rbind(output,matrix(NA,2^k,5))#2^k??????? originally output=1*5
    output[,1]<-1L:(2^(k+1)-1)# this does not equal to the number of rows(2^k+1)
    if(k==1L){
      # apply lookahead to the first split, the most important split
      # only to first split so as to save computation time
      # use a larger minsplit for the first split
      if(lookahead){
        best.H.1<-best.H.lh(H=H,A=A,mus.hat=mus.hat,minsplit=0.15*n)
      } else{
        best.H.1<-best.H(H=H,A=A,mus.hat=mus.hat,minsplit=minsplit)
      }
      if(is.null(best.H.1)==F && best.H.1$mEy.opt1>Ey0+lambda){#meet the split criteria
        output[k,-1]<-c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1, NA)
        I.node[I.node==k & H[,best.H.1$X] <= best.H.1$X.subset]<-2*k
        output[2*k,-1]<-c(NA,NA,NA,best.H.1$trt.L)
        I.node[I.node==k & H[,best.H.1$X] > best.H.1$X.subset]<-2*k+1
        output[2*k+1,-1]<-c(NA,NA,NA,best.H.1$trt.R)
      } else{
        output[k,4:5]<-c(root$Ey.opt1,root$trt.opt1)
        break
      }
    } else{
      for(j in (2^(k-1)):(2^k-1)){
        if(!is.na(output[trunc(j/2),2])){
          best.H.j<-best.H(H=H[I.node==j,],A=A[I.node==j],mus.hat=mus.hat[I.node==j,],minsplit=minsplit)
          if(is.null(best.H.j)==F && best.H.j$mEy.opt1>output[trunc(j/2),4]+lambda){
            output[j,-1]<-c(best.H.j$X, best.H.j$X.subset, best.H.j$mEy.opt1, NA)
            I.node[I.node==j & H[,best.H.j$X] <= best.H.j$X.subset]<-2*j
            output[2*j,-1]<-c(NA,NA,NA,best.H.j$trt.L)
            I.node[I.node==j & H[,best.H.j$X] > best.H.j$X.subset]<-2*j+1
            output[2*j+1,-1]<-c(NA,NA,NA,best.H.j$trt.R)
          }
        }
      }
      if(sum(is.na(output[(2^(k-1)):(2^k-1),2]))==2^(k-1)) break
    }
  }
  output<-output[!is.na(output[,2]) | !is.na(output[,5]),]
  return(output)
}
