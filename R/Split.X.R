#' Split a parent node into two child nodes by a covariate.
#' Split the parent node in order to maximize the outcome of interest by signing different treatment in child nodes. Calculate new outcome means for each child node.
#' @param X Covariate matrix.
#' @param A A vector of observed treatments.
#' @param mus.hat Estimated conditional mean outcome.
#' @param minsplit Minimal node size.
#'  
Split.X<-function(X,A,mus.hat,minsplit=20){#X is a covariate for patients; minsplit is the minimum cases in each node
  n<-length(X)
  X.val<-unique(X)
  n.X.val<-length(X.val)
  class.A<-sort(unique(A))
  
  if(n < 2*minsplit || n.X.val<2L || length(class.A)<2L) return(NULL)
  
  if(is.numeric(X)==TRUE || is.ordered(X)==TRUE || n.X.val==2L){# is.ordered=ordered factorial feature
    X.val<-sort(X.val)
    # reduce computation by using quantiles
    if(n.X.val>100L){
      X.val<-quantile(X,1:100/100)#change the unique x to quantiles of x(only test 100 possible x as candidates)
      n.X.val<-100L
    }
    Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.val-1)#initalize E(Y|optX)
    for(i in 1L:(n.X.val-1)){
      left<-which(X<=X.val[i])#left<- index of X's that is less than the evaluated optX candidate
      if(length(left)>=minsplit && length(left)<=n-minsplit){#Make sure after split the resulting two nodes has cases more than minsplit
        
        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
        
        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1 #Average the optimum Ey for each cutoff point
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L){#check if one of the EY for candidate cutoffs is non NA
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]#take the minimum cutoff
      X.cutoff<-X.val[cutoff1]
      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]
      
      output<-data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }
    
  }
  
  if(is.numeric(X)==F && is.ordered(X)==F && n.X.val>2L){
    #this is the case when X is not numerical or ordered catagorical or catagorical with 2 classes== catagorical with more than 2 classes
    n.X.combo<-2^(n.X.val-1)-1#Assume there are c classes for X, each class we can either pick or not pick. -1 for not pick any class. -1 for power for the ????
    X.combo<-combn(X.val,1)#???????unique(X.val)?
    if(n.X.val>3L && n.X.val%%2==1L){#EVEN AND ODD CASES
      for(k in 2L:(n.X.val-1)/2) X.combo<-combine.mat(X.combo,combn(X.val,k))
    }
    if(n.X.val>3L && n.X.val%%2==0L){
      for(k in 2L:(n.X.val/2)){
        if(k<(n.X.val/2)) X.combo<-combine.mat(X.combo,combn(X.val,k))
        if(k==(n.X.val/2)){
          temp.mat<-combn(X.val[-1],k-1)
          first.row<-rep(X.val[1],ncol(temp.mat))
          X.combo<-combine.mat(X.combo,rbind(first.row,temp.mat))
        }
      }
    }
    
    Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.combo)
    for(i in 1L:n.X.combo){
      left<-which(X %in% X.combo[,i])
      if(length(left)>=minsplit && length(left)<=n-minsplit){
        
        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])
        
        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L){
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
      X.subset<-X.combo[,cutoff1]
      # change a vector into a single string while removing NA's
      X.subset<-paste(X.subset[!is.na(X.subset)], collapse=" ")
      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]
      
      output<-data.frame(X.subset, mEy.opt1, trt.L, trt.R)#data.frame only has one row
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }
  }
  return(output)
  #RETURN all the avg Y*|cutoff=each X
}