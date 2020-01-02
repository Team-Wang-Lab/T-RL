#' Search for best split by lookahead
#' Compare single best split with further split to calculate the final cutoff.
#' @param X Covariate matrix.
#' @param A A vector of observed treatment options.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param mus.hat Estimated conditional mean outcome.
#' @param minsplit Minimal node size.
Split.X.lh<-function(X,A,H,mus.hat,minsplit=20){
  n<-length(X)
  X.val<-unique(X)
  n.X.val<-length(X.val)
  class.A<-sort(unique(A))
  
  if(n < 2*minsplit || n.X.val<2L || length(class.A)<2L) return(NULL)
  
  if(is.numeric(X)==T || is.ordered(X)==T || n.X.val==2L){
    X.val<-sort(X.val)
    # reduce computation by using 1% - 100% quantiles
    if(n.X.val>100L){
      X.val<-quantile(X,1:100/100)
      n.X.val<-100L
    }
    Ey.opt1.X<-Ey.opt2.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.val-1)
    for(i in 1L:(n.X.val-1)){
      left<-which(X<=X.val[i])
      if(length(left) >= minsplit && length(left) <= n-minsplit){
        # further split left and right
        # lookahead one step
        best.H.L<-best.H.R<-NULL
        if(length(left)>= 2*minsplit){#After the split left still have room for the next split
          best.H.L<-best.H(H=H[left,],A=A[left],
                           mus.hat=mus.hat[left,which(class.A %in% unique(A[left]))],minsplit=minsplit)#output the best variable to split output cutoff, mean,left, right optimal Y
        }
        if(n-length(left)>= 2*minsplit){#After the split right still have room for the next split
          best.H.R<-best.H(H=H[-left,],A=A[-left],
                           mus.hat=mus.hat[-left,which(class.A %in% unique(A[-left]))],minsplit=minsplit)
        }
        # if unable to further split, then pick 1 treatment
        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])#optimal Y in all left area
        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])#optimal Y in all right area
        
        trt.left.X[i]<-left.est$trt.opt1 #option
        trt.right.X[i]<-right.est$trt.opt1
        
        if(!is.null(best.H.L)){
          left.Ey<-best.H.L$mEy.opt1
        } else{
          left.Ey<-left.est$Ey.opt1 # if cannot further split assign to the current value.
        }
        if(!is.null(best.H.R)){
          right.Ey<-best.H.R$mEy.opt1
        } else{
          right.Ey<-right.est$Ey.opt1
        }
        # the max Ey after assigning 2 treatments to two child nodes is used for later splits
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
        
        # the max Ey if lookahead is used to choose the optimal X and its threshould
        Ey.opt2.X[i]<-length(left)/n*left.Ey+(1-length(left)/n)*right.Ey
      }
    }
    # pick the best split cutoff of X
    if(sum(!is.na(Ey.opt1.X))>0L & sum(!is.na(Ey.opt2.X))>0L){
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      mEy.opt2<-max(Ey.opt2.X, na.rm=T)
      # based on lookahead
      cutoff<-which(Ey.opt2.X==mEy.opt2)[1]
      #If can further split use the second split as the final cutoff.
      X.cutoff<-X.val[cutoff]
      trt.L<-trt.left.X[cutoff]
      trt.R<-trt.right.X[cutoff]
      
      output<-data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }
  }else{
    stop("Lookahead currently only supports numerical or ordinal covariates!")
  }
  
  return(output)
}
