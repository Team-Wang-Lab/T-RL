#' Choose the optimal treatment with given pseudo outcome matrix in one given stage.
#' @param A Treatment options given in this stage
#' @param mu.hat The final counterfactural outcome according to the treatment.
Opt.A<-function(A,mus.hat){
  class.A<-sort(unique(A))
  if(length(class.A)==1){
    trt.opt1<-class.A
    Ey.opt1<-mean(mus.hat)
  } else{
    if(length(A)!= nrow(mus.hat) || length(class.A)!= ncol(mus.hat)){
      stop("Treatment options and mean matrix dimension do not match!")
    }
    
    # pick a single best treatment for all patients
    c.means<-apply(mus.hat,2,mean)
    Ey.opt1<-max(c.means)
    trt.opt1<-class.A[which(c.means==Ey.opt1)]###Problem: if the order of class.A is different from mus.hat becuase of the sort() this will cause problems
  }
  outs<-list(Ey.opt1,trt.opt1)
  names(outs)<-c("Ey.opt1","trt.opt1")
  return(outs)
}