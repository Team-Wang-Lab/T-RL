#' Search for best covariate for splitting based on lookahead.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param A A vector of observed treatment options.
#' @param mus.hat Estimated conditional mean outcome.
#' @param minsplit Minimal node size.
best.H.lh<-function(H,A,mus.hat,minsplit=20){
  p<-ncol(H)
  output<-as.data.frame(matrix(NA,p,5))
  output[,1]<-1:p
  colnames(output)<-c("X","X.subset","mEy.opt1","trt.L","trt.R")
  
  for(i in 1:p){
    split.i<-Split.X.lh(X=H[,i],A=A,H=H,mus.hat=mus.hat,minsplit=minsplit)
    
    if(!is.null(split.i)) output[i,-1]<-split.i
  }
  if(sum(!is.na(output$mEy.opt1))>0L){
    max.p<-which(output$mEy.opt1==max(output$mEy.opt1,na.rm=T))[1]
    opt.output<-output[max.p,]
    if(opt.output$trt.L==opt.output$trt.R){
      return(NULL)
    } else{
      return(opt.output)
    }
  } else{
    return(NULL)
  }
}