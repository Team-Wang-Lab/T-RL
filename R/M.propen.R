#' Estimate propensity score by a multinomial model.
#' @param A Treatment vector.
#' @param Xs Covariate matrix.
#' @export
M.propen<-function(A,Xs){
  if(ncol(as.matrix(A))!=1) stop("Cannot handle multiple stages of treatments together!")
  if(length(A)!= nrow(as.matrix(Xs))) stop("A and Xs do not match in dimension!")
  if(length(unique(A))<=1) stop("Treament options are insufficient!")
  class.A<-sort(unique(A))#class.A=Unique treatments
  
  require(nnet)
  s.data<-data.frame(A,Xs)
  # multinomial regression with output suppressed
  model<-capture.output(mlogit<-multinom(A ~., data=s.data))
  s.p<-predict(mlogit,s.data,"probs")#using model to predic the traning data and get probabilities of each outcome
  if(length(class.A)==2){
    s.p<-cbind(1-s.p,s.p)
  }
  colnames(s.p)<-paste("pi=",class.A,sep="")
  
  return(s.p)#matrix of N*length(unique(A))
}