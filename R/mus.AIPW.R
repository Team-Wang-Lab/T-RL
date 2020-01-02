#' Calculate weighted mean outcome with AIPW.
#' @param Y A vector of outcome of interest.
#' @param A Treatment vector.
#' @param pis.hat Estimated propensity score matrix. 
#' @param mus.reg Regression-based conditional mean outcome.
#' @export

mus.AIPW<-function(Y,A,pis.hat,mus.reg){
  class.A<-sort(unique(A))
  K<-length(class.A)
  N<-length(A)
  if(K<2 | N<2) stop("No multiple treatments or samples!")
  if(ncol(pis.hat)!=K | ncol(mus.reg)!=K | nrow(pis.hat)!=N | nrow(mus.reg)!=N) stop("Treatment, propensity or conditional means do not match!")
  
  #AIPW estimates
  mus.a<-matrix(NA,N,K)
  for(k in 1L:K){
    mus.a[,k]<-(A==class.A[k])*Y/pis.hat[,k]+(1-(A==class.A[k])/pis.hat[,k])*mus.reg[,k]
  }
  return(mus.a)
  #AIPW Y* that weights the original Y and the Y* depend on model
}