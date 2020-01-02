#' Combine matrices with different dimensions, add NA for additional rows/columns
#' Used in outputting the DTR tree.
#' @param m1 Matrix 1.
#' @param m2 Matrix 2.
combine.mat<-function(m1,m2,by="column"){
  nrow1<-nrow(m1);ncol1<-ncol(m1)
  nrow2<-nrow(m2);ncol2<-ncol(m2)
  if(by=="column"){
    combine<-matrix(NA,max(nrow1,nrow2),ncol1+ncol2)
    combine[1:nrow1,1:ncol1]<-m1
    combine[1:nrow2,(ncol1+1):(ncol1+ncol2)]<-m2
  }
  if(by=="row"){
    combine<-matrix(NA,nrow1+nrow2,max(ncol1,ncol2))
    combine[1:nrow1,1:ncol1]<-m1
    combine[(nrow1+1):(nrow1+nrow2),1:ncol2]<-m2
  }
  return(combine)
}