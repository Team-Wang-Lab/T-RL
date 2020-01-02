#' Predit optimal treatment using the output from DTRtree.
#' 
#' @param treeout An object outputted from DTRtree function.
#' @param newdata New data containing H (history).
#' @export

predict.DTR<-function(treeout,newdata){
  n<-nrow(newdata)
  predicts<-rep(NA,n)
  
  # treeout is supposed to be a matrix
  # if there is no split
  if(length(treeout)==5){
    predicts<-rep(treeout[5],n)
  } else{ # if there are splits
    treeout<-as.data.frame(treeout)
    newdata<-as.data.frame(newdata)
    
    for(i in 1:n){
      nd<-1
      while(is.na(treeout$trt[treeout$node==nd])){
        if(newdata[i,treeout$X[treeout$node==nd]] <= treeout$cutoff[treeout$node==nd]){#if the node<= cutoff
          nd=2*nd #yes proceed first
        } else{
          nd=2*nd+1#then no
        }
      }
      predicts[i]<-treeout$trt[treeout$node==nd]
    }
  }
  return(predicts)
}
