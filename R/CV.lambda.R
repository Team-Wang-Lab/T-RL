#' Cross-validation to decide minimum purity improvement for splitting.
#' @param Y A vector of outcome of interest.
#' @param A Treatment vector.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @param g.opt Optimal dynamic treatment regime.
#' @param K Number of folds.
#' @param m.method Method for calculating estimated conditional mean.
#' @param depth Maximum tree depth.
#' @param minsplit Minimal node size.
#' @export

CV.lambda<-function(Y,A,H,g.opt,lambda.pct,K=10,m.method=c("AIPW","randomForest"),depth=5,minsplit=20){
  N<-length(Y)
  # function randomly split the data into K folds, return the group indicators
  split.sample<-function(K,N){ #K=number of fold N=number of samples
    rns<-runif(N,0,1)
    group<-rep(1,N)
    i=1
    while(i < K){
      group<-group+(rns>=quantile(rns,i/K))
      i=i+1
    }
    return(group)
  }
  I.group<-split.sample(K,N)
  ppower<-0
  
  
  for(i in 1L:K){
    tree1<-DTRtree(Y[I.group != i],A[I.group != i],H[I.group != i,],m.method=m.method,
                   depth=depth,lambda.pct=lambda.pct,minsplit=minsplit)
    g.tree<-predict.DTR(tree1,H[I.group==i,])
    ppower<-ppower + mean(g.tree==g.opt[I.group==i])#This gives the final percent of correctness of the tree.
  }
  return(ppower/K)
}