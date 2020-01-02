#' Estimate conditional means for multiple stages
#' 
#' @param Y A continous outcome of interest.
#' @param As (A1, A2, ...) a matrix of treatments at multiple stages; Stage t has treatment K_t options labeled as 0, 1, ..., K_t-1.
#' @param H A matrix of covariates before assigning final treatment, excluding previous treatment variables.
#' @export
Reg.mu<-function(Y,As,H){
  if(nrow(as.matrix(As))!=nrow(as.matrix(H))) stop("Treatment and Covariates do not match in dimension!")
  Ts<-ncol(as.matrix(As)) # number of stages
  N<-nrow(as.matrix(As))
  if(Ts<0 | Ts>3) stop("Only support 1 to 3 stages!")#As as maxixmum 3 colmns
  H<-as.matrix(H)
  
  if(Ts==1L){ #one stage
    A1<-as.matrix(As)[,1]
    A1<-as.factor(A1)
    KT<-length(unique(A1)) # treatment options at last stage
    if(KT<2) stop("No multiple treatment options!")
    
    RegModel<-lm(Y ~ H*A1)
    print(RegModel)
    #Build linear model between outcome Y~historyX+trt+historyX*trt 
    #Eg C1<-c(2,4,5,2,4);A1<-c(1,2,1,1,2);B1<-matrix(c(1,2,3,10,2,5,3,9,8,2),ncol = 2);lm(C1~B1*A1);lm(C1~B1+A1%in%B1)
    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1=factor(rep(sort(unique(A1))[k],N))))
    # use model to predict when treatment =each treatment what is the predicted outcomeY
    # Thus, mus.reg is a N*KT matrix where each column is the countradictory outcome value Y* given the input H and aT.
  }
  if(Ts==2L){
    A1<-as.matrix(As)[,1];A2<-as.matrix(As)[,2]
    A1<-as.factor(A1);A2<-as.factor(A2)#require A1 A2 to be ordered by the same way
    KT<-length(unique(A2))
    if(KT<2) stop("No multiple treatment options!")
    
    RegModel<-lm(Y ~ (H + A1)*A2)
    
    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1,A2=factor(rep(sort(unique(A2))[k],N))))#Has only change the last trt
  }
  if(Ts==3L){
    A1<-as.matrix(As)[,1];A2<-as.matrix(As)[,2];A3<-as.matrix(As)[,3]
    A1<-as.factor(A1);A2<-as.factor(A2);A3<-as.factor(A3)
    KT<-length(unique(A3))
    if(KT<2) stop("No multiple treatment options!")
    
    RegModel<-lm(Y ~ (H + A1 + A2)*A3)
    
    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1,A2,A3=factor(rep(sort(unique(A3))[k],N))))
  }
  
  output<-list(mus.reg, RegModel)
  names(output)<-c("mus.reg","RegModel")
  return(output)
  #this is the Y*'s for all the possible treatments in stage K_T(last stage) given history and fixed previous trt
}