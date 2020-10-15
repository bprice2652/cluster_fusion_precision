### This is RDA code from the original RidgeFusion paper.
## THe implementation is for validation likelihood selection


## Calculates an MLE for covraiance
Cov<-function(x){(dim(x)[1]-1)*cov(x)/(dim(x)[1])}

## Calculates a common covariance matrix across classes based on the MLE 
ComCov<-function(X){
  nc=lapply(X,function(x){dim(x)[1]})
  n=sum(unlist(nc))
  CL=0
  Cov2=lapply(X,Cov)
  for(l in 1:length(X)){
    CL=nc[[l]]*Cov2[[l]]
  }
  return(CL/n)
}



MyRDA<-function(X,lambda,gamma){
  p=dim(X[[1]])[2]
  ClassCov=lapply(X,Cov)
  nc2=lapply(X,function(x){dim(x)[1]})
  Common=ComCov(X)
  CovLam=lapply(ClassCov,function(x){
    return((1-lambda)*x+lambda*Common)
  })
  
  CovGamma=lapply(CovLam,function(x){
    return((1-gamma)*x+gamma*mean(diag(x))*diag(1,p,p))
  })
  Conv=TRUE
  for(l in 1:length(X)){
    if(is.finite(1/det(CovGamma[[l]]))==FALSE){Conv=FALSE}
  }
  List=list(Cov=CovGamma, converge=Conv)
  return(List)
  
}

## A quick RDA Claculation
MyRDA2<-function(Cov,Common,lambda,gamma){
  
  CovLam=lapply(Cov,function(x){
    return((1-lambda)*x+lambda*Common)
  })
  
  CovGamma=lapply(CovLam,function(x){
    return((1-gamma)*x+gamma*mean(diag(x))*diag(1,p,p))
  })
  Conv=TRUE
  List=list(Cov=CovGamma, converge=Conv)
  return(List)
  
}



##  A cross validation function for selection of RDA

myCVrda<-function(X,lambda,gamma,Folds){
  
  XTrain=list(0,0)
  XTest=list(0,0)
  S=list(0,0)
  TV=list(0,0)
  TV2=list(0,0)
  for(i in 1:length(Folds)){
    XTrain[[i]]=list(0,0)
    XTest[[i]]=list(0,0)
    for(j in 1:length(X)){
      XTrain[[i]][[j]]=X[[j]][-Folds[[i]][[j]],]
      XTest[[i]][[j]]=X[[j]][Folds[[i]][[j]],]
    }
    S[[i]]=lapply(XTrain[[i]],Cov)
    TV[[i]]=lapply(XTest[[i]],Cov)
  }
  TV2=lapply(XTrain,ComCov)
  
  ResMat=matrix(0,length(lambda),length(gamma))
  for(l in 1:length(lambda)){
    for(m in 1:length(gamma)){
      L=list(0,0)
      for(i in 1:length(Folds)){	
        L=lapply(MyRDA2(S[[i]],TV2[[i]],lambda[l],gamma[m])$Cov,solve)
        
        for(j in 1:length(X)){
          ResMat[l,m]=ResMat[l,m]+(length(XTest[[i]][[j]])/2)*(sum(diag(TV[[i]][[j]]%*%L[[j]]))-determinant(L[[j]])$mod)
        }
      }
      ResMat[l,m]=ResMat[l,m]/length(Folds)
    }
  }
  
  Max=which(ResMat==min(ResMat),arr.ind=TRUE)
  Lam=lambda[Max[1,1]]
  Gam=gamma[Max[1,2]]
  Ret=list(Lam=Lam,Gam=Gam,ResMat)
  return(Ret)
}



