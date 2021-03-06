## Run from terminal in batch mode
## p is 20,50 or 100
## r is .4,.47,or .5

args <- commandArgs(trailingOnly = TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

print(p)
print(r4)
#out_name<-paste("crf_1_eps",eps,"p",p,"q1",q1,"q2",q2,".Rdata",sep="_")


out_name<-paste("rda_3_p",p,"r4",r4,".Rdata",sep="_")


### This is RDA code from the original RidgeFusion paper.
## THe implementation is for validation likelihood selection
source("~/Dropbox/cluster_inv_cov/r_code/rda.R")
setwd("~/Dropbox/cluster_inv_cov_sims/sim3")

evs<- sapply(1:p,function(m){(p-m+1)/p})
mp<-rep(1,p)
mp[1:6]=q1*100
mp[7:11]=q1*10
mp2=mp
mp2[1:6]=(q1*100)-1
mp2[7:11]=(q1*10)-1



mean1=20*log(p)/p
mean2=-10*log(p)/p
mean3=10*log(p)/p
mean4=-20*log(p)/p
set.seed(616)
library(mvtnorm)
am=vector("list",100)
best_mod=vector("list",100)
prd<-NULL
prdT=NULL
prd_like=NULL
prdT_like=NULL
S1=matrix(0,p,p)
S2=S1
S3=S1
S4=S1



for(j in 1:100){
  print(j)
  Eig<-matrix(rnorm(100*p),100,p)
  Eig=eigen(t(Eig)%*%Eig)$vec
  S1=Eig%*%(diag(q1*mp))%*%t(Eig)  
  S2=Eig%*%(diag(q1*mp2))%*%t(Eig)
  #S3=Eig%*%(diag(q2*mp))%*%t(Eig)
  #S4=Eig%*%(diag(q2*mp2))%*%t(Eig)
  for(l in 1:p){
    for(k in 1:p){
      if(abs(k-l)<=1){
        S3[l,k]=.45^(abs(l-k))
        S4[l,k]=r4^(abs(l-k))
      }
      
    }
  }
  
  
  Omeg=list(solve(S1),solve(S2),solve(S3),solve(S4))
  Means=list(mean1,mean2,mean3,mean4)
  
  X<-list()
  X[[1]]=rmvnorm(25,mean=rep(mean1,p),S1)
  X[[2]]=rmvnorm(25,mean=rep(mean2,p),S2)
  X[[3]]=rmvnorm(25,mean=rep(mean3,p),S3)
  X[[4]]=rmvnorm(25,mean=rep(mean4,p),S4)
  
  
  fold=vector("list",5)
  samp=sample(1:25)
  for(k in 1:5){
    fold[[k]]=lapply(vector("list",4),function(m){samp[((k-1)*5+1):(k*5)]})
  }
  
  
  am[[j]]<-myCVrda(X,10^seq(-10,0,.1),10^seq(-10,0,.1),Folds = fold)
  X_test=list()
  X_test[[1]]=rmvnorm(500,mean=rep(mean1,p),S1)
  X_test[[2]]=rmvnorm(500,mean=rep(mean2,p),S2)
  X_test[[3]]=rmvnorm(500,mean=rep(mean3,p),S3)
  X_test[[4]]=rmvnorm(500,mean=rep(mean4,p),S4)
  Z1=Reduce("rbind",X_test)

best_mod[[j]]<-MyRDA(X,am[[j]]$Lam,am[[j]]$Gam)
SigG=best_mod[[j]]$Cov
Omeg=lapply(SigG,solve)
M=lapply(X,function(x){apply(x,2,mean)})
Class2=rep(0,dim(Z1)[1])
Result=matrix(0,dim(Z1)[1],4)
for(i in 1:dim(Z1)[1]){
  for(j in 1:4){
    Result[i,j]=-.5*(t(Z1[i,]-M[[j]])%*%Omeg[[j]]%*%(Z1[i,]-M[[j]])-determinant(Omeg[[j]])$mod)+log(1/15)
  }
  Class2[i]=which(Result[i,]==max(Result[i,]))
}

prd=c(prd,sum(Class2==c(rep(1,500),rep(2,500),rep(3,500),rep(4,500))))
ma=t(apply(Result,1,function(m){exp(m)/sum(exp(m))}))
minesa=0
for(j in 1:4){
  minesa=minesa+sum(log(ma[(((j-1)*500+1):(j*500)),j]))
}
prd_like=c(prd_like,minesa)

#table(preds1$Predicted,c(rep(1,500),rep(2,500),rep(3,500),rep(4,500)))

newdata=Z1
K=4
n=dim(newdata)[1]
Result2=matrix(nrow=n,ncol=K)
Class2=rep(0,n)
for(i in 1:n){
  for(k in 1:K){
    Result2[i,k]=-.5*(t(newdata[i,]-mean(X[[k]]))%*%Omeg[[k]]%*%(newdata[i,]-mean(X[[k]]))-determinant(Omeg[[k]])$mod)+log(.25)
  }
  Class2[i]=which(Result2[i,]==max(Result2[i,]))	

}
ma2=t(apply(Result,1,function(m){exp(m)/sum(exp(m))}))
minesa2=0
for(j in 1:4){
  minesa2=minesa2+sum(log(ma2[(((j-1)*500+1):(j*500)),j]))
}
prdT_like=c(prdT_like,minesa2)
prdT=c(prdT,sum(Class2==c(rep(1,500),rep(2,500),rep(3,500),rep(4,500))))
}

save.image(out_name)
