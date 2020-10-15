## Run in batch mode from terminal

## p is 20,50 or 100
## r is .4,.47,or .5
args <- commandArgs(trailingOnly = TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

print(p)
print(r4)
#out_name<-paste("crf_1_eps",eps,"p",p,"q1",q1,"q2",q2,".Rdata",sep="_")


source("cluster_da.R")


out_name<-paste("oracle_3_p",p,"r4",r4,".Rdata",sep="_")
#setwd("~/Dropbox/cluster_inv_cov_sims/sim3")


library(parallel)
library(cluster)
library(RidgeFusion)
library(JGL)
library(compiler)


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


#am[[j]]<-clust_inv_cov_cv(X = X,lambda1 = 10^seq(-3,1,.5),lambda2=10^seq(5,7,.5),Q_set=c(2,3),method="ridgeFusion",Fold=fold,ncores_Q = 3,ncores_cv = 5)

X_test=list()
X_test[[1]]=rmvnorm(500,mean=rep(mean1,p),S1)
X_test[[2]]=rmvnorm(500,mean=rep(mean2,p),S2)
X_test[[3]]=rmvnorm(500,mean=rep(mean3,p),S3)
X_test[[4]]=rmvnorm(500,mean=rep(mean4,p),S4)
X_test1=Reduce("rbind",X_test)

#best_mod[[j]]<-clust_qda(X,lambda1=am[[j]]$best_triplet[2],lambda2=am[[j]]$best_triplet[2],Q=am[[j]]$best_triplet[1],method = "ridgeFusion")
#preds1=predict.RidgeFusedQDA(best_mod[[j]],newdata = X_test1)
#prd=c(prd,sum(preds1$Predicted==c(rep(1,500),rep(2,500),rep(3,500),rep(4,500))))
#prd_like=c(prd_like, logLike.val(best_mod[[j]],newdata=X_tests,1,response=c(rep(1,500),rep(2,500),rep(3,500),rep(4,500))))
#table(preds1$Predicted,c(rep(1,500),rep(2,500),rep(3,500),rep(4,500)))

newdata=X_test1
K=4
n=dim(newdata)[1]
Result=matrix(nrow=n,ncol=K)
Class=rep(0,n)
for(i in 1:n){
  for(k in 1:K){
    Result[i,k]=-.5*(t(newdata[i,]-Means[[k]])%*%Omeg[[k]]%*%(newdata[i,]-Means[[k]])-determinant(Omeg[[k]])$mod)+log(.25)
  }
  Class[i]=which(Result[i,]==max(Result[i,]))	
  
}
prdT=c(prdT,sum(Class==c(rep(1,500),rep(2,500),rep(3,500),rep(4,500))))
}

save.image(out_name)
