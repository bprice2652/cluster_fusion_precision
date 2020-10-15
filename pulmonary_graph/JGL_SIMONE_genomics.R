### This code implemnts the cooperative lasso and 

library(simone)
library(igraph)
#setwd("~/Dropbox/cluster_inv_cov/genomics_example")
load("gds5499_seed_full.RData")

X1=Reduce("rbind",x_split)
lapply(x_split,dim)
taska=c(rep(1,41),rep(2,30),rep(3,42),rep(4,19))
k1=simone(X1,tasks = as.factor(taska))
k1$networks[[12]]

omegs_outs<-k1$networks[[12]]



library(glasso)

aic_pcen<-function(S,Omega,nc){
  a1=0
  for(m in 1:length(S)){
    zeros=length(which(Omega[[m]]!=0))
    
    a1=a1+nc[m]*(sum(diag(S[[m]]%*%Omega[[m]]))-determinant(Omega[[m]])$mod)+2*zeros
  }
  return(a1)
}

lambda1=10^seq(-2,0,.1)
nc1=sapply(x_split,function(m){dim(m)[1]})
S1=lapply(x_split,function(m){(dim(m)[1]-1)*cov(m)/dim(m)[1]})
library(glasso)
Res_AIC_G=NULL
for(h in 1:length(lambda1)){
  
      OMS=vector("list",4)
      for(j in 1:4){
        OMS[[j]]=glasso(S1[[j]],rho = lambda1[h])$wi
      }
      Res_AIC_G=rbind(Res_AIC_G,c(lambda1[h],aic_pcen(S1,OMS,nc1)))
  
}

OMS=vector("list",4)
for(j in 1:4){
  OMS[[j]]=glasso(S1[[j]],rho = lambda1[8])$wi
}


omegs_outs<-OMS

