#!/usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-d", "--pvar"), type="numeric", default=1,
              help="number of variables", metavar="number"),
  make_option(c("-l", "--l1"),type="numeric",default=1,
              help="l1", metavar="number"),
  make_option(c("-s", "--setting"),type="numeric",default=1,
              help="setting_number", metavar="number"),
  make_option(c("-m","--method"),default="JGL",type="character",metavar = "character")
);




opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#### The above code reads in the parameters from the array job on the HPC system.
### p is the dimension of the problem, i.e. number of variables.
### l1 is the fixed value for lambda_2  in our simulation we use 10^c(-3,1,3)
### method is the name of the method, which is PCEN-2, PCEN-3, JGL, GLASSO, LASICH1,
### LASICH2, and SIMONE

### Setting is settings 1-4
## This file runs replications 26-50

p=opt$pvar
l1=opt$l1
method=opt$method
setting=opt$setting
clusters_est=NULL

##read in p
## read in l1
## read in method

library(RidgeFusion)
library(JGL)
library(mvtnorm)
library(parallel)
library(LASICH)
library(simone)
source("cluster_da.R")
#source("/scratch/bprice5/cluster_da.R")
#setwd("/scratch/bprice5/graph_sim")
ar1<-function(p,rho){
  times=1:p
  H=abs(outer(times,times,"-"))
  mat=rho^H
  return(mat)
}



erg<-function(p,con,diff){
  g1<-erdos.renyi.game(p,con,"gnm")
  g1a<-as.matrix(as_adj(g1,type="Upper"))
  k1=which(g1a!=0)
  k1s=which(g1a==0)
  too=sample(k1s,diff,replace=TRUE)
  change=sample(1:length(k1),diff)
  g1ab=g1a
  g1ab[k1[change]]=0
  #g1ab[too]=1
  
  keep=sapply(1:length(k1),function(b){
    y=rbinom(1,size = 1,prob = .5)
    return(y*runif(1, min=-.7,max=-.5)+(1-y)*runif(1,.5,.7))
  })
  keep2=sapply(1:length(too),function(b){
    y=rbinom(1,size = 1,prob = .5)
    return(y*runif(1, min=-.7,max=-.5)+(1-y)*runif(1,.5,.7))
  })
  
  g2a=g1a
  g2a[k1]=keep
  kz=apply(abs(g2a),1,sum)
  kz[which(kz==0)]=1
  g2a=apply(g2a,2,function(z){z/(1.5*kz)})
  diag(g2a)=1
  g2a=(g2a+t(g2a))/2
  #diag(g2a)=apply(abs(g2a),1,sum)+.1
  
  k2=c(k1[-change])
  keep3=c(keep[-change])
  g2ab=g1ab
  g2ab[k2]=keep3+c(runif(length(keep3),-.01,.01))
  kz2=apply(abs(g2ab),1,function(z){(sum(z))})
  kz2[which(kz2==0)]=1
  g2ab=apply(g2ab,2,function(z){z/(1.5*kz2)})
  diag(g2ab)=1
  g2ab=(t(g2ab)+g2ab)/2
  diag(g2ab)=1
  #diag(g2ab)=apply(abs(g2ab),1,sum)+.1
  g2a=cov2cor(solve(g2a))
  g2ab=cov2cor(solve(g2ab))
  diag(g2a)=1
  diag(g2ab)=1
  g2a=round(solve(g2a),10)
  g2ab=round(solve(g2ab),10)
  
  
  return(list(p1=g2a,p2=g2ab))
}

erg_block<-function(p,con,diff){
  m1=diag(1,p,p)
  m2=diag(1,p,p)
  b1<-erg(p/2,con,diff)
  b2<-erg(p/2,con,diff)
  m1[(1:(p/2)),(1:(p/2))]=b1$p1
  m2[(1:(p/2)),(1:(p/2))]=b1$p2
  m1[(((p/2)+1):p),(((p/2)+1):p)]=b2$p1
  m2[(((p/2)+1):p),(((p/2)+1):p)]=b2$p2
  
  return(list(p1=m1,p2=m2))
}



erg_block_set<-function(p,con,diff){
  m1=diag(1,p,p)
  m2=diag(1,p,p)
  set=sample(1:p,replace=FALSE)
  b1<-erg(p/2,con,diff)
  b2<-erg(p/2,con,diff)
  m1[set[(1:(p/2))],set[(1:(p/2))]]=b1$p1
  m2[set[(1:(p/2))],set[(1:(p/2))]]=b1$p2
  m1[set[(((p/2)+1):p)],set[(((p/2)+1):p)]]=b2$p1
  m2[set[(((p/2)+1):p)],set[(((p/2)+1):p)]]=b2$p2
  
  return(list(p1=m1,p2=m2))
}


erg_block_upper<-function(p,con,diff){
  m1=diag(1,p,p)
  m2=diag(1,p,p)
  b1<-erg(p/2,con,diff)
  m1[(1:(p/2)),(1:(p/2))]=b1$p1
  m2[(1:(p/2)),(1:(p/2))]=b1$p2
  return(list(p1=m1,p2=m2))
  
}

erg_block_lower<-function(p,con,diff){
  m1=diag(1,p,p)
  m2=diag(1,p,p)
  b1<-erg(p/2,con,diff)
  m1[(((p/2)+1):p),(((p/2)+1):p)]=b1$p1
  m2[(((p/2)+1):p),(((p/2)+1):p)]=b1$p2
  return(list(p1=m1,p2=m2))
  
}

cov_gen<-function(num,p,rho){
  if(num==1){
    return(erg_block(p,p,rho))
  }
  if(num==2){
    return(erg(p,floor(p*(p-1)*.20),rho))
  }
  if(num==3){
    return(erg_block_set(p,p,rho))
  }
  if(num==4){
    return(erg_block_upper(p,p,rho))
  }
  if(num==5){
    return(erg_block_lower(p,p,rho))
  }
  
}

rho1=.9
rho2=.91
res_mat=NULL
tp_mat=NULL
fp_mat=NULL

for(h in 26:50){
set.seed(h+120)
if(setting==1){
  rho3=4
  sb=cov_gen(1,p,rho3)
  s1=solve(sb$p1)
  s2=solve(sb$p2)
  sa=cov_gen(3,p,rho3)
  s3=solve(sa$p1)
  s4=solve(sa$p2)
  
  
  ss=list(s1,s2,s3,s4)  
}
if(setting==2){
  rho3=.2*p
  sb=cov_gen(1,p,rho3)
  s1=solve(sb$p1)
  s2=solve(sb$p2)
  sa=cov_gen(1,p,rho3)
  s3=solve(sa$p1)
  s4=solve(sa$p2)
  ss=list(s1,s2,s3,s4)  
}

if(setting==3){
  rho3=4
  sb=cov_gen(4,p,rho3)
  s1=solve(sb$p1)
  s2=solve(sb$p2)
  sa=cov_gen(5,p,rho3)
  s3=solve(sa$p1)
  s4=solve(sa$p2)
  ss=list(s1,s2,s3,s4)  
  
}
  
if(setting==4){
  rho3=4
  sb=cov_gen(4,p,rho3)
  s1=solve(sb$p1)
  s2=solve(sb$p2)
  sa=cov_gen(5,p,rho3)
  s3=solve(sa$p1)
  s4=solve(sa$p2)
  rho3=.2*p
  sb=cov_gen(1,p,rho3)
  s5=solve(sb$p1)
  s6=solve(sb$p1)
  ss=list(s1,s2,s3,s4,s5,s6)  
  
}

true01=lapply(ss,function(b){round(solve(b),10)})
x1<-rmvnorm(200,rep(0,p),s1)
x2<-rmvnorm(200,rep(0,p),s2)
x3<-rmvnorm(200,rep(0,p),s3)
x4<-rmvnorm(200,rep(0,p),s4)
X1=list(x1,x2,x3,x4)
nc=rep(200,4)
C=4
if(setting==4){
  x5<-rmvnorm(200,rep(0,p),s5)
  x6<-rmvnorm(200,rep(0,p),s6)
  X1=list(x1,x2,x3,x4,x5,x6)
  nc=rep(200,6)
  C=6
}

g1=10^seq(-5,3,.2)
gridme=expand.grid(g1,l1)
colnames(gridme)=c("Lambda1","Lambda2")
if(method=="JGL"){
outs=mclapply(1:dim(gridme)[1],function(b){JGL(X1,lambda1=gridme[b,1],lambda2=gridme[b,2],return.whole.theta = TRUE)},mc.cores = 20)
omegs=lapply(outs,function(z){z$theta})
l2=sapply(omegs,function(z){
  keep=0
  for(j in 1:C){
    keep=sum((z[[j]]-true01[[j]])^2)+keep
  }
  return(keep)
})

l2_tp_fp=sapply(omegs,function(z){
  tps=0
  fps=0
  for(j in 1:C){
    tps=tps+sum(which(z[[j]]!=0)%in%which(true01[[j]]!=0))
    fps=fps+sum(which(z[[j]]!=0)%in%which(true01[[j]]==0))
  }
  return(list(tps=tps,fps=fps))
})
}

if(method=="PCEN-2"){
  
outs2=mclapply(1:dim(gridme)[1],function(b){cluster_inv_cov(X1,lambda1 = gridme[b,1],method="pcen",lambda2=gridme[b,2],Q=2,scale=FALSE)},mc.cores = 20)
k1=lapply(outs2,function(z){z$clusters})
clusters_est=rbind(cbind(h,Reduce("rbind",k1)),clusters_est)
omegs=lapply(outs2,function(z){z$omega})
l2=sapply(omegs,function(z){
  keep=0
  for(j in 1:C){
    keep=sum((z[[j]]-true01[[j]])^2)+keep
  }
  return(keep)
})

l2_tp_fp=sapply(omegs,function(z){
  tps=0
  fps=0
  for(j in 1:C){
    tps=tps+sum(which(z[[j]]!=0)%in%which(true01[[j]]!=0))
    fps=fps+sum(which(z[[j]]!=0)%in%which(true01[[j]]==0))
  }
  return(list(tps=tps,fps=fps))
})
}

if(method=="PCEN-3"){
outs2=mclapply(1:dim(gridme)[1],function(b){cluster_inv_cov(X1,lambda1 = gridme[b,1],method="pcen",lambda2=gridme[b,2],Q=3,scale=FALSE)},mc.cores =20)
k1=lapply(outs2,function(z){z$clusters})
clusters_est=rbind(cbind(h,Reduce("rbind",k1)),clusters_est)
omegs=lapply(outs2,function(z){z$omega})
l2=sapply(omegs,function(z){
  keep=0
  for(j in 1:C){
    keep=sum((z[[j]]-true01[[j]])^2)+keep
  }
  return(keep)
})

l2_tp_fp=sapply(omegs,function(z){
  tps=0
  fps=0
  for(j in 1:C){
    tps=tps+sum(which(z[[j]]!=0)%in%which(true01[[j]]!=0))
    fps=fps+sum(which(z[[j]]!=0)%in%which(true01[[j]]==0))
  }
  return(list(tps=tps,fps=fps))
})
}



if(method=="LASICH-1"){
S1s=lapply(X1,cov)

Sigs=array(0,dim = c(p,p,C))
t0=array(0,dim=c(p,p,C))
true01a=array(0,dim=c(p,p,C))
for(j in 1:C){
  Sigs[,,j]=S1s[[j]]
  t0[,,j]=ss[[j]]
  true01a[,,j]=true01[[j]]
}

k2=matrix(0,C,C)

k2=diag(1,C)
k2[1,2]=-1
k2[1,3]=-1/1000
k2[2,4]=-1/1000
k2[3,4]=-1
if(setting==4){
  k2[5,6]=-1
  k2[1,5]=-1/1000
  k2[3,5]=-1/1000
  k2[2,6]=-1/1000
  k2[4,6]=-1/1000
}

k2=k2+t(k2)
library(LASICH)
outs3=mclapply(1:dim(gridme)[1],function(b){lasich(Sigmas=Sigs,Sigma0=t0,Omega0=true01a,ns=rep(200,C),L=k2,lambda1=gridme[b,1],lambda2=gridme[b,2],tol=10^-5,rho=100,thr=FALSE)},mc.cores = 20)

l_ests=lapply(1:length(outs3),function(m){ 
  lapply(1:C, function(d){outs3[[m]]$B[,,d]})
  })


l2=sapply(1:length(l_ests),function(g){
keep=0
for(j in 1:C){
  keep=sum((l_ests[[g]][[j]]-true01[[j]])^2)+keep
}
keep
})


l2_tp_fp=sapply(1:length(l_ests),function(z){
  tps=0
  fps=0
  for(j in 1:C){
    tps=tps+sum(which(l_ests[[z]][[j]]!=0)%in%which(true01[[j]]!=0))
    fps=fps+sum(which(l_ests[[z]][[j]]!=0)%in%which(true01[[j]]==0))
  }
  return(list(tps=tps,fps=fps))
})
}

if(method=="LASICH-2"){

    
    S1s=lapply(X1,cov)

    Sigs=array(0,dim = c(p,p,C))
    t0=array(0,dim=c(p,p,C))
    true01a=array(0,dim=c(p,p,C))
    for(j in 1:C){
      Sigs[,,j]=S1s[[j]]
      t0[,,j]=ss[[j]]
      true01a[,,j]=true01[[j]]
    }

    k2=matrix(0,C,C)
 
k2=matrix(0,C,C)

k2=diag(1,C)
k2[1,2]=-1
k2[1,3]=-1
k2[2,4]=-1
k2[3,4]=-1
if(setting==4){
  k2[5,6]=-1
  k2[1,5]=-1
  k2[3,5]=-1
  k2[2,6]=-1
  k2[4,6]=-1
}

k2=k2+t(k2)
diag(k2)=3

outs4=mclapply(1:dim(gridme)[1],function(b){lasich(Sigmas=Sigs,Sigma0=t0,Omega0=true01a,ns=rep(200,C),L=k2,lambda1=gridme[b,1],lambda2=gridme[b,2],tol=10^-5,rho=100,thr=FALSE)},mc.cores = 20)

l_ests2=lapply(1:length(outs4),function(m){ 
  lapply(1:C, function(d){outs4[[m]]$B[,,d]})
})


l2=sapply(1:length(l_ests2),function(g){
  keep=0
  for(j in 1:C){
    keep=sum((l_ests2[[g]][[j]]-true01[[j]])^2)+keep
  }
  keep
})

  
  l2_tp_fp=sapply(1:length(l_ests2),function(z){
    tps=0
    fps=0
    for(j in 1:C){
      tps=tps+sum(which(l_ests2[[z]][[j]]!=0)%in%which(true01[[j]]!=0))
      fps=fps+sum(which(l_ests2[[z]][[j]]!=0)%in%which(true01[[j]]==0))
    }
    return(list(tps=tps,fps=fps))
  })
}
 
 if(method=="GLASSO"){
 outs_glasso=mclapply(1:dim(gridme)[1],function(b){JGL(X1,lambda1=gridme[b,1],lambda2=0,return.whole.theta = TRUE)},mc.cores = 20)
  omegs_glasso=lapply(outs_glasso,function(z){z$theta}) 
  l2=sapply(omegs_glasso,function(z){
    keep=0
    for(j in 1:C){
      keep=sum((z[[j]]-true01[[j]])^2)+keep
    }
    return(keep)
  })
  
  l2_tp_fp=sapply(omegs_glasso,function(z){
    tps=0
    fps=0
    for(j in 1:C){
      tps=tps+sum(which(z[[j]]!=0)%in%which(true01[[j]]!=0))
      fps=fps+sum(which(z[[j]]!=0)%in%which(true01[[j]]==0))
    }
    return(list(tps=tps,fps=fps))
  })
 }
if(method=="SIMONE"){
  g1=10^seq(-10,3,.25)
  gridme=expand.grid(g1,l1)
  colnames(gridme)=c("Lambda1","Lambda2")
  if(setting==4){
    taska=c(rep(1,200),rep(2,200),rep(3,200),rep(4,200),rep(5,200),rep(6,200))
  }else{
    taska=c(rep(1,200),rep(2,200),rep(3,200),rep(4,200)) 
  }
  Xg=Reduce("rbind",X1)
  #outs_simone=mclapply(1:dim(gridme)[1],function(b){simone(X=Xg,tasks = as.factor(taska),control=setOptions(penalties=gridme[b,1],n.penalties=1))$networks[[1]]},mc.cores=20)
  v3=setOptions(normalize=FALSE,verbose=FALSE, penalties=sort(gridme[,1],decreasing = TRUE))
  outs_simone=simone(X=Xg,tasks=as.factor(taska),control=v3)$networks
  
  l2=sapply(outs_simone,function(z){
    keep=0
    for(j in 1:C){
      keep=sum((z[[j]]-true01[[j]])^2)+keep
    }
    return(keep)
  })
  if(length(l2)<length(gridme[,1])){
    l2=c(l2,rep(l2[length(l2)],length(gridme[,1])-length(l2)))
  }
  
  l2_tp_fp=sapply(outs_simone,function(z){
    tps=0
    fps=0
    for(j in 1:C){
      tps=tps+sum(which(z[[j]]!=0)%in%which(true01[[j]]!=0))
      fps=fps+sum(which(z[[j]]!=0)%in%which(true01[[j]]==0))
    }
    return(list(tps=tps,fps=fps))
  })
  if(length(l2_tp_fp[1,])<length(gridme[,1])){
    l2_tp_fp=cbind(l2_tp_fp,rbind(rep(l2_tp_fp[1,dim(l2_tp_fp)[2]],length(gridme[,1])-dim(l2_tp_fp)[2]),rep(l2_tp_fp[2,dim(l2_tp_fp)[2]],length(gridme[,1])-dim(l2_tp_fp)[2])))
  }
}
res_mat=rbind(res_mat,c(method,p,l1,h,l2))
tp_mat=rbind(tp_mat,c(method,p,l1,h,l2_tp_fp[1,]))
fp_mat=rbind(fp_mat,c(method,p,l1,h,l2_tp_fp[2,]))
print(h)
}

p1=paste(method,p,l1,setting,"b.RData",sep="")
save.image(p1)

