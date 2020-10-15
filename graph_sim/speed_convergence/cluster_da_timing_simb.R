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

library(parallel)
library(cluster)
library(RidgeFusion)
library(JGL)
library(compiler)

## A basic common covariance calcuation function
## S is a list of covariance matrices (MLE)
## nc is a vector of sample sizes
MLEcov<-function(S,nc){
  R=0
  n=sum(nc)
  for(i in 1:length(S)){
    R=R+nc[i]*S[[i]]
  }
  return(R/n)
}


## A quick function for a Ridge Calcualtion of a single covariance matrix
## S is a MLE of a covariance matrix
## nc is the sample size
## lam is a tuning parameter for a ridge penalty
RidgeSample<-function(S,nc,lam){
  p=dim(S[[1]])[1]
  L=eigen(MLEcov(S,nc),symmetric=TRUE)
  NewEVS=(-L$val+sqrt(L$val^2+4*lam))/(2*lam)
  NewOmeg=tcrossprod(L$vec*rep(NewEVS,each=p),L$vec)
  return(NewOmeg)
}

#' c_rf_horse
#' @param X is a list of n by p matrices where the jth element is the data for class j
#' @param lamb1 is the ridge tuning parameter and is a number >= 0
#' @param lamb2 is the ridge tuning parameter and is a number >= 0
#' @param nc is a samole size vector of length # of classes
#' @param D is a list of size Q where each element contains the index of classes that are associated with a cluster 
#' @param ncores_horse The parallel operation of the workhorse function 
#' @return omega the estimates of C covariance matrices
#' @author Stephen Banick <stbanick@mix.wvu.edu>, Brad Price

c_rf_horse <- function(X,lamb1,lamb2,nc,D,ncores_horse=1,...){
  
  S=lapply(X,function(z){cov(z)*((dim(z)[1]-1)/dim(z)[1])})
  Q=length(D)
  #om1=vector("list",Q)
  #for(j in 1:Q){
  #  keep=vector("list",length(D[[j]]))
  #  for(k in 1:length(D[[j]])){
  #    keep[[k]]=S[[D[[j]][k]]]
  #  }
  #  om1[[j]]=RidgeFused(keep,lambda1=lamb1,lambda2=lamb2/length(D[[j]]),nc)$Omega
  #}
  
  om1=vector("list",Q)
  keep_me=vector("list",Q)
  for(j in 1:Q){
    keep_me[[j]]$keep=vector("list",length(D[[j]]))
    keep_me[[j]]$size=length(D[[j]])
    keep_me[[j]]$nca=NULL
    for(k in 1:length(D[[j]])){
      keep_me[[j]]$keep[[k]]=S[[D[[j]][k]]]
      keep_me[[j]]$nca=c(keep_me[[j]]$nca,dim(X[[D[[j]][[k]]]])[1])
    }
  }
  
  ##MAKE MCLAPPLY
  om1=mclapply(keep_me,function(h){RidgeFused(h$keep,lambda1=lamb1,lambda2=lamb2/h$size,nc,scale=FALSE)$Omega},mc.cores=min(ncores_horse,length(D)))
  #om1=lapply(keep_me,function(h){RidgeFused(h$keep,lambda1=lamb1,lambda2=lamb2/h$size,nc)$Omega})
  
  omega=vector("list",length(S))
  for(d in 1:Q){
    for(k in 1:length(D[[d]])){
      omega[[D[[d]][k]]]=om1[[d]][[k]]
    }
  }
  return(list(omega=omega))
}

## complied function
c_rf_horse=cmpfun(c_rf_horse)


#' ridge_update
#' @param S an covariance matrix estimate (MLE)
#' @param rho1 a tuning paramter >=0
#' @return omeg an estimate of an inverse covariance matrix
#' @author Brad Price <brad.price@mail.wvu.edu>
ridge_update<-function(S,rho1){
  p=dim(S)[1]
  L2=eigen(S,symmetric=TRUE)
  NewEVS=(-L2$val+sqrt(L2$val^2+4*rho1))/(2*rho1)
  NewOmeg=tcrossprod(L2$vec*rep(NewEVS,each=p),L2$vec)
  return(list(omeg=NewOmeg))
}


### For Brad Local:
#setwd("~/Dropbox/cluster_inv_cov/r_code")
#dyn.load("prox_grad.so")


## Aaron's Prox Gradient Descent Update
## S is the covariance matrix
## alpha is the Ridge tuning parameter
## lambda is the lasso tuning parameter
## t0 is the initali stepsize
## Omega.init is the initlaiziation of Oemga.  If true it initalizes at a ridge solution with tuning alpha
## tol is the tolerance
## max.iter is max iterations
## gamma is the backtrack step
## quiet does not print the objective function.
PGD_EN <- function(S, alpha, lambda, t0 = 1, Omega.init = NULL, tol = 1e-8, max.iter = 1e4, gamma = .5, quiet=TRUE){
  
  
  # -----------------------------------------
  # inner functions 
  # -----------------------------------------
  grad.eval <- function(S, Omega, alpha){
    return(S - chol2inv(chol(Omega)) + 2*alpha*Omega)
  }
  
  grad.eval.chol <- function(S, Omega, Omega.chol, alpha){
    return(S - chol2inv(Omega.chol) + 2*alpha*Omega)
  }
  
  f.eval <- function(S, Omega, alpha){
    return(sum(Omega*(t(S + alpha*Omega))) - determinant(Omega, logarithm=TRUE)$modulus[1])
  }
  
  f.eval.chol <- function(S, Omega, Omega.chol, alpha){
    return(sum(Omega*(t(S + alpha*Omega))) - 2*sum(log(diag(Omega.chol))))
  }
  
  # -----------------------------------------
  # preliminaries and initializers
  # -----------------------------------------
  k.iter <- 1
  p <- dim(S)
  if(is.null(Omega.init)){
    Omega.init = diag(1, p)
  }
  if(isTRUE(Omega.init)){
    
    RidgeSample<-function(S,lam){
      p=dim(S[[1]])[1]
      L=eigen(S,symmetric=TRUE)
      NewEVS=(-L$val+sqrt(L$val^2+4*lam))/(2*lam)
      NewOmeg=tcrossprod(L$vec*rep(NewEVS,each=p),L$vec)
      return(NewOmeg)
    }
    Omega.init=RidgeSample(S,(alpha+lambda))
  }
  
  Omega.k <- Omega.init
  objfunc.val <- rep(0, max.iter)
  old.f.val <- f.eval(S, Omega.init, alpha)
  orig.objfunc.val <- old.f.val + lambda*sum(abs(Omega.init))
  old.objfunc.val <- orig.objfunc.val
  old.grad <- grad.eval(S, Omega.k, alpha)
  converged <- FALSE
  a1=eigen(S)$val[1]
  a=1/(.5*((a1+lambda*dim(S)[1]+sqrt((a1+lambda*dim(S)[1])^2+8*alpha))))
  t0=a^2/((2*a^2*alpha)+1)
  # -------------------------------------
  # iterations
  # -------------------------------------
  for(kk in 1:max.iter){
    
    line.search <- TRUE
    # --- proximal step & line search
    t <- t0
    
    
    while(line.search){
      
      # -- candidate iterate
      Omega.temp <- pmax(abs(Omega.k - t*old.grad) - lambda*t, 0)*sign(Omega.k - t*old.grad)
      
      # -- check PD-ness
      temp <- try(chol(Omega.temp),silent = TRUE)
      
      if(class(temp)!="try-error"){
        
        # -- check majorization condition
        f.val.temp <- f.eval.chol(S, Omega.temp, temp, alpha)
        grad.temp <- grad.eval.chol(S, Omega.temp, temp, alpha)
        RHS <- old.f.val + sum((Omega.temp - Omega.k)*(t(old.grad))) + (1/(2*t))*sum((Omega.temp - Omega.k)^2)
        
        
        if(f.val.temp <= RHS){
          line.search = FALSE
        } else {
          t <- t*gamma
        }
        
      } else {
        
        t <- t*gamma
        if(t < 1e-20){
          line.search <- FALSE
        }
      }
      
    }
    
    Omega.kp1 <- Omega.temp
    new.objfunc.val <- f.val.temp + lambda*sum(abs(Omega.kp1))
    objfunc.val[kk] <- new.objfunc.val
    objfunc.change <- (old.objfunc.val - new.objfunc.val)/(abs(orig.objfunc.val))
    if(objfunc.change < tol){
      converged <- TRUE
      break
    } else {
      
      old.objfunc.val <- new.objfunc.val
      old.grad <- grad.temp
      old.f.val <- f.val.temp
      Omega.k <- Omega.kp1
      
    }
    
    if(!quiet){
      cat(objfunc.val[kk], "\n")
    }
  }
  
  # --- check first order conditions
  suboptim1 <- max(abs((abs(old.grad) - lambda)[Omega.k!=0]))
  suboptim2 <- sum((abs(old.grad)<=lambda)[Omega.k==0])/sum(Omega.k==0)
  
  
  return(list("Omega" = Omega.k, "objfunc" = objfunc.val[1:kk], "subopt1" = suboptim1, "subopt2" = suboptim2))
}




#' pcen
#' Note that this is used for one set of joint estimations and does not consder clusters, much like the RidgeFusion appraoch
#' @param X a list of nc x p data matrices, each element represents one class
#' @param nc a vector of sample sizes corresponding to X
#' @param lambda1 tuning parameter for the ridge penalty
#' @param lambda2 tuning parameter for cluster fusion penalty
#' @param rho ADMM tuning parameter
#' @param tol Coordinate descent convergence tolerance
#' @param tol1 absolute tolerance
#' @param tol2 relative tolerance
#' @param maxiter max iterations for CD alogrithm
#' @param gista_max_iter maximum number of iterations for the ADMM algorithm
#' @param backtrack  the backtrack step
#' @param scale Logical if TRUE shrinkage is done on correlation scale

pcen<-function(X,nc,lambda1,lambda2,rho=10,tol=10^-8,tol1=10^-7,tol2=10^-7,maxiter=100,gista_max_iter=1000, backtrack=.2, scale=FALSE){
  S=lapply(X,function(z){cov(z)*((dim(z)[1]-1)/dim(z)[1])})
  nc=sapply(X,function(x){dim(x)[1]})
  iter=0
  obj_step=NULL
  J=length(S)
  p=dim(S[[1]])[1]
  n=sum(nc)
  # if(is.null(warm.start)){
  Omeg1=lapply(S,function(x){return((diag(1/diag(x),p,p)))})
  #}else{Omeg1=warm.start}
  
  
  getsumabs=function(om){return(sum(abs(diag(om))))}
  tolb=tol*sum(sapply(Omeg1,getsumabs))
  #tolb=tol
  
  
  if(lambda2>=10^4){
    Init=RidgeSample(S,nc,lambda1)
    for(i in 1:J){
      Omeg1[[i]]=Init
    }
  }
  
  
  
  Diffiter=tolb+1
  #Z=0
  #for(m in 1:J){
  #  Z=Z+Omeg1[[m]]
  #}
  D=lapply(S,function(x){return(sqrt(1/diag(x)))})
  if(scale==TRUE){S1=lapply(S,cov2cor)}else{S1=S}	
  iter=0
  W=Reduce("+",Omeg1)
  while(Diffiter>tolb && iter<=maxiter ){
    Diffiter=0
    for(j in 1:J){
      Omeg=Omeg1
      #new_omeg<-admm_update(S1[[j]]-((lambda2)*(W-Omeg[[j]])/nc[j]),nc = nc[j],l1 = lambda1,l2=lambda2,rho=rho,eps_abs = eps_abs,eps_rel = eps_rel,iter_max = admm_max_iter)$omeg
     # new_omeg<-matrix(.C("G_ISTA_C",p=as.integer(p),S=as.double(as.vector(S1[[j]]-((lambda2)*(W-Omeg[[j]])/nc[j]))),rho=as.double(lambda1/nc[j]),rho2=as.double(lambda2*(length(S1)/nc[j])),gap_tol=as.double(tol1),term_tol=as.double(tol2),
    #               backtrack= as.double(backtrack),max_iter=as.integer(gista_max_iter),Omega=as.double(as.vector(diag(1,p))),Omega_inv=as.double(as.vector(diag(1,p))),opt_val=as.double(0),
    #               time_elapsed=as.double(0), cpu_time=as.double(0),total_it=as.integer(0), gap=as.double(0))$Omega,byrow = TRUE,nrow = p,ncol = p)
       check<-PGD_EN(S = S1[[j]]-((lambda2/nc[j])*((W-Omeg[[j]]))),alpha=(lambda2/(2*nc[j]))*(length(S1)-1),lambda=lambda1/nc[j],t0 = 1,max.iter = 1e4,quiet = TRUE,tol = 1e-8)
       #print(check$subopt1)
        new_omeg=check$Omega
       Diffiter=Diffiter+sum(abs(new_omeg-Omeg1[[j]]))
      W=W-Omeg1[[j]]+new_omeg
      Omeg1[[j]]=new_omeg

      #obj_step=c(obj_step,v1)
    }
    iter=iter+1
    #print(sum(sapply(1:length(S1),function(b){nc[b]*(sum(diag(S1[[b]]%*%Omeg1[[b]]))-determinant(Omeg[[b]])$mod)+lambda1*sum(abs(Omeg1[[j]]))+lambda2*sum(Omeg[[1]]^2)})))
   
    }
  
  k1=0
  k2=0
  k3=0
  for(b in 1:length(S1)){
    k1=k1+nc[b]*(sum(diag(S1[[b]]%*%Omeg1[[b]]))-determinant(Omeg1[[b]])$mod)
    k2=k2+lambda1*sum(abs(Omeg1[[b]]))
    k3a=0
    for(h in (b):length(S1)){
      k3a=k3a+sum(abs(Omeg1[[b]]-Omeg1[[h]]))
    }
    k3=k3+k3a
  }
  v1=k1+k2+(lambda2/(2*(length(Omeg1)-1)))*k3
  
  if(scale==TRUE){
    Omeg2=vector("list",J)
    for(k in 1:J){
      #Omeg1[[k]][which(abs(Omeg1[[k]])<=10^-8)]=0
      Omeg2[[k]]=rep(D[[k]],each=p)*Omeg1[[k]]*rep(D[[k]],each=p)
    }
  }else{
    Omeg2=Omeg1
  }
  
  Rest=list(Omega=Omeg2,Ridge=lambda1,FusedRidge=lambda2,iter=iter,obj=v1)
  class(Rest)="RidgeFusion"	
  return(Rest)
}


#' c_pcen_horse 
#' @param X a list of nc x p data matrices corresponding to each class
#' @param lamb1 tuning parameter for lasso penalty
#' @param lamb2 tuning parameter for ridge fusion penalty
#' @param rho tuning paramter for ADMM
#' @param nc a vector of sample sizes corresponding to X
#' @param D a list of corresponding clusters
#' @param ncores_horse number of cores to use if parallel needed/wanted maximum will be length of D. 
#' @return omega a list of inverse covariance estimates
c_pcen_horse <- function(X,lamb1,lamb2,rho,nc,D,ncores_horse=1,...,scale=FALSE){
  
  S=lapply(X,function(z){cov(z)*((dim(z)[1]-1)/dim(z)[1])})
  Q=length(D)
  #om1=vector("list",Q)
  #for(j in 1:Q){
  #  keep=vector("list",length(D[[j]]))
  #  for(k in 1:length(D[[j]])){
  #    keep[[k]]=S[[D[[j]][k]]]
  #  }
  #  om1[[j]]=RidgeFused(keep,lambda1=lamb1,lambda2=lamb2/length(D[[j]]),nc)$Omega
  #}
  
  om1=vector("list",Q)
  keep_me=vector("list",Q)
  for(j in 1:Q){
    keep_me[[j]]$keep=vector("list",length(D[[j]]))
    keep_me[[j]]$size=length(D[[j]])
    keep_me[[j]]$nca=NULL
    for(k in 1:length(D[[j]])){
      keep_me[[j]]$keep[[k]]=X[[D[[j]][k]]]
      keep_me[[j]]$nca=c(keep_me[[j]]$nca,dim(X[[D[[j]][[k]]]])[1])
    }
  }
  
  ##MAKE MCLAPPLY
  
#  om1=vector("list",length(om1))
# for(j in 1:length(keep_me)){
#    om1[[j]]=pcen(X=keep_me[[j]]$keep,lambda1=lamb1,lambda2=lamb2/(keep_me[[j]]$size),nc=keep_me[[j]]$nca)$Omega
#    print(j)
#  }
  om2=mclapply(keep_me,function(h){pcen(X = h$keep,nc=h$nca,lambda1=lamb1,lambda2=lamb2/(h$size),scale=scale)},mc.cores=min(ncores_horse,length(D)))
  obj_vec=Reduce("+",lapply(om2,function(z){z$obj}))
  om1=lapply(om2,function(z){z$Omega})
  #om1=lapply(keep_me,function(h){RidgeFused(h$keep,lambda1=lamb1,lambda2=lamb2/h$size,nc)$Omega})
  
  omega=vector("list",length(S))
  for(d in 1:Q){
    for(k in 1:length(D[[d]])){
      omega[[D[[d]][k]]]=om1[[d]][[k]]
    }
  }
  return(list(omega=omega,obj_fun=obj_vec))
}

## Compiled version of c_pen_horse for testing. 
c_pcen_horse=cmpfun(c_pcen_horse)

#' Ridge_mine
#' @param S covariance matrix
#' @param nc sample size
#' @param lam tuning parameter for non-fusion penalty
#' @return inverse covariance estimate
#' @author Brad Price <brad.price@mail.wvu.edu>, Stephen Banick
#' 
Ridge_mine<-function(S,nc,lam){
  p=dim(S)[1]
  S=S*(nc-1)/(nc)
  L=eigen(S,symmetric=TRUE)
  NewEVS=(-L$val+sqrt(L$val^2+4*lam))/(2*lam)
  NewOmeg=tcrossprod(L$vec*rep(NewEVS,each=p),L$vec)
  return(NewOmeg)
}

#' SetEq
#' @param set1 is a vector of positive integers
#' @param set2 is a vector of positive integers
#' @return is a logical if the two sets are equivalent 
#' @author Stephen Banick <stbanick@mix.wvu.edu>, Brad Price

SetEq <-function(set1,set2){
  s1<-vector("list",max(set1))
  for(i in 1:max(set1)){
    s1[[i]]=which(set1==i)
  }
  Sets=sapply(s1,function(z){
    s2=vector("list",max(set2))
    for(i in 1:max(set2)){
      s2[[i]]=which(set2==i)
    }
    sum(sapply(s2,function(w){identical(z,w)}))
  })
  sum(Sets)==max(set1)
}

#' clus_prep
#' @param theta is a list of matrices
#' @return matrix, where each row is a vectorized form of a corresponding matrix.
#' @author Stephen Banick <stbanick@mix.wvu.edu>, Brad Price
clus_prep<-function(theta){
  sapply(theta,function(x)as.vector(x))
}




#' inv_clust
#' @theta is a list of matrices
#' @Q is the number of clusters
#' @return the cluster assignments
#' @author Stephen Banick <stbanick@mix.wvu.edu>, Brad Price
#'

inv_clust<-function(theta,Q){
  mine<-t(clus_prep(theta))
    D=kmeans(mine,centers = Q,nstart=100)$cluster
  return(list(clust=D))
}

inv_clust<-cmpfun(inv_clust)

#' prep_me
#' @param m an n times p data matrix
#' @return a list of an inverse covariance matrix and the dimension of the matrix. 
prep_me<-function(m){
  return(list(sig=cov(m), n=dim(m)[1]))
  ## THIS WAS solve(cov(m))
}


identicalValue<-function(ests){
  v1=sapply(1:length(ests),function(m){
    sapply(1:length(ests),function(z){
      all.equal(ests[[m]],ests[[z]])
    }
    )})
  sum((v1==TRUE)/length(v1))==1
}



calc_obj_fun<-function(S,inv,clust,lambda2,lambda1){
  Q=max(clust)
  vk=0
  for(g in 1:Q){
  S1=which(clust==g)  
  k1=0
  k2=0
  k3=0
  for(b in S1){
    k1=k1+nc[b]*(sum(diag(S[[b]]%*%inv[[b]]))-determinant(inv[[b]])$mod)
    k2=k2+lambda1*sum(abs(inv[[b]]))
    k3a=0
    if(length(S1)>1){
    for(h in (b):length(S1)){
      k3a=k3a+sum(abs(inv[[b]]-inv[[h]])^2)
    }
    }
    k3=k3+k3a
  }
  if(length(S1)>1){
  v1=k1+k2+(lambda2/(2*(length(S1)-1)))*k3
  }else{
    v1=k1+k2
  }
  print(v1)
  vk=vk+v1
  }
return(vk)
}


#' cluster_inv_cov
#' @param X a list of n times p matrices where each element of a list corresponds to a class
#' @param lambda1 tuning parameter for non-fusion penalty
#' @param lambda2 tuning parameter for fusion penalty
#' @param rho admm parameter for pcen
#' @param Q nubmer of clusters
#' @param method defining of inverse covariance estimation can take values "ridgeFusion" "fgl" or "hybrid"
#' @param ncores the number of cores to be used in the parallelization of the method. 
#' @param max_iter the maximum number of iterations for the large routine
#' @return inverse covariance estimates and estimates of the clusters
#' @author Brad Price <brad.price@mail.wvu.edu>, Stephen Banick
#' Note that this is the workhorse for the estimation of the precision matrix
cluster_inv_cov<-function(X,lambda1,lambda2,rho=10, Q, method="ridgeFusion",ncores=1,max_iter=20,scale=FALSE){
  
  ##Put SCale in here
  S_n=lapply(X,prep_me)
  D=lapply(S_n,function(x){return(sqrt(1/diag(x$sig)))})
  if(scale==TRUE){S_n=lapply(S_n,function(d){sig=cov2cor(d$sig)
  return(list(sig=sig,n=d$n))})}else{S_n=S_n}
  
  
  
  
  
  #inits<-lapply(S_n,function(m){Ridge_mine(m$sig,nc=m$n,lam=lambda1/m$n)})
  if(lambda1<10^3){
    inits=pcen(X,nc = sapply(X,function(m){dim(m)[1]}),lambda1 = lambda1,lambda2 = 0)$Omega}
  else{
    inits=lapply(lapply(S_n,function(z){z$sig}),function(g){diag(1/diag(g),dim(g)[1],dim(g)[1])})
  }
  clust_new<-inv_clust(theta = inits,Q =Q)$clust
  obj_fun=calc_obj_fun(lapply(S_n,function(g){g$sig}),inits,clust_new,lambda2,lambda1)
  
  
  conv=FALSE
  nc=sapply(X,function(m){dim(m)[1]})
  iter=0
  while(conv==FALSE && iter<=max_iter){
    clust_old=clust_new
    clust_list<-lapply(1:Q,function(m){which(clust_old==m)})
    if(method=="ridgeFusion"){
      ests=c_rf_horse(X,lamb1=lambda1,lamb2 =lambda2,nc=nc,D=clust_list,ncores_horse = ncores)$omega
  }
    if(method=="pcen"){
      ests=c_pcen_horse(X,lamb1=lambda1, lamb2=lambda2,rho=rho,nc=nc,D=clust_list,ncores_horse=ncores)
      ests=ests$omega
      obj_fun=c(obj_fun,calc_obj_fun(lapply(S_n,function(g){g$sig}),ests,clust_old,lambda2,lambda1))
    }
    if(isTRUE(identicalValue(ests))){
      conv=TRUE
    }else{
      clust_new<-inv_clust(theta = ests,Q =Q)$clust
      obj_fun=c(obj_fun,calc_obj_fun(lapply(S_n,function(g){g$sig}),ests,clust_new,lambda2,lambda1))
      conv=SetEq(clust_old,clust_new)
      iter=iter+1
      print(iter)
    }
  }
  
  
  if(scale==TRUE){
    J=length(S_n)
    Omeg2=vector("list",J)
    for(k in 1:J){
      #Omeg1[[k]][which(abs(Omeg1[[k]])<=10^-8)]=0
      Omeg2[[k]]=rep(D[[k]],each=p)*ests[[k]]*rep(D[[k]],each=p)
    }
    ests=Omeg2
  }
  
  return(list(omega=ests,clusters=clust_new,obj_fun=obj_fun)) 
}

cluster_inv_cov<-cmpfun(cluster_inv_cov)

#' inv_cov_cv
#' @param X a list of n times p data matrices where each element represents a class
#' @param lambda1 vector of possible tuning parameters for non-fusion penalty
#' @param lambda2 vector of possible tuning parameters for fusion penalty
#' @param Q an integer defining the nubmer of clusters
#' @param rho the admm penalty parameter
#' @param method inverse covariance matrix estimation procedure
#' @param Fold a list of size of the number of folds, where each element is a list of size number of classes containing the indices are contained in that cluster.
#' @param tol Convergence tolerance
#' @param ncores.cv Number of cores for parallel processing.  Can be bounded by number of folds.
#' @return Best tuning parameters and cross validation scores for all parameters.
inv_cov_cv<-function(X,lambda1,lambda2,Q,rho,method="ridgeFusion",Fold,tol=10^-6,ncores_cv=1,ncores_horse=1,...){
  #print(paste("working on inv_cov_cv with lambda1", lambda1, "lambda2", lambda2, "Q", Q))
  #S2=lapply(X,function(x){(dim(x)-1)*cov(x)/dim(x)})
  L1=length(lambda1)
  L2=length(lambda2)
  nc=lapply(X,function(x){dim(x)[1]})
  J=length(X)
  nc=lapply(X,function(x){dim(x)[1]})
  p=dim(X[[1]])[2]
  if(length(Fold)==1){stop("Need more than 1 fold")}
  
  if(length(Fold)>1){
    K=length(Fold)
    XValid=list(0,0)
    XTest=list(0,0)
    S=XValid
    TV2=XValid	
    
    for(k in 1:K){
      XValid[[k]]=list(0,0)
      XTest[[k]]=list(0,0)
      #S[[k]]=list(0,0)
      for(j in 1:J){
        XValid[[k]][[j]]=X[[j]][-Fold[[k]][[j]],]
        XTest[[k]][[j]]=X[[j]][Fold[[k]][[j]],]
      }
      # S[[k]]=lapply(XValid[[k]],function(x){return(list((dim(x)[1]-1)*cov(x)/(dim(x)[1]),dim(x)[1]))})
      TV2[[k]]=lapply(XTest[[k]],function(x){(dim(x)[1]-1)*cov(x)/(dim(x)[1])})
    }
    
    #Inits=vector("list",L1*L2)
    #for(i in 1:L1){
    #  for(j in 1:L2){
    #    Inits[[L2*(i-1)+j]]=0
    #    Inits[[L2*(i-1)+j]]=mclapply(XValid,function(m){cluster_inv_cov(m,lambda1 = lambda1[j],lambda2=max(lambda2), rho=10, Q=Q,method = method,clus_method = clust_method)$omeg},mc.cores = min(ncores_cv,K))
    
    
    #cluster_inv_cov(XValid[[k]],lambda1 = lambda1,lambda2=lambda2,rho = rho, Q=Q)
    #RidgeFusedL(S[[k]],lambda1[i],Inf,##lambda2[j],
    #tole=tol,ws=NULL,scaleL=scaleCV)
    #  }
    #}
    Rest=matrix(0,length(lambda1),length(lambda2))
    for(m in 1:length(lambda1)){
      OmegEst=0
      for(l in 1:length(lambda2)){
        M=0
        OmegEst=list(0,0)		
        OmegEst2=OmegEst
        OmegEst=mclapply(XValid,function(h){cluster_inv_cov(h,lambda1 = lambda1[m],lambda2 = lambda2[l],rho = rho,Q=Q,method = method,ncores=ncores_horse,...)$omeg},mc.cores=min(K,ncores_cv))
        
        for(k in 1:K){		
          #if(warm.start==TRUE){
          #OmegEst[[k]]=lapply()
          #RidgeFusedL(S[[k]],lambda1[m],lambda2[l],tole=tol,ws=Inits[[L2*(m-1)+l]],
          #scaleL=scaleCV)
          #    }else{
          #RidgeFusedL(S[[k]],lambda1[m],lambda2[l],tole=tol,ws=NULL,scaleL=scaleCV)
          for(j in 1:J){
            M=M+(length(Fold[[k]][[j]])/2)*(sum(diag(TV2[[k]][[j]]%*%OmegEst[[k]][[j]]))-determinant(OmegEst[[k]][[j]])$mod)
          }
        }
        Rest[m,l]=M/K
      }
    }
    #if(INF==TRUE){
    #  InfRest=rep(0,length(lambda1))
    #  for(m in 1:length(lambda1)){
    #    M=0
    #    InfEst=mapply(RidgeSampleL,S,lambda1[m],SIMPLIFY=FALSE)
    #    for(k in 1:k){
    #      for(j in 1:J){
    #        M=M+(length(Fold[[k]][[j]])/2)*(sum(diag(TV2[[k]][[j]]%*%InfEst[[k]]))-log(det(InfEst[[k]])))
    #     }
    #    }
    #   InfRest[m]=M/K
    # }
    # Rest=cbind(Rest,InfRest)
    #  lambda2=c(lambda2,Inf)
    #}
    colnames(Rest)=lambda2
    row.names(Rest)=lambda1
    Min=which(Rest==min(Rest),arr.ind=TRUE)[1,]
    Lam1=lambda1[Min[1]]
    Lam2=lambda2[Min[2]]
  }
  
  Ret=list(best_p1=Lam1,best_fused=Lam2,CV=Rest)
  return(Ret)
}

inv_cov_cv<-cmpfun(inv_cov_cv)

#' clust_inv_cov_cv
#' @param X a list of n times p data matrices where each element represents a class
#' @param lambda1 vector of possible tuning parameters for non-fusion penalty
#' @param lambda2 vector of possible tuning parameters for fusion penalty
#' @param Q_set a vector containing possible number of clusters
#' @param rho the admm penalty parameter
#' @param method inverse covariance matrix estimation procedure
#' @param Fold a list of size of the number of folds, where each element is a list of size number of classes containing the indices are contained in that cluster.
#' @param tol Convergence tolerance
#' @param ncores_Q Number of cores for parallel processing based on Q.  Can be bounded by maximum number of clusters.
#' @param ncores_cv Number of cores used in parallel processing for cross validtion.  Can be bounded by number of folds.
#' @return Best tuning parameters for eaching possible cluster, and the best overall. 
clust_inv_cov_cv<-function(X,lambda1,lambda2,Q_set,rho,method="ridgeFusion", Fold,tol=10^-6,ncores_Q=1,ncores_cv=1,ncores_horse=1){
  
  keep_me=mclapply(Q_set,function(h){
    inv_cov_cv(X,lambda1 = lambda1, lambda2=lambda2, Q = h, rho=rho,method=method,Fold=Fold,tol=tol,ncores_cv=ncores_cv,ncores_horse=ncores_horse)
  },mc.cores = min(ncores_Q,length(Q_set)))
  # keep_me=vector("list",length(Q_set))
  #for(j in 1:length(Q_set)){
  #  keep_me[[j]]=inv_cov_cv(X,lambda1=lambda1, lambda2=lambda2, Q=Q_set[j],rho=rho,method=method,clust_method=clust_method, Fold=Fold,ncores=ncores_horse)
  #  }
  
  
  get_res<-sapply(keep_me,function(v){
    
    keep_n=which(v$CV==min(v$CV),arr.ind = TRUE)[1,]
    c(as.numeric(rownames(v$CV)[keep_n[1]]),as.numeric(colnames(v$CV)[keep_n[2]]),v$CV[keep_n[1],keep_n[2]])})
  # a1=which(as.numeric(rownames(v$CV))==v$best_p)[1]
  #b1=which(as.numeric(colnames(v$CV))==v$best_f)[1]
  #c(a1,b1,v$best_p1,v$best_f,v$CV[a1,b1])})  
  get_res=cbind(Q_set,t(get_res))
  colnames(get_res)=c("Clusters","lambda1","lambda2","VL")
  keep_out<-which(get_res[,4]==min(get_res[,4]))[1]
  out_put=list(top_in_cluster=get_res,best_triplet=get_res[keep_out,-4])#,everything=keep_me)
  return(out_put)
}

clust_inv_cov_cv<-cmpfun(clust_inv_cov_cv)


#' clust_qda
#' @param X is a list of n by p matrices where the jth element is the data for class j
#' @param lambda1 is the ridge tuning parameter and is a number >= 0
#' @param lambda2 is the ridge tuning parameter and is a number >= 0
#' @param Q is the number of clusters
#' @param rho is a tuning parameter only applicable for method hybrid
#' @param method takes characters ridgeFusion, fgl, and hybrid
#' @param clus_method either PAM or kmeans
#' @return omega the estimates of C covariance matrices
#' @author Stephen Banick <stbanick@mix.wvu.edu>, Brad Price
#' 
clust_qda<-function(X,lambda1,lambda2,rho=10, Q, method="ridgeFusion"){
  
  Means=lapply(X,function(x){apply(x,2,mean)})
  SampCov=lapply(X,function(x){(dim(x)-1)[1]*cov(x)/(dim(x)[1])})
  nc=sapply(X,function(x){dim(x)[1]})
  pi=unlist(nc)/(sum(unlist(nc)))
  Precision=cluster_inv_cov(X,lambda1=lambda1,lambda2=lambda2, Q=Q, method=method)
  clusts=Precision$clust
  Ret=list(Omeg=Precision$omeg, Means=Means, Pi=pi,lambda1=lambda1,lambda2=lambda2,clusters=clusts)
  class(Ret)="RidgeFusedQDA"
  return(Ret)	
}


## just brining in the RidgeFused QDA prediction this won't go in the software. 
predict.RidgeFusedQDA<-function(object,newdata,class=TRUE,...){
  x=object
  K=length(x$Omeg)
  n=dim(newdata)[1]
  Result=matrix(nrow=n,ncol=K)
  Class=rep(0,n)
  for(i in 1:n){
    for(k in 1:K){
      Result[i,k]=-.5*(t(newdata[i,]-x$Means[[k]])%*%x$Omeg[[k]]%*%(newdata[i,]-x$Means[[k]])-determinant(x$Omeg[[k]])$mod)+log(x$Pi[[k]])
    }
    Class[i]=which(Result[i,]==max(Result[i,]))	
  }
  
  if(class==TRUE){
    return(list(Predicted=Class))
  }
  if(class==FALSE){
    return(list(Predicted=Result))
  }
}



#' aic_pen
#' @param S is a list of p by p covariance matrices matrices where the jth element is the data for class j
#' @param Omega is a list of p by p inverse covariance matrix estimates, where each must be positive definite.  The jth element is the estimate of class j.
#' @param nc is a vector of sample sizes associated with class.  
#' @return omega the estimates of C covariance matrices
#' @author Brad Price <brad.price@mail.wvu.edu>
#' 
aic_pcen<-function(S,Omega,nc){
  a1=0
  for(m in 1:length(S)){
    zeros=length(which(Omega[[m]]!=0))
    a1=a1+nc[m]*(sum(diag(S[[m]]%*%Omega[[m]]))-determinant(Omega[[m]])$mod)+(2*zeros)
  }
  return(a1)
}

bic_pcen<-function(S,Omega,nc){
  a1=0
  for(m in 1:length(S)){
    zeros=length(which(Omega[[m]]!=0))
    a1=a1+nc[m]*(sum(diag(S[[m]]%*%Omega[[m]]))-determinant(Omega[[m]])$mod)+((log(nc[m])/nc[m])*zeros)
  }
  return(a1)
}


clust_inv_cov_aic<-function(X,lambda1,lambda2,Q_set,rho,tol=10^-6,ncores_Q=1,ncores_cv=1,ncores_horse=1,returnAll=FALSE){
  S=lapply(X,function(z){cov(z)*((dim(z)[1]-1)/dim(z)[1])})
  nc=sapply(X,function(z){dim(z)[1]})
  keeps=expand.grid(lambda1,lambda2,Q_set)
  colnames(keeps)=c("lambda1","lambda2","Clusters")
  mes=mclapply(1:dim(keeps)[1],function(m){cluster_inv_cov(X ,lambda1=keeps[m,1],lambda2=keeps[m,2],Q=keeps[m,3],method="pcen",ncores = ncores_Q)},mc.cores=ncores_cv)
  aics=unlist(mclapply(1:dim(keeps)[1],function(m){aic_pcen(S = S,nc=nc,Omega = mes[[m]]$omega)}))
  best=which(aics==min(aics))
  aic_grid=cbind(keeps,aics)
  if(returnAll==FALSE){
  return(list(best_pcen=mes[[best]],best_tune=aic_grid[best,],aic_grid=aic_grid))
  }else{
    return(list(grids=aic_grid,ests=mes))
  }
}

compare_test<-function(S,Omega,nc){
  a1=0
  for(m in 1:length(S)){
    zeros=length(which(Omega[[m]]!=0))
    a1=a1+nc[m]*(sum(diag(S[[m]]%*%Omega[[m]]))-determinant(Omega[[m]])$mod)
  }
  return(a1)
}


library(RidgeFusion)
library(JGL)
library(mvtnorm)
library(parallel)
library(LASICH)
library(simone)

setwd("/scratch/bprice5/timing_sim")
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
clusters_list=vector("list",25)
obj_fun_l=vector("list",25)
times=NULL

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
    
    outs2=mclapply(1:dim(gridme)[1],function(b){k1<-system.time(b1<-cluster_inv_cov(X1,lambda1 = gridme[b,1],method="pcen",lambda2=gridme[b,2],Q=2,scale=FALSE))
    return(list(run_time=k1[3],output=b1))},mc.cores =20)
    times1=sapply(outs2,function(g){g$run_time})
    times=rbind(times,times1)
    obj_fun_l[[h]]=lapply(outs2,function(m){m$output$obj_fun})
    clusters_list[[h]]=lapply(outs2,function(m){m$output$clusters})
  }
  
  if(method=="PCEN-3"){
    outs2=mclapply(1:dim(gridme)[1],function(b){k1<-system.time(b1<-cluster_inv_cov(X1,lambda1 = gridme[b,1],method="pcen",lambda2=gridme[b,2],Q=3,scale=FALSE))
    return(list(run_time=k1[3],output=b1))},mc.cores =20)
    times1=sapply(outs2,function(g){g$run_time})
    times=rbind(times,times1)
    obj_fun_l[[h]]=lapply(outs2,function(m){m$output$obj_fun})
    clusters_list[[h]]=lapply(outs2,function(m){m$output$clusters})
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
      for(j in 1:4){
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
    if(setting==4){
      taska=c(rep(1,200),rep(2,200),rep(3,200),rep(4,200),rep(5,200),rep(6,200))
    }else{
      taska=c(rep(1,200),rep(2,200),rep(3,200),rep(4,200)) 
    }
    Xg=Reduce("rbind",X1)
    outs_simone=mclapply(1:dim(gridme)[1],function(b){simone(X=Xg,tasks = as.factor(taska),control=setOptions(penalties=gridme[b,1]))$networks[[1]]},mc.cores=20)
    l2=sapply(outs_simone,function(z){
      keep=0
      for(j in 1:C){
        keep=sum((z[[j]]-true01[[j]])^2)+keep
      }
      return(keep)
    })
    
    l2_tp_fp=sapply(outs_simone,function(z){
      tps=0
      fps=0
      for(j in 1:C){
        tps=tps+sum(which(z[[j]]!=0)%in%which(true01[[j]]!=0))
        fps=fps+sum(which(z[[j]]!=0)%in%which(true01[[j]]==0))
      }
      return(list(tps=tps,fps=fps))
    })
    
  }

}

p1=paste(method,p,l1,setting,"timing_b.RData",sep="")
save.image(p1)
