library(parallel)
library(cluster)
library(RidgeFusion)
library(JGL)
library(compiler)
library(klaR)

###################
## Need to make sure you source the cluster_da file from the first folder of this system.
##########
source("cluster_da.R")
#source("~/Dropbox/cluster_inv_cov/r_code/cluster_da.R")






file=read.csv(url("https://archive.ics.uci.edu/ml/machine-learning-databases/libras/movement_libras.data"),header=FALSE)
news=split(file[,-91],file$V91)
news2=lapply(news,function(m){m[-(1:4),]})
news_test=lapply(news,function(m){m[1:4,]})
news_test_1=Reduce("rbind",news_test)
#table(out<-cluster_inv_cov(news2,10^-7,10^-3,rho=10,Q=2,method="ridgeFusion")$cluster)

fold=vector("list",5)
set.seed(2020)
samp=sample(1:20)
for(k in 1:5){
  fold[[k]]=lapply(vector("list",15),function(m){samp[((k-1)*4+1):(k*4)]})
}

system.time(a3<-clust_inv_cov_cv(news2,lambda1 = 10^seq(-7,1,.25),lambda2 = 10^seq(-7,1,.25),rho=10,Q_set=2:5,method="ridgeFusion",Fold=fold,ncores_Q = 4,ncores_cv = 5,ncores_horse = 2))
save.image("libras_run_2a.Rdata")
#load(file = "libras_run.Rdata")
a3_out=a3$best_triplet
#cluster_inv_cov(news2,10^-7,10^-3,rho=10,Q=2,method = "ridgeFusion",ncores = 10)
a_now=clust_qda(news2,a3_out[2],a3_out[3],rho=10,Q=a3_out[1],method="ridgeFusion")
outs<-predict(a_now,newdata=as.matrix(news_test_1,nrow=60,ncol=90))$Predicted
mine<-as.vector(sapply(1:15,function(m){rep(m,4)}))
sum(mine==outs)
## Result is 47
## Classifies 9 of the 15 classes perfectly based on this estimator.  


Abnow<-RidgeFusedCV(news2, 10^seq(-7,2,.25),10^seq(-7,2,.25),Fold=fold,warm.start = TRUE, scale=TRUE)
a_nowb=FusedQDA(news2,Abnow$BestRidge,Abnow$BestFusedRidge,scaleC=TRUE)
outsb<-predict(a_nowb,newdata=as.matrix(news_test_1,nrow=60,ncol=90))$Predicted
mine<-as.vector(sapply(1:15,function(m){rep(m,4)}))
sum(mine==outsb)

Abnow2<-RidgeFusedCV(news2,10^seq(-7,2,.25),0, Fold=fold, warm.start=TRUE,scale=TRUE)
a_nowb=FusedQDA(news2,Abnow2$BestRidge,0,scale=TRUE)
outsb<-predict(a_nowb,newdata=as.matrix(news_test_1,nrow=60,ncol=90))$Predicted
mine<-as.vector(sapply(1:15,function(m){rep(m,4)}))
sum(mine==outsb)
## Result is 9

news_rda<-data.frame(Reduce("rbind",news2))
news_rda=cbind(news_rda,as.factor(file$V91[as.numeric(rownames(news_rda))]))
names(news_rda)[91]="class"

mins_rda=0
alpha=seq(.001,.999,.001)
lambda=alpha
for(j in 1:length(lambda)){
for(k in length(alpha)){
x <- rda(class ~ ., data = news_rda, gamma = alpha[k], lambda = lambda[j])
mins_rda=max(sum(mine==predict(x, news_test_1)$class),mins_rda)
}
}
mins_rda
## result is 36
