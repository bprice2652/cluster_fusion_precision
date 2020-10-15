
### This code runs the pulmonary hypertension data provided in the manuscript.
## We first load the files do data cleaning, and then relate it use our AIC metric.
## We note we had tried using cross validaiton so the fold selection is still in the code. 

###############
## Need to source cluster_da.R from file system
source("cluster_da.R")
###################
#below two lines are needed for KU server
library(RevoUtilsMath)
setMKLthreads(1)
#source("~/dropbox/cluster_inv_cov/r_code/cluster_da.R")
#setwd("~/dropbox/cluster_inv_cov/genomics_example")
gds5499 <- readRDS("full_gds5499data.RDS")

genes_sd <- apply(gds5499$expression_data,1,sd)

bad_genes <- which(is.na(genes_sd))

x <- gds5499$expression_data[-bad_genes,]

#simple check if they data is on wildly different scales and that is not the case
# mins <- sort(apply(x,1,min))
# maxs <- sort(apply(x,1,max))
# mins[1]
# mins[length(mins)]
# maxs[1]
# maxs[length(maxs)]

x <- t(x) #puts it on the standard n x p format
scale_x <- scale(x) # if you so desire. 
x <- scale_x

y <- gds5499$y_labels
xtabs(~y)

small_level <- levels(y)[4]
remove_samples <- which(y==small_level)

y <- y[-remove_samples]
x <- x[-remove_samples,]

y <- factor(paste(y))
y <- as.numeric(y)

get_anova_p_value <- function(x,y){
	a <- aov(x ~ y)
	summary(a)[[1]][1,5]
}

get_top_genes <- function(p_values, k){
	ord <- order(p_values)
	good_spots <- ord[1:k]
	good_spots
}

randomly_assign <- function(n,k){
#randomly assign n samples into k groups
   small_set <- floor(n/k)
   group_assign <- NULL
  if(n %% k == 0){
     group_assign <-  rep(seq(1,k),n/k)
   } else{
     remainder <- n %% k
     for(i in 1:remainder){
        group_assign <- c(group_assign, rep(i,small_set+1))
     }
     group_assign <- c(group_assign, rep(seq((i+1),k),small_set))
   }
   sample(group_assign)
}

generate_folds <- function(x_list,fold_num){
	fold=vector("list",fold_num)
	r <- length(x_list)
	
	assignments <- vector("list",r)
	index_adjust <- 0
	
	for(i in 1:r){
		r_n <- dim(x_list[[i]])[1]
		r_index <- ( randomly_assign(r_n,fold_num) + index_adjust) %% fold_num
		r_index[which(r_index==0)] <- fold_num
		group_n <- xtabs(~r_index)
		if(sd(group_n) != 0){
			#this group was not balanced. Adjustment made to balance the overall training/testing split.
			most_vals <- which(group_n==max(group_n))
			index_adjust <- most_vals[length(most_vals)]
		}
		
		assignments[[i]] <- r_index
	}
	
	
	for(i in 1:fold_num){
		fold[[i]] <- vector("list",r)
		for(j in 1:r){
			fold[[i]][[j]] <- which(assignments[[j]]==i)
		}
	}
	fold
}

#run_analysis <- function(seed,x,y,train_p=.8,fold_num=5){
	seed <- 1
	fold_num <- 5
	set.seed(seed)
	n <- length(y)
	r <- length(unique(y))
	
	
	anova_p_values <- apply(x,2,get_anova_p_value,y)
	keeper_genes <- get_top_genes(anova_p_values, n)  # made it 10 to test code
	
	x <- x[,keeper_genes]
	
	x_split <- split(data.frame(x),y)
	cv_folds <- generate_folds(x_split,fold_num)
lambda1=2^seq(1,3,.1)
lambda2=10^seq(1,2,.1)
pcen_aic=clust_inv_cov_aic(X=x_split,lambda1 =lambda1,lambda2=lambda2,Q_set=c(2,3),rho,tol=10^-8,ncores_Q=1,ncores_cv=20,ncores_horse=1,returnAll = TRUE)



#load("bp_pcen_aic_422.RData")
omegs_outs<-pcen_aic$ests[[250]]$omega

