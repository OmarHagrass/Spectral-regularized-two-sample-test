
library(rdist)
library(tidyverse)
library(expm) 
library(parallel)
library(movMF)


generate_samples<- function(n,typ,k,d){    #we can change the null distribution here
  
  mu <- rep(1,d)
  mu <- mu/sqrt(sum(mu^2))
  if(typ==1){
    samples_x <-mclapply(rep(1,n),rmovMF,theta=0*mu,mc.cores = 1) # draw samples
    samples_x <- t(matrix(unlist(samples_x),d,n))
  }else{
    
    samples_x <-mclapply(rep(1,n),rmovMF,theta=k*mu,mc.cores = 1) # draw samples
    samples_x <- t(matrix(unlist(samples_x),d,n))
  }
  
  
}





RNGkind("L'Ecuyer-CMRG")
set.seed(1247)


iter<-  200

n <-  500

dim_arr<-  c(10,20,50,100)

sim_iter <- function(itr,k_indx,d){
  set.seed(23+itr*5+7*k_indx)
  samples_x <- generate_samples(n=n,k=k_indx,typ=1,d=d)
  return(t(samples_x))
}

for(idx in dim_arr){
  
  data <- mclapply(seq(1,iter),sim_iter,k_indx=0,d=idx,mc.cores = 20)
  data_arr <- t(matrix(unlist(data),idx,n*iter))
  write.csv(data_arr,file=paste0("samples_x_dim",idx,".csv"))
}



