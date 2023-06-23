library(rdist)
library(tidyverse)
library(expm) 
library(parallel)
library(movMF)
library(energy)
library(fasano.franceschini.test)
######################################## functions for spectral regularized test

kernel_fun <- function(dist,h){
  if (kernel_type==1){#gaussian
    ret <- exp((-1/2)*(1/h)*dist^2)
  }
  if (kernel_type==2){#laplace
    ret <- exp((-1/h)*dist)  
  }
  return(ret) 
}

get_metric <- function(idx){
  if (idx==1){ #gaussian
    ret <- "euclidean"
  }
  if (idx==2){ #laplace
    ret <- "manhattan"
  }
  return(ret) 
}
get_h_Ks <- function(samples_x,samples_y, samples_z){
  ZdistZ <- pdist(t(samples_z),metric = get_metric(kernel_type))
  agg_samples <- cbind(samples_x,samples_y)
  agg_dist <- pdist(t(agg_samples),metric = get_metric(kernel_type))
  
  
  #band width
  if(kernel_type==1){
    h<- median(agg_dist)^2
  }
  if(kernel_type==2){
    h <-  median(agg_dist)
  }
  K_s <- kernel_fun(ZdistZ,h) 
  ret <- list(h=h,K_s=K_s)
}

get_cov_eigs <- function(Ks,h_mult){
  
  K_s<- Ks^(1/h_mult)
  # other matrices
  Hs <- diag(s) - (1/s) * rep(1,s)%*%t(rep(1,s))
  Hsd <- (s/(s-1)) * Hs
  
  
  #computing eigen vectors and values
  H1 <- Re(sqrtm(Hsd))
  rank <- rankMatrix(H1%*%K_s%*%H1)[1]
  eig<-eigen(H1%*%K_s%*%H1)
  eig_val <- eig$values[1:rank]
  eig_val<-Re(eig_val)/s
  
  eig_vec <- Re(eig$vectors[,1:rank])
  ret <- list(H1=H1,eig_val=eig_val,eig_vec=eig_vec)
  return(ret)
  
  
}
evaluate_matrices<-function(samples_y,samples_x,samples_z,h){
  
  
  XdistX <- pdist(t(samples_x),metric = get_metric(kernel_type))
  
  
  
  XdistZ <- cdist(t(samples_x),t(samples_z),metric = get_metric(kernel_type))
  YdistZ <- cdist(t(samples_y),t(samples_z),metric = get_metric(kernel_type))
  YdistX <- cdist(t(samples_y),t(samples_x),metric = get_metric(kernel_type))
  YdistY <- pdist(t(samples_y),metric = get_metric(kernel_type))
  
  
  #Kernel choice
  
  K_n <- kernel_fun(XdistX,h) 
  K_mn <- kernel_fun(YdistX,h) 
  K_ns <- kernel_fun(XdistZ,h) 
  K_ms <- kernel_fun(YdistZ,h) 
  K_m <- kernel_fun(YdistY,h) 
  
  
  ret <- list(K_n=K_n,K_ms=K_ms,K_mn=K_mn,K_ns=K_ns,K_m=K_m)
  return(ret)
}

G_ld<-function(x,ld){
  
  ret1 <- ((x+ld)^(-1))
  ret2 <-  (if_else(x!=0,((1-exp(-x/ld))/x),1/ld))
  return(switch(method,ret1,ret2))
}
evaluate_stat_ld<- function(ld,H1,K_n,K_m,K_ms,K_mn,K_ns,eig_val,eig_vec,n,m){
  
  Gld_sig = eig_vec %*% diag( (1/eig_val)*(G_ld(eig_val,ld)-G_ld(0,ld) ) ) %*% t(eig_vec)
  
  mat1 <- (G_ld(0,ld)*K_n)+(1/s)*K_ns %*% H1 %*% Gld_sig %*% H1 %*%  t(K_ns)
  mat2 <- (G_ld(0,ld)*K_m)+(1/s)*K_ms %*% H1 %*% Gld_sig %*% H1 %*%  t(K_ms)
  mat3 <- t(rep(1,m)) %*% (G_ld(0,ld)*K_mn+(1/s)*K_ms %*% H1 %*% Gld_sig %*% H1 %*%  t(K_ns)) %*% rep(1,n)
  
  
  eta <- (n*(n-1))^(-1) * (sum(mat1) - sum(diag(mat1)) ) + (m*(m-1))^(-1) * (sum(mat2)-sum(diag(mat2))) - 2*(n*m)^(-1)*mat3
  
  
  
  ret <- eta
  
  return(ret)
}

get_Lambda <- function(lower,upper){
  
  
  #ld_low <- (n/log2(log2(n)))^-1
  ld_low <- 10^(-lower)
  ld_up <- upper
  
  card <- ceiling(log2(ld_up/ld_low))
  Lambda <- 2^seq(0,card)
  Lambda <- ld_low*Lambda
  return(Lambda)
}

eval_stat<- function(n, m, cov_eigs, samples_x, samples_y, samples_z,Lambda, perm_id,h){
  
  
  
  
  if(d==1){
    samples_x<-t(samples_x)
    samples_y<-t(samples_y)
    
  }
  perm_samples_x <- samples_x
  perm_samples_y <- samples_y
  
  if(perm_id>0){
    #set.seed(perm_id*7+23)
    agg_samples <- cbind(samples_x,samples_y)
    len_n <- ncol(samples_x)
    len_m <- ncol(samples_y)
    ind_perm <- seq(1,len_n+len_m)
    ind_perm <-sample(ind_perm)
    agg_samples <- agg_samples[,ind_perm]
    if(d==1){
      agg_samples<-t(agg_samples)
    }
    perm_samples_x <- agg_samples[,1:len_n]
    perm_samples_y <- agg_samples[,(len_n+1):(len_m+len_n)]
    if(d==1){
      perm_samples_x<-t(perm_samples_x)
      perm_samples_y<-t(perm_samples_y)
    }
  }
  
  mat <- evaluate_matrices(samples_z=samples_z,samples_y=perm_samples_y,samples_x=perm_samples_x,h=h)
  
  
  
  stat_results <- c()
  for(i in 1:length(h_mult_arr)){
    bb<- 1/(h_mult_arr[i])
    stat_results <- c(stat_results,unlist(lapply(Lambda,evaluate_stat_ld,H1=cov_eigs[[i]]$H1,K_n=(mat$K_n)^bb,K_m=(mat$K_m)^bb,K_ms=(mat$K_ms)^bb,K_mn=(mat$K_mn)^bb,K_ns=(mat$K_ns)^bb,eig_val=cov_eigs[[i]]$eig_val,eig_vec=cov_eigs[[i]]$eig_vec,n=n,m=m)))
  }
  
  #print(paste0("adp_stat=",adp_stat," perm_id=",perm_id))
  return(stat_results)
}

compute_test <- function(samples_x,samples_y,n,m,s,num_perm, h_arr, Lambda_arr,alpha){
  
  bern <- sample(c(0,1),s,replace=TRUE)
  samples_z <- t(bern*(t(samples_x[,1:s])) + (1-bern)*(t(samples_y[,1:s])))
  if(d==1){
    samples_z <- t(samples_z)
  }
  arr_h_ks <- get_h_Ks(samples_x = samples_x, samples_y= samples_y, samples_z=samples_z)
  cov_eigs <- lapply(h_mult_arr,get_cov_eigs,Ks=arr_h_ks$K_s)
  
  perm_stats <- mclapply(seq(1,num_perm),eval_stat,n=n-s,m=m-s, h=arr_h_ks$h,cov_eigs=cov_eigs, samples_x=samples_x[,(s+1):n], samples_y=samples_y[,(s+1):m], samples_z=samples_z,Lambda=Lambda,mc.cores = 1)
  
  perm_stats <- t(matrix(unlist(perm_stats),length(Lambda)*length(h_mult_arr),num_perm))
  
  stat <- eval_stat(n=n-s,m=m-s, h=arr_h_ks$h ,cov_eigs=cov_eigs, samples_x=samples_x[,(s+1):n], samples_y=samples_y[,(s+1):m], samples_z=samples_z,Lambda=Lambda,perm_id = 0)
  
  thres <- apply(perm_stats,2,sort)[ceiling((1-alpha/(length(Lambda)*length(h_mult_arr)))*num_perm),]
  
  
  err_flg <- min(as.numeric(stat<thres))
  return (err_flg)
}
###################################################



sim_iter <- function(itr){   #s : number of samples to estimate the covariance
  
  set.seed(23+itr*5)
  
  x <- matrix(rnorm(n*d),nrow=n,ncol=d)
  y <- matrix(rnorm(m*d,mean = 0,sd = 1*shift),nrow=m,ncol=d)
  
  booled<- rbind(x,y)
  
  #energy_error <- as.numeric(eqdist.etest(booled,c(n,m),R=199)$p.value > alpha)
  
  if(s==20){
    KS_error<-as.numeric(fasano.franceschini.test(x,y,verbose=FALSE)$p.value > 0.05)
  }
  if(s>20){
    KS_error <- 5  
  }
  spectral_error <- compute_test(samples_x=t(x),samples_y = t(y), n = n,m=m,s=s,num_perm=num_perm,h_arr = h_mult_arr,Lambda_arr = Lambda,alpha = alpha)
  #print(paste0("energy_error=",energy_error," spectral_error=", spectral_error," itr=",itr))
  
  #return(list(energy_error=energy_error,spectral_error=spectral_error))
  return(list(KS_error=KS_error,spectral_error=spectral_error))
}


RNGkind("L'Ecuyer-CMRG")
set.seed(1247)
args <- commandArgs(trailingOnly=T)

#args<- c(3,5,20)


iter<- as.numeric(args[1])
s <- as.numeric(args[2])
n <- as.numeric(args[3])
method <- 2 # 1 is Tikhonov , 2 is shwalter
kernel_type <- as.numeric(args[4]) #1 is gaussian , 2 is laplace
m <- n
num_perm <- 60


Lambda <- 10^seq(-6,1,0.75)

h_low <- 2
h_up <- 2
h_mult_arr <- 10^seq(-h_low,h_up,0.5)            #band_width array
alpha <- 0.05

start_time = Sys.time()
results <- data.frame(method=c(),sd_scale=c(),d=c(),sample_size=c(),s_size=c(),power=c())

dim_array <-  c(1,10,20,50,100)

scales_array <- c(1,10^seq(0.005,0.6,0.1))



#dim_array<- c(1,2)
#shifts_array<-c(1,0.7)
for(d in dim_array){
  for(shift in scales_array){
    
    output <- mclapply(1:iter,sim_iter,mc.cores =20 )
    output<- matrix(unlist(output),2,iter)
    
    error_KS <- sum(output[1,])/iter  
    
    error_spectral <- sum(output[2,])/iter
    
    df1 <- data.frame(method="KS_test",sd_scale=shift,d=d,sample_size=n,s_size=0,power=1-error_KS)
    df2 <- data.frame(method="spectral_showalter",sd_scale=shift,d=d,sample_size=n,s_size=s,power=1-error_spectral)
    
    print(paste0("GaussianSdscales_s_",s,"_d_",d,"_n_",n,"sd_scale",shift))
    
    print(paste0("error_KS=",error_KS,"error_spectral_showalter=",error_spectral))
    
    
    results<-rbind(results,df1,df2)
  }
  end_time = Sys.time()
  
  print(end_time-start_time)  
}
saveRDS(results,paste0("Gaussian_sd_scale_results_s_",s,"_method_",method,"_kerneltyp_",kernel_type,".rds"))


