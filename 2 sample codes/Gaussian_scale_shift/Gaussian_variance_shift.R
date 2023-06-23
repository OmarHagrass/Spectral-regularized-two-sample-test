library(rdist)
library(tidyverse)
library(expm) 
library(parallel)
library(movMF)
library(energy)
library(fasano.franceschini.test)

source("spectral_test.R")

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


