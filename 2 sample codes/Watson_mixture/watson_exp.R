
library(rdist)
library(tidyverse)
library(expm) 
library(parallel)
library(movMF)
library(energy)
library(fasano.franceschini.test)
source("spectral_test.R")

sim_iter <- function(itr,indx){ #indx is the concentration parameter in von-mises distributon  
  
  set.seed(23+itr*5+7*indx)
  fr <- 1+n*(itr-1)
  to <- n*itr
  x <- matrix(unlist(data_x[fr:to,]),n,d)
  y <- matrix(unlist(data_y[fr:to,]),n,d) 
  
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

#args<- c(50,10,500,1)
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
results <- data.frame(method=c(),kernel=c(),k_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())

dim_array <-  c(10,20)

kindx_array <- c(8,6,4,2,1)

for(d in dim_array){
  data_x <- read_csv(paste0("/storage/work/oih3/2-samples_codes/directional_watson_mixture/data/samples_x_watson_dim",d,".csv"))[-1]
  for(k_indx in kindx_array){
    
    data_y <- read_csv(paste0("/storage/work/oih3/2-samples_codes/directional_watson_mixture/data/samples_y_watson_dim",d,"_k_indx_",k_indx,".csv"))[-1]
    output <- mclapply(1:iter,sim_iter,indx=k_indx,mc.cores =20 )
    output<- matrix(unlist(output),2,iter)
    
    error_KS <- sum(output[1,])/iter  
    
    error_spectral <- sum(output[2,])/iter
    
    df1 <- data.frame(method="KS_test",kernel=0,k_indx=k_indx,d=d,sample_size=n,s_size=0,power=1-error_KS)
    df2 <- data.frame(method="spectral_showalter",kernel=kernel_type,k_indx=k_indx,d=d,sample_size=n,s_size=s,power=1-error_spectral)
    
    print(paste0("watson_s_",s,"_d_",d,"_n_",n,"k_indx",k_indx))
    
    print(paste0("error_KS=",error_KS,"error_spectral=",error_spectral))
    
    
    results<-rbind(results,df1,df2)
  }
  end_time = Sys.time()
  
  print(end_time-start_time)  
}
saveRDS(results,paste0("Directional_watson_results_s_",s,"_method_",method,"_kerneltyp_",kernel_type,".rds"))
