library(rdist)
library(tidyverse)
library(expm) 
library(parallel)

fun_G <- function(x){
  
  if ((-1 < x) & (x < -0.5)){
    return(exp(-1 / (1 - (4 * x + 3)^2)))
  }
  if ((-0.5 < x) & (x < 0)){
    return (-exp(-1 / ( 1 - (4 * x + 1)^2)))   
  }
  return(0)  
  
}
pert_uni <- function(x,p,ss,pert_mult,theta,V){
  output <-0
  for(i in 1:(p^d)){
    v <- V[i,]
    prod <- 1
    for(j in 1:d){
      prod <- prod * theta[i]*fun_G(p*x[j]-v[j])
    }
    output <- output + prod
  }
  output <- output * pert_mult * p^(-ss)
  if( (min(x)>=0) & (max(x)<=1)){
    output=output+1
  }
  return(output)
}

rej_samp <- function(n,p,ss,pert_mult,theta,V){
  samp <- c()
  cnt <-0
  density_max = 1 + pert_mult * (p^(-ss)) *exp(-d)
  while(cnt < n){
    unif <- runif(1,0,density_max)
    Y <- runif(d,0,1)
    
    if(unif <= pert_uni(Y,p,ss,pert_mult,theta,V)){
      cnt<-cnt+1
      #print(cnt)
      samp <- cbind(samp,Y)
      
    }
  }
  return(unname(samp))
}
generate_samples<- function(n,p,typ,theta,V){    #we can change the null distribution here
  a=0
  b=1
  if(d==1){
    pert_mult <- 2.7
  }else{
    pert_mult <- 7.3
  }
  if(typ==1){
    samples_x <- mclapply(rep(d,n), runif,min=a,max=b,mc.cores=1) # draw samples
    samples_x <- matrix(unlist(samples_x),d,n)
  }else{
    
    #theta <- rep(1,p^d)
    #theta <- c(1,-1)
    samples_x<-rej_samp(n,p,1,pert_mult,theta,V)
    
  }
  
  return(samples_x)
  
}




d=2   #d=1
n=2000  #n=500
m=n
iter=200
flg=0
for(indx in 1:6){
  
  p<- indx
  V <- expand.grid(rep(list(1:p),d))
  V <- array(unlist(V),dim=c(p^d,d))
  samples_x <- c()
  samples_y <- c()
  for (itr in 1:iter){
    set.seed(23+itr*5+7*indx)
    p <-indx
    theta <- sample(c(1,-1),p^d,replace=TRUE)
    if(flg==0){
    x <- t(generate_samples(n=n,p=indx,typ=1,theta=theta,V=V))
    samples_x  <- rbind(samples_x,x)
    }
    y <- t(generate_samples(n=m,p=indx,typ=2,theta=theta,V=V))
    samples_y <- rbind(samples_y,y)
  }
  if(flg==0){
  write.csv(samples_x,file=paste0("samples_x_pert_uni_dim",d,".csv"))
  flg=1
  }
  write.csv(samples_y,file=paste0("samples_y_pert_uni_dim",d,"_indx_",indx,".csv"))
}

