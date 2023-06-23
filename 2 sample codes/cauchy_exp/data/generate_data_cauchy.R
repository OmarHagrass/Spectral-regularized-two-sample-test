library(watson)

set.seed(1247)
#d_array=c(1,10)
d_array=c(20,50)
#shift_array= c(1,0.7,0.5,0.3,0.1,0.05,0)

shift_array= c(1,0.7,0.5,0.3,0.1,0.05)
for(d in d_array){
  iter=200
  n=500
  samples_x<-matrix(rcauchy(n=n*iter*d,location=0,scale=1),nrow=n*iter,ncol = d)
  
  write.csv(samples_x,file=paste0("samples_x_cauchy_dim",d,".csv"))
  
  
  for(shift in shift_array){
    samples_y<-matrix(rcauchy(n=n*iter*d,location=shift,scale=1),nrow=n*iter,ncol = d)
    
    write.csv(samples_y,file=paste0("samples_y_cauchy_dim",d,"_shift_",shift,".csv"))
  } 
  
}



