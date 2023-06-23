library(watson)

set.seed(1247)
d_array=c(10,20,50,100)
k_indx_array= c(8,6,4,2,1)
for(d in d_array){
iter=200
n=500
samples_x<-rmwat(n = n*iter, weights = c(0.5,0.5), kappa = c(0,0),
                mu =matrix(c(rep(1,d),-1,rep(1,d-1)),nrow=d))

write.csv(samples_x,file=paste0("samples_x_watson_dim",d,".csv"))


for(k_indx in k_indx_array){
  samples_y<-rmwat(n = n*iter, weights = c(0.5,0.5), kappa = c(k_indx,k_indx),
                   mu =matrix(c(rep(1,d),-1,rep(1,d-1)),nrow=d))
  
  write.csv(samples_y,file=paste0("samples_y_watson_dim",d,"_k_indx_",k_indx,".csv"))
} 

}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

dat <- circleFun(c(0,0),2,npoints = 500)
d=2 

samples<-rmwat(n = 500, weights = c(0.5,0.5), kappa = c(10,10),
                 mu =matrix(c(rep(1,d),-1,rep(1,d-1)),nrow=d))
df <- data_frame(samples_x=samples[,1],samples_y=samples[,2])

df%>%ggplot()+geom_point(aes(x=samples_x,y=samples_y,colour="black"))+xlim(-1,1)+geom_path(aes(x=dat$x,y=dat$y,colour = "blue"))+theme(legend.position = "none")




