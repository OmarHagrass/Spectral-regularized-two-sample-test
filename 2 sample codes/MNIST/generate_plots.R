library(tidyverse)
# Experiment : pert_uni

methods_array <- c(1,2)
kernel_array <- c(1,2) 

d=49
s_array <-  c(20,50,100,200)
sample_size <- 500
results <- data.frame(method=c(),kernel=c(),Q_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())



    
  
  for(method_id in methods_array){
    for(kernel in kernel_array){
      
      for(s in s_array){
        if(file.exists(paste0("Mnist_data_results_d_",d,"_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size_",sample_size,".rds"))){
        d1<- readRDS(paste0("Mnist_data_results_d_",d,"_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size_",sample_size,".rds"))
        }
        if(method_id==1){
          d1<-d1%>%mutate(method=ifelse(method=='spectral' | method=='spectral_showalter','Tikhonov','Energy'))
          if(s!=20){
            d1 <- d1%>%filter(method=='Tikhonov')  
          }
        }
        if(method_id==2){
          d1<-d1%>%mutate(method=ifelse(method=='spectral' | method=='spectral_showalter','Showalter','KS'))
          if(s!=20){
            d1 <- d1%>%filter(method=='Showalter')  
          }
        }
        
        results <- rbind(results,d1)
      }
    }
  }
results <- results %>% filter(power >=0)
########################## Remove
# d1<- readRDS(paste0("Mnist_data_results_d_",d,"_s_",20,"_method_",2,"_kerneltyp_",1,"sample_size_",500,".rds"))   #remove this later
# d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
# results <-rbind(results,d1)

###########################





#add MMAgg results 
power_MMD_gaussian <- c(1,0.992,	0.454,0.464,	0.056)
power_MMD_laplace <- c(1,1,0.594,	0.58,	0.094) 

d1<-data.frame(method='MMDAgg',kernel=1,Q_indx=c(1,2,3,4,5),d=d,sample_size=500,s_size=0,power=power_MMD_gaussian)
d2<-data.frame(method='MMDAgg',kernel=2,Q_indx=c(1,2,3,4,5),d=d,sample_size=500,s_size=0,power=power_MMD_laplace)
results <- rbind(results,d1,d2)
##########




saveRDS(results,'results/combined_results_Mnist.rds')
# Comparisons versus other tests 
results <-readRDS('results/combined_results_Mnist.rds') 

df_test <- results %>% filter(s_size==20 | s_size==0)
#df_test <- df_test %>% filter(kernel!=2)
df_test <- df_test%>%mutate(kernel=ifelse(kernel==0 | kernel==1,'Gaussian','Laplace'))
myplot<-df_test %>% ggplot(aes(x=as.factor(Q_indx),y=power,group=method,color=method))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~kernel,ncol = 2,labeller = "label_both",scales = "free_x")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Alternative set")+ylab("Power")+theme(text = element_text(size = 20))



#Different values of S

df_s <- results%>% filter(method=='Showalter')
#df_s <- df_s %>% filter(kernel!=2)
df_s <- df_s %>% mutate(method_comb=paste0(method,"_s_",s_size))
df_s <- df_s%>%mutate(kernel=ifelse(kernel==0 | kernel==1,'Gaussian','Laplace'))
myplot<-df_s %>% ggplot(aes(x=as.factor(Q_indx),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~kernel,ncol= 2,labeller = "label_both",scales="free_x")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Alternative set")+ylab("Power")+theme(text = element_text(size = 18))+ylim(0,1)


