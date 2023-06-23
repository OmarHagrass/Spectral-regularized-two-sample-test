library(tidyverse)
# Experiment : pert_uni

methods_array <- c(1,2)
kernel_array <- c(1,2) 



results <- data.frame(method=c(),kernel=c(),P_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())

for(d in c(1,2)){
  if(d==1){
    sample_size = 500
    s_array <-  c(20,50,100)
  }
  if(d==2){
    sample_size =2000
    s_array <-  c(20,50,100,200)
  }
  for(method_id in methods_array){
    for(kernel in kernel_array){
      
      for(s in s_array){
        if(file.exists(paste0("Perturbed_uniform_results_d_",d,"_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size_",sample_size,".rds"))){
        d1<- readRDS(paste0("Perturbed_uniform_results_d_",d,"_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size_",sample_size,".rds"))  
        }
        if(method_id==1){
          d1<-d1%>%mutate(method=ifelse(method=='spectral'| method=='spectral_showalter','Tikhonov','Energy'))
          if(s!=20){
            d1 <- d1%>%filter(method=='Tikhonov')  
          }
        }
        if(method_id==2){
          d1<-d1%>%mutate(method=ifelse(method=='spectral'| method=='spectral_showalter','Showalter','KS'))
          if(s!=20){
            d1 <- d1%>%filter(method=='Showalter')  
          }
        }
        
        results <- rbind(results,d1)
      }
    }
  }
} 
results <- results %>% filter(power >=0)
########################## Remove
# d1<- readRDS(paste0("Perturbed_uniform_results_d_",1,"_s_",20,"_method_",2,"_kerneltyp_",1,"sample_size_",500,".rds"))   #remove this later
# d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
# results <-rbind(results,d1)
# d2<- readRDS(paste0("Perturbed_uniform_results_d_",2,"_s_",20,"_method_",2,"_kerneltyp_",1,"sample_size_",2000,".rds"))   #remove this later
# d2<-d2%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
#results <-rbind(results,d2)
###########################





#add MMAgg results 
power_MMD_gaussian <- c(1,0.99,0.422,0.064,0.006,0.004) #d=1
power_MMD_laplace <- c(1,0.998,0.358,0.03,0.006,0) #d=1

d1<-data.frame(method='MMDAgg',kernel=1,P_indx=c(1,2,3,4,5,6),d=1,sample_size=500,s_size=0,power=power_MMD_gaussian)
d2<-data.frame(method='MMDAgg',kernel=2,P_indx=c(1,2,3,4,5,6),d=1,sample_size=500,s_size=0,power=power_MMD_laplace)
results <- rbind(results,d1,d2)
power_MMD_gaussian <- c(1,0.874,0.014) #d=2
power_MMD_laplace <- c(1,0.62,0.004) #d=2

d1<-data.frame(method='MMDAgg',kernel=1,P_indx=c(1,2,3),d=2,sample_size=2000,s_size=0,power=power_MMD_gaussian)
d2<-data.frame(method='MMDAgg',kernel=2,P_indx=c(1,2,3),d=2,sample_size=2000,s_size=0,power=power_MMD_laplace)
results <- rbind(results,d1,d2)




saveRDS(results,'results/combined_results_pert_uni.rds')
# Comparisons versus other tests 
results <-readRDS('results/combined_results_pert_uni.rds') 

df_test <- results %>% filter(s_size==20 | s_size==0)
df_test <- df_test%>%mutate(kernel=ifelse(kernel==0 | kernel==1,'Gaussian','Laplace'))
#df_test <- df_test %>% filter(kernel!=2)

myplot<-df_test %>% ggplot(aes(x=as.factor(P_indx),y=power,group=method,color=method))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(kernel~d,ncol = 2,labeller = labeller(d=label_both,kernel=label_value,.multi_line=FALSE),scales = "free_x")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Perturbations")+ylab("Power")+theme(text = element_text(size = 20))



#Different values of S

df_s <- results%>% filter(method=='Showalter')
df_s <- df_s%>%mutate(kernel=ifelse(kernel==0 | kernel==1,'Gaussian','Laplace'))
df_s <- df_s %>% mutate(method_comb=paste0(method,"_s_",s_size))

myplot<-df_s %>% ggplot(aes(x=as.factor(P_indx),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(kernel~d,ncol= 2,labeller = labeller(d=label_both,kernel=label_value,.multi_line=FALSE),scales="free_x")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Perturbations")+ylab("Power")+theme(text = element_text(size = 18))



#####################################################

# increasing sample_size experiments
d=2 
P_indx_array = c(2,3)
method_id=2
kernel=1
results <- data.frame(method=c(),kernel=c(),P_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())
s_array <-  c(20,50,100,200)

for (P_indx in P_indx_array){
  if(P_indx==2){
    sample_size_array = c(500,800,1100,1400,1700)
  }
  if(P_indx==3){
    sample_size_array = c(2500,3000,3500,4000)
  }
  for(sample_size in sample_size_array){
    for(s in s_array){
      if(file.exists(paste0("Sample_size_exp/Perturbed_uniform_varying_sampleSize_results_d_",d,"_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size_",sample_size,"Pindx_",P_indx,".rds"))){
        d1<- readRDS(paste0("Sample_size_exp/Perturbed_uniform_varying_sampleSize_results_d_",d,"_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size_",sample_size,"Pindx_",P_indx,".rds"))  
        d1<-d1%>%mutate(method='Showalter')
        results <- rbind(results,d1)
      
      }else{
        print(paste0(" s= ",s," sample_size= ",sample_size," d= ",d))
      }
  }
  }
  
}
saveRDS(results,'results/combined_results_pert_uni_increasing_sampleSize.rds')
df_plot <- results %>% mutate(method_comb=paste0(method,"_s_",s_size))

myplot<-df_plot %>% ggplot(aes(x=as.factor(sample_size),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~P_indx,ncol= 2,labeller = "label_both",scales="free_x")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Sample_size")+ylab("Power")+theme(text = element_text(size = 18))









