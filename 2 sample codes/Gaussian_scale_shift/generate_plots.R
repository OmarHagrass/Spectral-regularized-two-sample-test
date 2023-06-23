
library(tidyverse)
# Experiment : Gaussian shifts

methods_array <- c(1,2)
kernel_array <- c(1) #here only use gaussian kernel 
s_array <-  c(20,50,100)

results <- data.frame(method=c(),shift=c(),d=c(),sample_size=c(),s_size=c(),power=c())
for(method_id in methods_array){
    for(kernel in kernel_array){
      
      for(s in s_array){
      
        d1<- readRDS(paste0("Gaussian_shifts_results_s_",s,"_method_",method_id,"_kerneltyp_",kernel,".rds"))  
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
########################## Remove
# d1<- readRDS(paste0("Gaussian_shifts_results_s_",20,"_method_",2,"_kerneltyp_",1,".rds"))   #remove this later
# d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
# results <-rbind(results,d1)
###########################



results <- results %>% filter(shift!=0)

#add MMAgg results 
results_mmd <- read.csv("/storage/work/o/oih3/MMDagg/user/raw/results_gaussian_shifts.csv")
d1<-data.frame(method='MMDAgg',shift=results_mmd$shift,d=results_mmd$d,sample_size=results_mmd$n,s_size=0,power=results_mmd$power)

results <- rbind(results,d1)

saveRDS(results,'results/combined_results_gaussian_shifts.rds')
# Comparisons versus other tests 
results <-readRDS('results/combined_results_gaussian_shifts.rds') 
df_test <- results %>% filter(s_size==20 | s_size==0)


myplot<-df_test %>% ggplot(aes(x=reorder(shift, -power),y=power,group=method,color=method))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,nrow = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Mean shift")+ylab("Power")+theme(text = element_text(size = 20))



#Different values of S

df_s <- results%>% filter(method=='Showalter')
df_s <- df_s %>% mutate(method_comb=paste0(method,"_s_",s_size))
myplot<-df_s %>% ggplot(aes(x=reorder(shift, -power),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,nrow = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Mean shift")+ylab("Power")+theme(text = element_text(size = 18))



#################################################
#Experiment: Gaussian_scale_variance

methods_array <- c(1,2)
kernel_array <- c(1) #here only use gaussian kernel 
s_array <-  c(20,50,100)

results <- data.frame(method=c(),sd_scale=c(),d=c(),sample_size=c(),s_size=c(),power=c())
for(method_id in methods_array){
  for(kernel in kernel_array){
    
    for(s in s_array){
      
      d1<- readRDS(paste0("Gaussian_sd_scale_results_s_",s,"_method_",method_id,"_kerneltyp_",kernel,".rds"))  
      if(method_id==1){
        d1<-d1%>%mutate(method=ifelse(method=='spectral','Tikhonov','Energy'))
        if(s!=20){
          d1 <- d1%>%filter(method=='Tikhonov')  
        }
      }
      if(method_id==2){
        d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
        if(s!=20){
          d1 <- d1%>%filter(method=='Showalter')  
        }
      }
      
      results <- rbind(results,d1)
    }
  }
}
########################## Remove
# d1<- readRDS(paste0("Gaussian_sd_scale_results_s_",20,"_method_",2,"_kerneltyp_",1,".rds"))   #remove this later
# d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
# results <-rbind(results,d1)
###########################



results <- results %>% filter(sd_scale!=1)


results <- results%>%mutate(sd_scale=log10(sd_scale))
#add MMAgg results 
results_mmd <- read.csv("/storage/work/o/oih3/MMDagg/user/raw/results_gaussian_scales.csv")
d1<-data.frame(method='MMDAgg',sd_scale=results_mmd$sd_scale,d=results_mmd$d,sample_size=results_mmd$n,s_size=0,power=results_mmd$power)

results <- rbind(results,d1)

saveRDS(results,'results/combined_results_gaussian_var_scale.rds')
# Comparisons versus other tests 
results <-readRDS('results/combined_results_gaussian_var_scale.rds') 
df_test <- results %>% filter(s_size==20 | s_size==0)

df_test <- df_test %>% filter(sd_scale>=0.1)
df_test <- df_test %>% mutate(sd_scale=sd_scale-0.005)


myplot<-df_test %>% ggplot(aes(x=reorder(sd_scale, -power),y=power,group=method,color=method))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,nrow = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Standard deviation (log scale)")+ylab("Power")+theme(text = element_text(size = 20))



#Different values of S

df_s <- results%>% filter(method=='Showalter')
df_s <- df_s %>% filter(sd_scale>=0.1)
df_s <- df_s %>% mutate(sd_scale=sd_scale-0.005)
df_s <- df_s %>% mutate(method_comb=paste0(method,"_s_",s_size))
myplot<-df_s %>% ggplot(aes(x=reorder(sd_scale, -power),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,nrow = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Standard deviation (log scale)")+ylab("Power")+theme(text = element_text(size = 18))








