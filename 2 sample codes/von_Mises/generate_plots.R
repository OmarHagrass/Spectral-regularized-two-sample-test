library(tidyverse)
# Experiment : Directional vonmises

methods_array <- c(1,2)
kernel_array <- c(1) #here only use gaussian kernel 
s_array <-  c(20,50,100,200)

results <- data.frame(method=c(),k_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())
for(method_id in methods_array){
  for(kernel in kernel_array){
    
    for(s in s_array){
      
      d1<- readRDS(paste0("Directional_von_mises_results_s_",s,"_method_",method_id,"_kerneltyp_",kernel,".rds"))  
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
# d1<- readRDS(paste0("Directional_von_mises_results_s_",20,"_method_",2,"_kerneltyp_",1,".rds"))   #remove this later
# d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
# results <-rbind(results,d1)
###########################





#add MMAgg results 
results_mmd <- read.csv("/storage/work/o/oih3/MMDagg/user/raw/results_directional.csv")
d1<-data.frame(method='MMDAgg',k_indx=results_mmd$kindx,d=results_mmd$d,sample_size=results_mmd$n,s_size=0,power=results_mmd$power,kernel=1)
d1 <- d1%>%filter(k_indx>0.5)


results <- rbind(results,d1)

saveRDS(results,'results/combined_results_von_mises.rds')
# Comparisons versus other tests 
results <-readRDS('results/combined_results_von_mises.rds') 
df_test <- results %>% filter(s_size==20 | s_size==0)

myplot<-df_test %>% ggplot(aes(x=reorder(k_indx, -power),y=power,group=method,color=method))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,nrow = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Concentration parameter")+ylab("Power")+theme(text = element_text(size = 20))



#Different values of S

df_s <- results%>% filter(method=='Showalter')
df_s <- df_s %>% mutate(method_comb=paste0(method,"_s_",s_size))
myplot<-df_s %>% ggplot(aes(x=reorder(k_indx, -power),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,nrow = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Concentration parameter")+ylab("Power")+theme(text = element_text(size = 18))
