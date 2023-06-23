library(tidyverse)
# Experiment : Directional watson

methods_array <- c(1,2)
kernel_array <- c(1) #here only use gaussian kernel 
s_array <-  c(20,50,100,200)

results <- data.frame(method=c(),k_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())
for(method_id in methods_array){
  for(kernel in kernel_array){
    
    for(s in s_array){
      
      d1<- readRDS(paste0("Directional_watson_results_s_",s,"_method_",method_id,"_kerneltyp_",kernel,".rds"))  
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
# d1<- readRDS(paste0("Directional_watson_results_s_",20,"_method_",2,"_kerneltyp_",1,".rds"))   #remove this later
# d1<-d1%>%mutate(method=ifelse(method=='spectral_showalter','Showalter','KS'))
# results <-rbind(results,d1)
###########################





#add MMAgg results 
results_mmd <- read.csv("/storage/work/o/oih3/MMDagg/user/raw/results_directional_watson.csv")
d1<-data.frame(method='MMDAgg',k_indx=results_mmd$kindx,d=results_mmd$d,sample_size=results_mmd$n,s_size=0,power=results_mmd$power,kernel=1)
d1 <- d1%>%filter(k_indx>0.5)


results <- rbind(results,d1)

saveRDS(results,'results/combined_results_watson.rds')
# Comparisons versus other tests 
results <-readRDS('results/combined_results_watson.rds') 
df_test <- results %>% filter(s_size==20 | s_size==0)
df_test <- df_test %>% filter(d==10 | d==20)

myplot<-df_test %>% ggplot(aes(x=reorder(k_indx, -power),y=power,group=method,color=method))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,ncol = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Concentration parameter")+ylab("Power")+theme(text = element_text(size = 20))



#Different values of S

df_s <- results%>% filter(method=='Showalter')
df_s <- df_s %>% mutate(method_comb=paste0(method,"_s_",s_size))
df_s <- df_s %>% filter(d==10 | d==20)
myplot<-df_s %>% ggplot(aes(x=reorder(k_indx, -power),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,ncol = 2,labeller = "label_both")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Concentration parameter")+ylab("Power")+theme(text = element_text(size = 18))


###############################################################
#Experiment varying sample size

# increasing sample_size experiments
d=20
k_indx_array = c(6)
method_id=2
kernel=1
results <- data.frame(method=c(),kernel=c(),k_indx=c(),d=c(),sample_size=c(),s_size=c(),power=c())
s_array <-  c(20,50,100,200)
sample_size_array <- c(800,1200,1600,2000)
for (k_indx in k_indx_array){
  for(sample_size in sample_size_array){
    for(s in s_array){
      if(file.exists(paste0("Sample_size_exp/Directional_watson_varying_sampleSize_results_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size",sample_size,"K_indx_",k_indx,".rds"))){
        d1<- readRDS(paste0("Sample_size_exp/Directional_watson_varying_sampleSize_results_s_",s,"_method_",method_id,"_kerneltyp_",kernel,"sample_size",sample_size,"K_indx_",k_indx,".rds"))  
        
        if(s==20){
          d1<-d1%>%mutate(method=ifelse(method=='spectral','Showalter','Energy'))
        }
        else if(s==50){
          d1<-d1%>%mutate(method=ifelse(method=='spectral','Showalter','KS'))
        }
        else{
          d1 <-d1%>%filter(method=='spectral')
          d1<-d1%>%mutate(method='Showalter')
        }
        
        
        results <- rbind(results,d1)
        
      }else{
        print(paste0(" s= ",s," sample_size= ",sample_size," d= ",d))
      }
    }
  }
  
}

#add data for sample_size 500

df2 <- data.frame(method='Energy',kernel=0,k_indx=6,d=20,sample_size=500,s_size=0,power=0.055)
df3 <- data.frame(method='KS',kernel=0,k_indx=6,d=20,sample_size=500,s_size=0,power=0.105)
df4 <- data.frame(method='Showalter',kernel=1,k_indx=6,d=20,sample_size=500,s_size=c(20,50,100,200),power=c(0.530,0.510,0.510,0.440))
results <- rbind(results,df2,df3,df4)


#add MMDAgg results
results_mmd <- read.csv("/storage/work/o/oih3/MMDagg/user/raw/results_directional_watson.csv")
d1<-data.frame(method='MMDAgg',k_indx=results_mmd$kindx,d=results_mmd$d,sample_size=results_mmd$n,s_size=0,power=results_mmd$power,kernel=1)
d1 <- d1%>%filter(k_indx>0.5)
results<-rbind(results,d1)
for(n_size in c(800,1200,1600,2000)){
  results_mmd <- read.csv(paste0("/storage/work/o/oih3/MMDagg/user/raw/results_directional_watson_varying_size",n_size,".csv"))
  d1<-data.frame(method='MMDAgg',k_indx=results_mmd$kindx,d=results_mmd$d,sample_size=results_mmd$n,s_size=0,power=results_mmd$power,kernel=1)
  d1 <- d1%>%filter(k_indx>0.5)
  results<-rbind(results,d1)
}



results <- results%>% filter(d==20)
results <- results%>% filter(k_indx==6)
saveRDS(results,'results/combined_results_watson_increasing_sampleSize.rds')
df_plot <- results %>% mutate(method_comb=paste0(method,"_s_",s_size))

df_plot <-df_plot%>% mutate(method_comb=ifelse(method=='Energy','Energy',method_comb))
df_plot <-df_plot%>% mutate(method_comb=ifelse(method=='KS','KS',method_comb))
df_plot <-df_plot%>% mutate(method_comb=ifelse(method=='MMDAgg','MMDAgg',method_comb))

myplot<-df_plot %>% ggplot(aes(x=as.factor(sample_size),y=power,group=method_comb,color=method_comb))+geom_point(size=2.5)+geom_line(size=1)+guides(color=guide_legend(title=NULL))  + facet_wrap(~d,ncol= 2,labeller = "label_both",scales="free_x")+guides(color=guide_legend(title=NULL))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=1),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="bottom")+xlab("Sample_size")+ylab("Power")+theme(text = element_text(size = 18))












