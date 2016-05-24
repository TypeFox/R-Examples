Prgcm_HW1 <-
function(vmc,p){    
                           v1<-ifelse(vmc[,1]==0 & vmc[,2]==0,1/(1+exp(p)),0)+
                           ifelse(vmc[,1]==0 & vmc[,2]==1,exp(p)/(1+exp(p)),0)+
                           ifelse(vmc[,1]==1 & vmc[,2]==0,1/(2*(1+exp(p))),0)+
                           ifelse(vmc[,1]==1 & vmc[,2]==1,1/2,0)+
                           ifelse(vmc[,1]==1 & vmc[,2]==2,exp(p)/(2*(1+exp(p))),0)+
                           ifelse(vmc[,1]==2 & vmc[,2]==1,1/(1+exp(p)),0)+
                           ifelse(vmc[,1]==2 & vmc[,2]==2,exp(p)/(1+exp(p)),0) 
                           return(v1)
                           }
