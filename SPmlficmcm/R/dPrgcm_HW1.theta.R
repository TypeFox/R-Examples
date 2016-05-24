dPrgcm_HW1.theta <-
function(vmc,p){    
                           v1<-ifelse(vmc[,1]==0 & vmc[,2]==0,-S1(p),0)+
                           ifelse(vmc[,1]==0 & vmc[,2]==1,S1(p),0)+
                           ifelse(vmc[,1]==1 & vmc[,2]==0,-S1(p)/2,0)+
                           ifelse(vmc[,1]==1 & vmc[,2]==2,S1(p)/2,0)+
                           ifelse(vmc[,1]==2 & vmc[,2]==1,-S1(p),0)+
                           ifelse(vmc[,1]==2 & vmc[,2]==2,S1(p),0) 
                           return(v1)
                           }
