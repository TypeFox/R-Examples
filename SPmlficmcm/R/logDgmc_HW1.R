logDgmc_HW1 <-
function(vmc,p){
                                  v<-(ifelse(vmc[,1]==0 & vmc[,2]==0,((-3)*exp(p))/(exp(p)+1),0)+ifelse(vmc[,1]==0 & vmc[,2]==1,(1-2*exp(p))/((exp(p)+1)),0)+
                                  ifelse(vmc[,1]==1 & vmc[,2]==0,(1-2*exp(p))/((exp(p)+1)),0)+ifelse(vmc[,1]==1 & vmc[,2]==1,(1-exp(p))/((exp(p)+1)),0)+
                                  ifelse(vmc[,1]==1 & vmc[,2]==2,(2-exp(p))/(exp(p)+1),0)+ifelse(vmc[,1]==2 & vmc[,2]==1,(2-exp(p))/(exp(p)+1),0)+
                                  ifelse(vmc[,1]==2 & vmc[,2]==2,3/(1+exp(p)),0)) 
                                 return(v) 
                                 }
