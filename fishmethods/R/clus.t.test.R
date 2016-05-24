clus.t.test<-function(popchar=NULL, cluster=NULL, group=NULL,rho=NULL,alpha=0.05,alternative = c("two.sided")){
  if(nlevels(as.factor(group))!=2) stop("Only two groups allowed in data")
  if(is.null(rho)) stop("rho is required.")
  if(is.null(group)) stop("group vector not specified.")
  if(is.null(cluster)) stop("cluster vector not specified.")
  if(is.null(popchar)) stop("popchar vector not specified.")                                
  if(any(c(length(group),length(popchar),length(cluster)) %in% length(group)=="FALSE"))
    stop("vectors lengths are different.")                                  
  x<-data.frame(cbind(popchar,cluster,group))                                        
  x$popchar<-as.numeric(as.character(x$popchar))
  x<-x[!is.na(x$popchar) & !is.na(x$cluster) & !is.na(x$group),]
  x<-x[x$popchar>0,]
  x$level<-as.integer(as.factor(x$group))-1 
  
  # First group
  x1<-x[x$level==0,] 
  summM1<-aggregate(x1$popchar,list(x1$cluster),length)
  names(summM1)<-c("cluster","M")                                        
  sM<-aggregate(x1$popchar,list(x1$cluster),mean)
  names(sM)<-c("cluster","mean")  
  stats1<-merge(sM,summM1,by.x="cluster",by.y="cluster")
  names(stats1)<-c("cluster","mean","M")
  R1<-mean(x1$popchar)
  SD1<-sqrt(var(x1$popchar)/length(x1$popchar))
  n1<-length(x1$popchar)                                         
  # Second group
   x2<-x[x$level==1,] 
   summM2<-aggregate(x2$popchar,list(x2$cluster),length)
   names(summM2)<-c("cluster","M")    
   sM2<-aggregate(x2$popchar,list(x2$cluster),mean)
   names(sM2)<-c("cluster","mean")                                                
   stats2<-merge(sM2,summM2,by.x="cluster",by.y="cluster")  
   names(stats2)<-c("cluster","mean","M")
   R2<-mean(x2$popchar)
   SD2<-sqrt(var(x2$popchar)/length(x2$popchar))
   n2<-length(x2$popchar)   
  # Calculate t stat adjustment and degrees of freedom
    M1<-n1
    M2<-n2
    A1<-((M1^2)*sum(stats1$M^2)+sum(stats1$M^2)^2-2*M1*sum(stats1$M^3))/(M1^2)
    A2<-((M2^2)*sum(stats2$M^2)+sum(stats2$M^2)^2-2*M2*sum(stats2$M^3))/(M2^2)
    M<-M1+M2
    d1<-exp(log(M1)+log(M))
    d2<-exp(log(M2)+log(M))
    nu<-(sum(stats1$M^2)/(2*M1))+(sum(stats2$M^2)/(2*M2))
    nb<-((M2*sum(stats1$M^2))/d1)+((M1*sum(stats2$M^2))/d2)
    hun<-((M-2)-2*(nu-1)*rho)^2
    hud<-(M-2)*(1-rho)^2+(A1+A2)*rho^2+2*(M-2*nu)*rho*(1-rho)
    hu<-hun/hud                                    
    cu<-sqrt(((M-2)-2*(nu-1)*rho)/((M-2)*(1+(nb-1)*rho))) 
  # T.test and correction 
   hyp<-NULL
  if(alternative=="two.sided"){
    tstat<-abs((R1-R2)/sqrt(SD1^2+SD2^2)) 
    tcrit2<-(qt(1-alpha/2,hu)*SD1^2+qt(1-alpha/2,hu)*SD2^2)/(SD1^2+SD2^2)
    hyp<-"mean 1 = mean 2"                                                             
    p<-(1-pt(tstat*cu,hu))*2                                   
    }  
  if(alternative=="less"){
    tstat<-(R1-R2)/sqrt(SD1^2+SD2^2) 
    tcrit2<-(qt(alpha,hu)*SD1^2+qt(alpha,hu)*SD2^2)/(SD1^2+SD2^2)
     p<-pt(tstat*cu,hu) 
    hyp<-"mean 1 <= mean 2"
  }
  if(alternative=="greater"){
    tstat<-(R1-R2)/sqrt(SD1^2+SD2^2)   
    tcrit2<-(qt(1-alpha,hu)*SD1^2+qt(1-alpha,hu)*SD2^2)/(SD1^2+SD2^2)
    p<-1-pt(tstat*cu,hu)
    hyp<-"mean 1 >= mean 2"
  }
ans<-NULL 
ans$Null<-c(paste("Ho:",hyp))                                   
ans$results<-matrix(NA,1L,6L)
ans$results<-cbind(round(R1,2),round(R2,2),rho,round(nu,1),round(nb,1), round(hu,1),round(tstat,3),round(cu,3),round(tstat*cu,3),round(tcrit2,3),round(p,5))                     
dimnames(ans$results)<-list(cbind(NULL),c("Mean 1","Mean 2","rho","ntilda","nu","df","t","cu","tadj","tcrit","p"))
return(ans) 
}
