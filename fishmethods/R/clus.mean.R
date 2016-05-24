clus.mean<-function(popchar=NULL, cluster=NULL, clustotal=NULL,rho=NULL,nboot=1000){
  if(is.null(cluster)) stop("cluster vector not specified.")
  if(is.null(popchar)) stop("popchar vector not specified.")
  if(is.null(clustotal)) stop("clustotal vector not specified.")                                  
  if(any(c(length(popchar),length(cluster),length(clustotal)) %in% length(popchar)=="FALSE"))
    stop("vectors lengths are different.")                                  
     x<-as.data.frame(cbind(popchar,cluster,clustotal)) 
     x$popchar<-as.numeric(as.character(x$popchar))
     x$clustotal<-as.numeric(as.character(x$clustotal))
     x<-x[!is.na(x$popchar) & !is.na(x$cluster) & !is.na(x$clustotal),]
     x<-x[x$popchar>0,]
     means<-aggregate(x$popchar,list(x$cluster),mean)
     names(means)<-c("cluster","mean")
     sumM<-aggregate(x$popchar,list(x$cluster,x$clustotal),length)
     names(sumM)<-c("cluster","M","m")
     stats<-merge(means,sumM,by.x=c("cluster"),by.y=c("cluster"))
     n<-length(stats$cluster)
     R<-sum(stats$mean*stats$M)/sum(stats$M)
     #Variance usual approximation
     varU<-sum(((stats$mean-R)^2*(stats$M/mean(stats$M))^2)/(n*(n-1)))
     #Variance Jackknife
      stats$wmean<-stats$mean*stats$M
      outs2<-jackknife(1:length(stats[,1]),function(x,stats){sum(stats[x,5])/sum(stats[x,3])},stats)
      varJ<-outs2$jack.se^2  
     #Bootstrap var
      ratio <- function(d, w) sum(d$wmean * w)/sum(d$M * w)
      outs<-boot(stats, ratio, R = nboot, stype = "w")
      varB<-apply(outs$t,2,var)  
     # Remaining terms
     rss<-aggregate((popchar-R)^2,list(cluster),sum)
     names(rss)<-c("cluster","rss")
     stats2<-merge(stats,rss,by.x="cluster",by.y="cluster")
     s2x<-sum((stats2$M/stats2$m)*stats2$rss)/(sum(stats2$M)-1)
     M<-sum(stats2$M);m<-sum(stats2$m);meffU<-round(s2x/varU,1)
     meffJ<-round(s2x/varJ,1);meffB<-round(s2x/varB,1)
     # Degrees of freedom based on Hedges 2007
    if(!is.null(rho)){
    M<-sum(stats2$m)
    mu<-sum(stats2$m^2)/sum(stats2$m)
    A<-((M^2)*sum(stats2$m^2)+(sum(stats2$m^2)^2)-2*M*sum(stats2$m^3))/M^2
    num<-((M-1)-2*(mu-1)*rho)^2
    den<-(M-1)*((1-rho)^2)+A*(rho^2)+2*(M-2*mu)*rho*(1-rho)
    df<-round(num/den,0)
    }
    if(is.null(rho)){
     ans<-matrix(NA,1L,7L)
     ans<-cbind(n,M,m,R,varU,varJ,varB,s2x,meffU,meffJ,meffB)
     dimnames(ans)[[1]][1]<-list(" ")
    }
    if(!is.null(rho)){
     ans<-matrix(NA,1L,8L)
     ans<-cbind(n,M,m,R,varU,varJ,varB,s2x,meffU,meffJ,meffB,df)
     dimnames(ans)[[1]][1]<-list(" ")
    }
     return(ans) 
}
