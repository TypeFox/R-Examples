clus.rho.g<-function(popchar=NULL, cluster=NULL, group=NULL){
  if(is.null(cluster)) stop("cluster vector not specified.")
  if(is.null(popchar)) stop("popchar vector not specified.")
  if(is.null(group)) stop("popchar vector not specified.")
  if(any(c(length(cluster),length(popchar),length(group)) %in% length(cluster)=="FALSE"))
    stop("vectors lengths are different.")
  if(!is.numeric(popchar)) stop("popchar must be numeric.")
  temp<-NULL
  temp<-as.data.frame(cbind(group,cluster,popchar))
  temp$popchar<-as.numeric(as.character(temp$popchar))
  temp$group<-as.factor(temp$group)
  temp$cluster<-as.factor(temp$cluster)
  names(temp)<-c("group","cluster","popchar")
  temp<-temp[!is.na(temp$group) & !is.na(temp$cluster) & !is.na(temp$popchar),]
  out<-matrix(NA,1,5L)
  dimnames(out)[1]<-list(rep(NA,1))
  colnames(out)<-c("value","Approx.95% LCI","Approx. 95% UCI","F","Prob")
  dimnames(out)[[1]][1]<-list("ANOVA rho")
 #m adjusted
  M<-length(temp$popchar)
  mi<-aggregate(temp$popchar,list(temp$group,temp$cluster),length)
  names(mi)<-c("group","cluster","m")
  mi$m2<-mi$m^2
  summ2<-aggregate(mi$m2,list(mi$group),sum)
  names(summ2)<-c("group","sm2")
  summ<-aggregate(mi$m,list(mi$group),sum)
  names(summ)<-c("group","sm")
  n<-aggregate(mi$cluster,list(mi$group),length)  
  ratio<-merge(summ2,summ,by.x="group",by.y="group",all.x=TRUE,all.y=TRUE)
  den<-sum(n$x-1)
  m0<-(M-sum(ratio$sm2/ratio$sm))/den
  dd<-aov(popchar~group/cluster,data=temp)
  msbna<-summary(dd)[[1]][2,3]
  msw<-summary(dd)[[1]][3,3]
  rho<-(msbna-msw)/(msbna+(m0-1)*msw)
  out[1,1]<-rho
  Zfg<-0.5*log((1+(m0-1)*rho)/(1-rho))
  s2zg<-0.5*((1/sum(n$x-1))+(1/(M-sum(n$x))))
  IL<-Zfg-qnorm(0.975)*sqrt(s2zg)
  IU<-Zfg+qnorm(0.975)*sqrt(s2zg)
  IZ<-(exp(2*Zfg)-1)/(exp(2*Zfg)+m0-1)
  out[1,2]<-(exp(2*IL)-1)/(exp(2*IL)+m0-1)
  out[1,3]<-(exp(2*IU)-1)/(exp(2*IU)+m0-1)
  out[1,4]<-msbna/msw
  ngroups<-nlevels(temp$group)
  df1<-sum(n$x)-ngroups
  df2<-M-sum(n$x)
  out[1,5]<-(1-pf(msbna/msw,df1,df2))  

  #Generate Information for Tests
  msws<-list(NULL)
  thom<-data.frame(Zfv=NA,Wh=NA,nih=NA,Zh=NA)
  tree<-NULL
  gr<-as.character(unique(temp$group))
  for(gg in 1:ngroups){
  tree<-temp[temp$group==gr[gg],]
  dd<-lm(popchar~cluster,data=tree)
  msws[[gg]]<-dd
  mc<-length(tree$popchar)-(sum(aggregate(tree$popchar,list(tree$cluster),length)$x^2)/length(tree$popchar)) 
  mc<-mc/(length(unique(tree$cluster))-1)
  rhoh<-(anova(dd)[1,3]-anova(dd)[2,3])/(anova(dd)[1,3]+(mc-1)*anova(dd)[2,3])
  thom[gg,1]<-0.5*log((1+(mc-1)*rhoh)/(1-rhoh))
  thom[gg,2]<-2/((1/(length(unique(tree$cluster))-1))+(1/(length(tree$popchar)-length(unique(tree$cluster)))))
  thom[gg,3]<-mc
  thom[gg,4]<-0.5*log((1+(mc-1)*rho)/(1-rho))
  }
  
  #BArtlett's Test of Homogeneity of Variances
  btest<-bartlett.test(msws)
  if(btest$p.value>0.05) ho<-"Accept"
  if(btest$p.value<=0.05) ho<-"Reject"
   btest1<-data.frame(Ksquared=btest$statistic[[1]][1],df=btest$parameter[[1]][1],Prob=btest$p.value)
   rownames(btest1)[1]<-""
    #STest for common rho 
    aZwv<-sum(thom[,1]*thom[,2])/sum(thom[,2])
    TH<-sum(((thom[,1]-aZwv)^2)*thom[,2])
    mZp<-sum(thom[,4]*thom[,2])/sum(thom[,2])
    THp<-TH-sum(((thom[,4]-mZp)^2)*thom[,2])
    cp<-data.frame(Tprime=THp,df=ngroups-1,Prob=1-pchisq(THp,ngroups-1))
    rownames(cp)[1]<-""
    outpt<-NULL
    outpt$icc<-out
    outpt$Bartlett<-btest1
    outpt$common_rho<-cp 
  return(outpt)  
}
