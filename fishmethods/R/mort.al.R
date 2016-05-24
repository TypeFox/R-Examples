mort.al<-function(relyr=NULL,tal=NULL,N=NULL, method=c(1,2,3),np=0,stper=NULL,nboot=500){
  if(is.null(tal)) stop("No times-at-large data")
  if(is.null(N)) stop("N data is missing")
  if(length(relyr)!=length(tal)) stop("Length of year and tal differ")
  if(length(N)!=length(tal)) stop("Length of N and tal differ")
  if(length(N)!=length(relyr)) stop("Length of N and tal differ")
  if(any(is.na(tal))) stop("Missing atlarge values are not allowed.")
  if(any(is.na(N))) stop("Missing N values are not allowed.")
  if(any(is.na(relyr))) stop("Missing relyr values are not allowed.")
  if(is.null(np)|is.na(np)) stop("np must be numeric")
  #Default setting
  datar<-as.data.frame(cbind(relyr,tal,N))
  datar$period<-"NA"
  if(np==0) datar$period<-as.character(datar$relyr)
  if(np>0){
    if(length(stper)!=np) stop("The number of speriods does not match np")
    if(length(stper)==np){
      if(min(datar$relyr)!=stper[1]) stop("The first year in the dataset must be the first year in speriods ")
      styrs<-c(stper,max(datar$relyr))
      for(y in 1:np){
       if(y<np){
         for(l in 1:c(length(datar[,1]))) datar$period[l]<-ifelse(datar$relyr[l]>=styrs[y] & datar$relyr[l]<styrs[y+1],
         c(paste(styrs[y],"-",c(styrs[y+1]-1),sep="")),datar$period[l])
       }
       if(y==np){
         for(l in 1:c(length(datar[,1]))) datar$period[l]<-ifelse(datar$relyr[l]>=styrs[y] & datar$relyr[l]<=styrs[y+1],
           c(paste(styrs[y],"-",styrs[y+1],sep="")),datar$period[l]) 
       }  
      }
    }
  }
   
  out<-data.frame(Method=NA,Period=NA,M=NA,MSE=NA,F=NA,FSE=NA,Z=NA,ZSE=NA)
  cnt<-0
  if(any(method==1)){
    cnt<-cnt+1
    getstats<-function(x){
      c(((length(x)-1)/length(x))/mean(x),
       ((((length(x)-1)/length(x))/mean(x))/sqrt(length(x)-2))^2) 
    }
    mcg<-aggregate(datar$tal,list(datar$period,datar$relyr),getstats)
  
    mcg<-cbind(mcg[,1],mcg[,2],as.data.frame(mcg[,3]))
    names(mcg)<-c("period","relyr","Z","VAR")
    mcg<-aggregate(cbind(mcg$Z,mcg$VAR),list(mcg$period),mean,na.rm=T)
    out[cnt:c(length(mcg[,1])),1]<-"McGarvey"
    out[cnt:c(length(mcg[,1])),2]<-as.character(mcg[,1])
    out[cnt:c(length(mcg[,1])),7]<-mcg$V1
    out[cnt:c(length(mcg[,1])),8]<-sqrt(mcg$V2) 
    cnt<-c(length(mcg[,1]))    
  }
  if(any(method==2)){
    cnt<-cnt+1
    bigN<-aggregate(datar$tal,list(datar$period,datar$relyr,datar$N),length)
    names(bigN)<-c("period","relyr","N","n")
    bigN$N<-as.numeric(as.character(bigN$N))
    N1<-aggregate(cbind(bigN$N,bigN$n),list(bigN$period,bigN$relyr),sum)
    names(N1)<-c("period","relyr","N","n")  
    npers<-length(unique(N1$period))
    pers<-unique(N1$period)      
    st<-aggregate(datar$tal,list(datar$period,datar$relyr),sum,na.rm=T)
    st$FF<-(N1$n)^2/(N1$N*st$x)
    st$M<-((N1$N-N1$n)*N1$n)/(N1$N*st$x)
    st$Z<-st$FF+st$M
    stats<-NULL
for(t in 1:nboot){
  temp3<-NULL
  for(y in 1:c(length(N1[,1]))){
  ntem<-data.frame(n=N1[y,"n"],N=N1[y,"N"]) 
  newn<-rbinom(1,size=ntem$N,prob=ntem$n/ntem$N)
  newn<-ifelse(newn==0,1,newn)  
  if(newn>0) {
    storep<-data.frame(tal=rexp(newn,rate=st$Z[y]),period=st[y,1],relyr=st[y,2])
    storep$tal<-ifelse(storep$tal<=0,min(datar$tal),storep$tal)
  }
  if(newn==0) storep<-data.frame(tal=0,period=st[y,1],relyr=st[y,2])
  temp3<-rbind(temp3,storep)
  }
bootst<-aggregate(temp3$tal,list(temp3$period,temp3$relyr),sum,na.rm=T)
 bootn<-aggregate(temp3$tal,list(temp3$period,temp3$relyr),length)
 bootst$FF<-(bootn$x)^2/(N1$N*bootst$x)
 bootst$M<-((N1$N-bootn$x)*bootn$x)/(N1$N*bootst$x)
 bootst$Z<-bootst$FF+bootst$M 
 
 stats<-rbind(stats,aggregate(cbind(bootst$FF,bootst$M,bootst$Z),list(bootst[,1],bootst[,2]),mean,na.rm=T))
}

bootse<-aggregate(cbind(stats[,3],stats[,4],stats[,5]),list(stats[,1],stats[,2]),sd,na.rm=T)
bootse<-aggregate(cbind(bootse[,3]^2,bootse[,4]^2,bootse[,5]^2),list(as.character(bootse[,1])),mean,na.rm=T)
bootse[,2:4]<-sqrt(bootse[,2:4])    
means<-aggregate(cbind(st$FF,st$M,st$Z),list(st[,1]),mean,na.rm=T)
out[cnt:c(npers+cnt-1),1]<-"Gulland"
out[cnt:c(npers+cnt-1),2]<-as.character(bootse[,1])
out[cnt:c(npers+cnt-1),3]<-means[,3]
out[cnt:c(npers+cnt-1),4]<-bootse[,3]
out[cnt:c(npers+cnt-1),5]<-means[,2]
out[cnt:c(npers+cnt-1),6]<-bootse[,2]
out[cnt:c(npers+cnt-1),7]<-means[,4]
out[cnt:c(npers+cnt-1),8]<-bootse[,4]
cnt<-(npers+cnt-1)   
}

  return(out)
}
 