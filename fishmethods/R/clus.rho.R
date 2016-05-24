clus.rho<-function(popchar=NULL, cluster=NULL, type=c(1,2,3),est=0,nboot=500){
  if(is.null(cluster)) stop("cluster vector not specified.")
  if(is.null(popchar)) stop("popchar vector not specified.")
  if(any(c(length(cluster),length(popchar)) %in% length(cluster)=="FALSE"))
    stop("vectors lengths are different.")
  if(!is.numeric(popchar)) stop("popchar must be numeric.")
  temp<-NULL
  temp<-as.data.frame(cbind(cluster,popchar))
  temp$popchar<-as.numeric(as.character(temp$popchar))
  names(temp)<-c("cluster","popchar")
  temp<-temp[!is.na(temp$cluster) & !is.na(temp$popchar),]
 
  if(est==0){
  n<-length(type)
  out<-matrix(NA,n,1L)
  dimnames(out)<-list(rep(NA,n),c("value"))
  rowcnt<-0
  if(any(type==1)){
    # ICC
    rowcnt<-rowcnt+1
    dimnames(out)[[1]][rowcnt]<-list("Lohr rho")
    ssto<-sum((temp$popchar-mean(temp$popchar))^2)
    myui<-aggregate(temp$popchar,list(temp$cluster),mean)
    names(myui)<-c("cluster","mean")
    newdats<-merge(temp,myui,by.x="cluster",by.y="cluster")
    ssw<-sum((newdats$popchar-newdats$mean)^2)
    catch<-aggregate(temp$popchar,list(temp$cluster),length)
    avgcatch<-mean(catch[,2])
    out[rowcnt,1]<-1-avgcatch/(avgcatch-1)*(ssw/ssto)
  }
  if(any(type==2)){
    # adjusted r2
    rowcnt<-rowcnt+1
    dimnames(out)[[1]][rowcnt]<-list("Adjusted r-square")
    ssto<-sum((temp$popchar-mean(temp$popchar))^2)
    myui<-aggregate(temp$popchar,list(temp$cluster),mean)
    names(myui)<-c("cluster","mean")
    newdats<-merge(temp,myui,by.x="cluster",by.y="cluster")
    ssw<-sum((newdats$popchar-newdats$mean)^2)
    mm<-aggregate(temp$popchar,list(temp$cluster),length)
    sm<-sum(mm[,2]-1)
    msw<-ssw/sm
    S2<-ssto/(sum(mm[,2])-1)
    out[rowcnt,1]<-1-msw/S2
  }
  if(any(type==3)){
    #ANOVA
    rowcnt<-rowcnt+1
    dimnames(out)[[1]][rowcnt]<-list("ANOVA rho")
    mi<-aggregate(temp$popchar,list(temp$cluster),length)
    names(mi)<-c("cluster","m")
    yi.<-aggregate(temp$popchar,list(temp$cluster),mean)
    names(yi.)<-c("cluster","yi")
    sumtab<-merge(yi.,mi,by.x="cluster",by.y="cluster",all.x=TRUE,all.y=TRUE)
    y..<-mean(temp$popchar)
    M<-sum(mi$m)
    ssg<-sum(sumtab$m*sumtab$yi^2)-M*y..^2
    sse<-sum(temp$popchar^2)-sum(sumtab$m*sumtab$yi^2)
    MSE<-sse/(M-length(unique(temp$cluster)))
    MSG<-ssg/(length(unique(temp$cluster))-1)
    mprime<-(sum(mi[,2])-(sum(mi[,2]^2)/sum(mi[,2])))*(1/(length(mi[,1])-1))
    out[rowcnt,1]<-(MSG-MSE)/(MSG+(mprime-1)*MSE)    
  }
  outpt<-NULL
  outpt$icc<-out
  return(outpt)
  } #est=0
  
 #### Estimate rho variance
  if(est==1){
   n<-length(type)
   out<-matrix(NA,n,4L)
   dimnames(out)<-list(rep(NA,n))
   colnames(out)<-c("Value","Var","P2.5%","P97.5%")
   # temp storage
    store<-as.data.frame(rbind(1:n))  
  for(b in 1:nboot){ # bootstrap 
     temp1<-data.frame(cluster=NA,popchar=NA)
  # select clusters
    clusters<-unique(temp$cluster)
    c1<-sample(clusters,length(clusters),replace=TRUE)
    for(k in 1:as.numeric(length(c1))){
      dodo<-temp[temp$cluster==c1[k],]
      dodo$cluster<-k
      temp1<-rbind(temp1,dodo) 
    }
    temp1<-temp1[!is.na(temp1$cluster),]
    colcnt<-0
    if(any(type==1)){
    # ICC
     colcnt<-colcnt+1
     ssto<-sum((temp1$popchar-mean(temp1$popchar))^2)
     myui<-aggregate(temp1$popchar,list(temp1$cluster),mean)
     names(myui)<-c("cluster","mean")
     newdats<-merge(temp1,myui,by.x="cluster",by.y="cluster")
     ssw<-sum((newdats$popchar-newdats$mean)^2)
     catch<-aggregate(temp1$popchar,list(temp1$cluster),length)
     avgcatch<-mean(catch[,2])
     store[b,colcnt]<-1-avgcatch/(avgcatch-1)*(ssw/ssto)
   }
   if(any(type==2)){
    # adjusted r2
    colcnt<-colcnt+1
    ssto<-sum((temp1$popchar-mean(temp1$popchar))^2)
    myui<-aggregate(temp1$popchar,list(temp1$cluster),mean)
    names(myui)<-c("cluster","mean")
    newdats<-merge(temp1,myui,by.x="cluster",by.y="cluster")
    ssw<-sum((newdats$popchar-newdats$mean)^2)
    mm<-aggregate(temp1$popchar,list(temp1$cluster),length)
    sm<-sum(mm[,2]-1)
    msw<-ssw/sm
    S2<-ssto/(sum(mm[,2])-1)
    store[b,colcnt]<-1-msw/S2
   }
   if(any(type==3)){
    #ANOVA
    colcnt<-colcnt+1
    mi<-aggregate(temp1$popchar,list(temp1$cluster),length)
    names(mi)<-c("cluster","m")
    yi.<-aggregate(temp1$popchar,list(temp1$cluster),mean)
    names(yi.)<-c("cluster","yi")
    sumtab<-merge(yi.,mi,by.x="cluster",by.y="cluster",all.x=TRUE,all.y=TRUE)
    y..<-mean(temp1$popchar)
    M<-sum(mi$m)
    ssg<-sum(sumtab$m*sumtab$yi^2)-M*y..^2
    sse<-sum(temp1$popchar^2)-sum(sumtab$m*sumtab$yi^2)
    MSE<-sse/(M-length(unique(temp$cluster)))
    MSG<-ssg/(length(unique(temp$cluster))-1)
    mprime<-(sum(mi[,2])-(sum(mi[,2]^2)/sum(mi[,2])))*(1/(length(mi[,1])-1))
    store[b,colcnt]<-(MSG-MSE)/(MSG+(mprime-1)*MSE)    
   }
  }# boot  
  rowcnt<-0
  colcnt<-0
  if(any(type==1)){
    rowcnt<-rowcnt+1
    colcnt<-colcnt+1
    dimnames(out)[[1]][rowcnt]<-list("Lohr rho")
    out[rowcnt,1]<-mean(store[,colcnt])
    out[rowcnt,2]<-var(store[,colcnt])
    out[rowcnt,3]<-quantile(store[,colcnt],probs=0.025)[[1]]
    out[rowcnt,4]<-quantile(store[,colcnt],probs=0.975)[[1]] 
    names(store)[colcnt]<-"Lohr"
  } 
  if(any(type==2)){
    rowcnt<-rowcnt+1
     colcnt<-colcnt+1
    dimnames(out)[[1]][rowcnt]<-list("Adjusted r-square")
    out[rowcnt,1]<-mean(store[,colcnt])
    out[rowcnt,2]<-var(store[,colcnt])
    out[rowcnt,3]<-quantile(store[,colcnt],probs=0.025)[[1]]
    out[rowcnt,4]<-quantile(store[,colcnt],probs=0.975)[[1]]
    names(store)[colcnt]<-"r-square"
   } 
  if(any(type==3)){
   rowcnt<-rowcnt+1
   colcnt<-colcnt+1
   dimnames(out)[[1]][rowcnt]<-list("ANOVA rho")
   out[rowcnt,1]<-mean(store[,colcnt])
   out[rowcnt,2]<-var(store[,colcnt])
   out[rowcnt,3]<-quantile(store[,colcnt],probs=0.025)[[1]]
   out[rowcnt,4]<-quantile(store[,colcnt],probs=0.975)[[1]]
  names(store)[colcnt]<-"ANOVA"
  } 
  outpt<-NULL
  outpt$icc<-out
  outpt$bootstraps<-store
  return(outpt)
  } #est=1
}
