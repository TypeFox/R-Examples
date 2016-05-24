alk<-function(age=NULL,size=NULL,binsize=NULL,type=1){
  if(is.null(age)) stop ("age vector does not exist")
  if(!is.numeric(age)) stop ("age vector is not numeric")
  if(is.null(size)) stop ("size vector does not exist")
  if(!is.numeric(size)) stop ("size vector is not numeric")
  d<-NULL;outs<-NULL;ll<-NULL;ul<-NULL;la<-NULL;ua<-NULL;lenlist<-NULL
  agelist<-NULL
  d<-as.data.frame(cbind(age,size))
  d<-d[!is.na(d$size),]
  if(binsize>0) d$lenclass<-trunc(d$size/binsize)*binsize+binsize/2
  if(binsize==0) d$lenclass<-d$size
  ll<-min(d$lenclass)
  ul<-max(d$lenclass)
  la<-min(d$age,na.rm=TRUE)
  ua<-max(d$age,na.rm=TRUE)
  if(binsize>0) lenlist<-seq(ll,ul,by=binsize)
  if(binsize==0) {
     cc<-sort(unique(d$lenclass))
     ints<-NULL
     for(i in 1:as.numeric(length(cc)-1)){
        ints[i]<-cc[i+1]-cc[i]
     }
     lenlist<-seq(ll,ul,by=min(ints))
  }

  agelist<-seq(la,ua,by=1)
  agelen<-matrix(0,nrow=length(lenlist),ncol=length(agelist))
  ff<-d[!is.na(d$age),]
  for(loops in 1:as.numeric(length(ff$age))){
    for(a in 1:length(agelist)){
      for(b in 1:length(lenlist)){
        if(ff$age[loops]==agelist[a] & ff$lenclass[loops]==lenlist[b]){
          agelen[b,a]<-agelen[b,a]+1
        }
      }
    }
  }
  if(type==1){
     if(binsize>0) lens<-trunc(d$size/binsize)*binsize+binsize/2
     if(binsize==0) lens<-d$size
     lenlist<-as.data.frame(lenlist)
     names(lenlist)<-"len"
     lenfreq<-as.data.frame(table(lens))
     names(lenfreq)<-c("len","number")
     final<-merge(lenlist,lenfreq,by.x="len",by.y="len",all.x=TRUE,all.y=TRUE)	
     final[,2]<-ifelse(is.na(final[,2]),0,final[,2])
     outs<-cbind(lenlist$len,final[,2],as.data.frame(agelen))
     names(outs)<-c("len","nl",paste("A",agelist,sep=""))
  }
  if(type==2){
    if(binsize>0) lens<-trunc(d$size/binsize)*binsize+binsize/2
     if(binsize==0) lens<-d$size
     lenlist<-as.data.frame(lenlist)
     names(lenlist)<-"len"
     lenfreq<-as.data.frame(table(lens))
     names(lenfreq)<-c("len","number")
     final<-merge(lenlist,lenfreq,by.x="len",by.y="len",all.x=TRUE,all.y=TRUE)  
     final[,2]<-ifelse(is.na(final[,2]),0,final[,2])
     outs<-cbind(lenlist$len,final[,2],as.data.frame(agelen)/rowSums(agelen))
      names(outs)<-c("len","nl",paste("A",agelist,sep=""))
    for(i in 3:ncol(outs)){
      for(j in 1:nrow(outs)){
        if(is.nan(outs[j,i])) outs[j,i]<-0
      }
    }
  }  
  return(outs)
}

