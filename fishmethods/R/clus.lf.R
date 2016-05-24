clus.lf<-function(group=NULL,haul=NULL,len=NULL, number=NULL, binsize=NULL, resamples=100){
     if(is.null(group)) stop("Comparison vector not specified.")
     if(is.null(haul)) stop("Haul vector not specified.")
     if(is.null(len)) stop("Length vector not specified.")
     if(is.null(number)) stop("Number vector not specified.")
     if(is.null(binsize)) stop("Bin size not specified.")
     if(any(c(length(group),length(haul),length(len),length(number)) %in% length(group)=="FALSE"))
         stop("vectors lengths are different.")
     haul<-as.character(haul)
     group<-as.character(group)
     if(binsize>0) lenclass<-as.numeric(ifelse(len>0,trunc(len/binsize)*binsize+binsize/2,0))
     if(binsize==0) lenclass<-len
     dater<-as.data.frame(cbind(group,haul,lenclass,number))
     dater$number<-as.numeric(as.character(dater$number))
     dater<-dater[dater$number>0,]
     names(dater)<-c("group","haul","len","number")
     dater$len<-as.numeric(as.character(dater$len))
     dater$group<-as.factor(dater$group)
     dater$region<-as.integer(dater$group)

  # count levels
  nl<-nlevels(dater$group)
  numcol<-factorial(nl)/(factorial(nl-2)*factorial(2))
  #looping variables
  rloop<-nl-1
  rrloop<-nl
  #create ouput dataframes
  B<-as.data.frame(array(rep(NA,resamples),dim=c(resamples,numcol)))
  Ds<-as.data.frame(rep(NA,numcol))
  names(Ds)[1]<-c("Ds")
  #Create length frequency master
    l1<-sort(unique(dater$len))
    l2<-min(abs(l1[1:length(l1)-1]-l1[2:length(l1)]))
    l2<-ifelse(l2==0,1,l2)
    lendist<-as.data.frame(seq(min(as.numeric(as.character(dater$len))),
              max(as.numeric(as.character(dater$len))),l2))
    names(lendist)<-c("len")
  
compcount<-0;r<-1;rr<-2
for (r in 1:rloop)
  {
  for (rr in (r+1):rrloop)
    {
    compcount<-compcount+1
    complab<-c(paste(as.character(unique(dater$group[dater$region==r])),"vs",  #comparison labels
              as.character(unique(dater$group[dater$region==rr]))))
    reg<-dater[dater$region %in% c(r,rr),] 
    reg<-reg[order(reg$group,reg$haul,reg$len),]
     #number hauls by region
    nhreg<-aggregate(reg[,2],list(reg$group),function(x){length(levels(factor(x)))})
    reg$count<-1
    for(i in 1:length(reg$number))					  #label haul number
      {
	ifelse(reg$haul[i]!=reg$haul[i+1],reg$count[i+1]<-reg$count[i]+1,
             reg$count[i+1]<-reg$count[i])
      }

#Calculate observed cumulative frequency distributons and KS
  
    obf1<-subset(reg,reg$region==r)
    obf1<-aggregate(obf1$number,list(obf1$len),sum)
    obf1$prop1<-(obf1$x/sum(obf1$x))
    obf1<-obf1[,c(1,3)]
    names(obf1)<-c("len","prop")
    
    obf2<-subset(reg,reg$region==rr)
    obf2<-aggregate(obf2$number,list(obf2$len),sum)
    obf2$prop2<-(obf2$x/sum(obf2$x))
    obf2<-obf2[,c(1,3)]
    names(obf2)<-c("len","prop")

    obf1<-merge(obf1,lendist,by.x="len",by.y="len",all.x=T,all.y=T)
    obf2<-merge(obf2,lendist,by.x="len",by.y="len",all.x=T,all.y=T)

    obf1$prop<-ifelse(is.na(obf1$prop),0,obf1$prop)
    obf2$prop<-ifelse(is.na(obf2$prop),0,obf2$prop)
    obf1<-obf1[order(obf1$len),]
    obf2<-obf2[order(obf2$len),]
    obf1$cumprop<-cumsum(obf1$prop)
    obf2$cumprop<-cumsum(obf2$prop)
    if(r==1 & rr==2) obs<-cbind(obf1[,c(1,3)],obf2[,c(3)])
    if(r==1 & rr>2)  obs<-cbind(obs,obf2[,c(3)])
 
    Ds[compcount,1]<-round(max(abs(obf1$cumprop-obf2$cumprop)),3)		#Observed KS statistic
    rm(obf1,obf2)

    tnh<-sum(nhreg$x)
#Randomization Test
   for(k in 1:resamples)
    {
       ss<-cbind(sample(1:tnh,nhreg[1,2],replace=FALSE),1)
       select<-merge(reg,ss,by.x=c("count"),by.y=ss[,2],all.x=TRUE)
       select$V2<-ifelse(is.na(select$V2),2,select$V2)

#segregate regions, calculate bin class, estimate proportions
    dodo<-aggregate(select$number,list(select$V2,select$len),sum)
    obf1<-subset(dodo,dodo$Group.1==1)
    obf1$prop1<-(obf1$x/sum(obf1$x))
    obf1<-obf1[,c(2,4)]
    names(obf1)<-c("len","prop")
    
    obf2<-subset(dodo,dodo$Group.1==2)
    obf2$prop2<-(obf2$x/sum(obf2$x))
    obf2<-obf2[,c(2,4)]
    names(obf2)<-c("len","prop")

    obf1<-merge(obf1,lendist,by.x="len",by.y="len",all.x=T,all.y=T)
    obf2<-merge(obf2,lendist,by.x="len",by.y="len",all.x=T,all.y=T)

    obf1$prop<-ifelse(is.na(obf1$prop),0,obf1$prop)
    obf2$prop<-ifelse(is.na(obf2$prop),0,obf2$prop)
    obf1<-obf1[order(obf1$len),]
    obf2<-obf2[order(obf2$len),]
    obf1$cumprop<-cumsum(obf1$prop)
    obf2$cumprop<-cumsum(obf2$prop)
#calculate maximum absolute difference 2-tailed KS
       KS<-round(max(abs(obf1$cumprop-obf2$cumprop)),3)
	 B[k,compcount]<-KS
	 names(B)[compcount]<-complab
     }#close randomization 
  }#rr loop
}#r loop
#calculate significance level
p<-as.data.frame(rep(NA,compcount))
names(p)[1]<-c("p")
for (j in 1:compcount){
   stat<-(B[,j])
    p[j,1]<-round(length(stat[stat>=Ds[j,1]])/resamples,2)
  }
ans<-NULL;ans$results<-matrix(NA,length(Ds[,1]),2L)
ans$results<-rbind(cbind(Ds,p)) 
dimnames(ans$results)<-list(cbind(names(B)),c("Ds","p"))
labels<-aggregate(dater$number,list(dater$group,dater$region),length)
names(obs)<-c("len",as.character(labels$Group.1))
ans$obs_prop<-obs
ans$Drandom<-B 
return(ans)
}
