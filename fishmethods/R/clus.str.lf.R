clus.str.lf<-function(group=NULL,strata=NULL, weights=NULL, haul=NULL,len=NULL, number=NULL, binsize=NULL, resamples=100){
     if(is.null(group)) stop("group vector not specified.")
     if(is.null(strata)) stop("strata vector not specified.")
     if(is.null(weights)) stop("weight vector not specified.")
     if(is.null(haul)) stop("Haul vector not specified.")
     if(is.null(len)) stop("Length vector not specified.")
     if(is.null(number)) stop("Number vector not specified.")
     if(is.null(binsize)) stop("Bin size not specified.")
     if(any(c(length(group),length(strata),length(weights),length(haul),length(len)) %in% length(group)=="FALSE"))
         stop("vectors lengths are different.")
     haul<-as.character(paste(strata,"-",haul,sep=""))
     group<-as.character(group)
     strata<-as.character(strata)
     if(binsize>0) lenclass<-as.numeric(ifelse(len>0,trunc(len/binsize)*binsize+binsize/2,0))
     if(binsize==0) lenclass<-len
     dater<-as.data.frame(cbind(group,strata,weights,haul,lenclass,number))  
     names(dater)[5]<-"len" 
     dater$number<-as.numeric(as.character(dater$number))
     n1<-aggregate(dater[,6],list(dater$group,dater$strata,dater$haul),length)
     nstratum<-aggregate(n1[,3],list(n1[,1],n1[,2]),length)
     names(nstratum)<-c("group","strata","nk")
     dater$len<-as.numeric(as.character(dater$len))
     dater<-dater[dater$len>0,]
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
  l1<-sort(unique(dater$len))
    l2<-min(abs(l1[1:length(l1)-1]-l1[2:length(l1)]))
    l2<-ifelse(l2==0,1,l2)
    lendist<-as.data.frame(seq(min(as.numeric(as.character(dater$len))),
              max(as.numeric(as.character(dater$len))),l2))
    names(lendist)<-c("len")
compcount<-0
for (r in 1:rloop){
  for (rr in (r+1):rrloop){
    compcount<-compcount+1
    complab<-c(paste(as.character(unique(dater$group[dater$region==r])),"vs",  #comparison labels
              as.character(unique(dater$group[dater$region==rr]))))
    reg<-dater[dater$region %in% c(r,rr),] 
    reg<-reg[order(reg$group,reg$strata,reg$haul,reg$len),]
   
# Create len and number at length  
    sumnum<-aggregate(reg$number,list(reg$group,reg$strata,reg$weights,reg$len),sum) 
    names(sumnum)<-c("group","strata","weights","len","numlen")
    sumnum$weights<-as.integer(as.character(sumnum$weights))

#Merge number of tows with summarized length data and calculate weight mean len per strata
    comb<-merge(sumnum,nstratum,by.x=c("group","strata"),by.y=c("group","strata"))
    comb$weights<-as.integer(as.character(comb$weights))
    comb$wmlen<-comb$weights*(comb$numlen/comb$nk)

#Sum weighted mean length over group and len category 
    sumlength<-aggregate(comb$wmlen,list(comb$group,comb$len),sum)
    names(sumlength)<-c("group","len","sumlen")
  
# Calculate mean number per stratum
  ax<-aggregate(comb$weight*(comb$numlen/comb$nk),list(comb$group),sum)
  names(ax)<-c("group","ax")

# Merge summed weighted mean length with sum strata areas and calculate stratifed mean length for group 
   both<-merge(sumlength,ax,by.x="group",by.y="group")
   both$prop<-both$sumlen/both$ax

# Add missing lengths and calculate proprortions
    both$region<-as.integer(both$group)
    obf1<-subset(both,both$region==r)
    obf2<-subset(both,both$region==rr)
     # Create length dist based on binsize or original binsize
   
    obf1<-merge(obf1,lendist,by.x="len",by.y="len",all.x=T,all.y=T)
    obf2<-merge(obf2,lendist,by.x="len",by.y="len",all.x=T,all.y=T)
    obf1$prop<-ifelse(is.na(obf1$prop),0,obf1$prop)
    obf2$prop<-ifelse(is.na(obf2$prop),0,obf2$prop)
    obf1<-obf1[order(obf1$len),]
    obf2<-obf2[order(obf2$len),]
    obf1$cumprop<-cumsum(obf1$prop)
    obf2$cumprop<-cumsum(obf2$prop)
   if(r==1 & rr==2) obs<-cbind(obf1[,c(1,7)],obf2[,c(7)])
   if(r==1 & rr>2)  obs<-cbind(obs,obf2[,c(7)])       
 
#Calculate observed KS
    Ds[compcount,1]<-round(max(abs(obf1$cumprop-obf2$cumprop)),3)		#Observed KS statistic
    rm(sumnum,comb,both,sumlength,obf1,obf2,ax,l2)
 nh<-aggregate(reg[,4],list(reg$group,reg$strata,reg$weights),function(x){length(levels(factor(x)))})
    names(nh)<-c("group","strata","weights","nk")
 picklen<-reg[,4:6]
 haul<-aggregate(rep(1,length(picklen$haul)),list(picklen$haul),length)[1]
    names(haul)<-"haul"
 strata<-rep(as.numeric(as.character(nh$strata)),nh$nk)
#Randomization Test
   for(k in 1:resamples){
   ss<-as.data.frame(sample(strata,length(strata),replace=FALSE))
   names(ss)<-"strata"
   cc<-cbind(ss,haul);cc<-cc[order(cc$strata),]
   dd<-merge(cc,picklen,by.x="haul",by.y="haul")
   ee<-merge(dd,nh,by.x="strata",by.y="strata")
   ee<-ee[order(ee$group,ee$strata),]
 # Create len and number at length  
    sumnum<-aggregate(ee$number,list(ee$group,ee$strata,ee$weights,ee$len),sum) 
    names(sumnum)<-c("group","strata","weights","len","numlen")
    sumnum$weights<-as.integer(as.character(sumnum$weights))

#Merge number of tows with summarized length data and calculate weight mean len per strata
    comb<-merge(sumnum,nstratum,by.x=c("group","strata"),by.y=c("group","strata"))
    comb$weights<-as.integer(as.character(comb$weights))
    comb$wmlen<-comb$weights*(comb$numlen/comb$nk)

#Sum weighted mean length over group and len category 
    sumlength<-aggregate(comb$wmlen,list(comb$group,comb$len),sum)
    names(sumlength)<-c("group","len","sumlen")
  
# Calculate mean number per stratum
  ax<-aggregate(comb$weight*(comb$numlen/comb$nk),list(comb$group),sum)
  names(ax)<-c("group","ax")

# Merge summed weighted mean length with sum strata areas and calculate stratifed mean length for group 
   both<-merge(sumlength,ax,by.x="group",by.y="group")
    both$prop<-both$sumlen/both$ax

# Add missing lengths and calculate proprortions
    both$region<-as.integer(both$group)
    obf1<-subset(both,both$region==r)
    obf2<-subset(both,both$region==rr)
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
       rm(sumnum,comb,both,sumlength,obf1,obf2,ax)
     }#close randomization 
  }#rr loop
}#r loop

#calculate significance level
p<-as.data.frame(rep(NA,compcount))
names(p)[1]<-c("p")
for (j in 1:compcount){
   stat<-B[,j]
    p[j,1]<-round((length(stat[stat>=Ds[j,1]]))/(resamples),3)
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
