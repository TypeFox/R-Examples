seerStats<-function(canc,popsa) {
  age=agedx=db=py=reg=cancer=over99=total=sex=year=seqnum=gte100=PY=NULL 
  #   cf=function (x) comma_format()(x)
  #seerStats(canc,popsae)
  L=NULL
  L$PDR=popsa%>%group_by(db,reg)%>%summarize(PY=round(sum(py)/1e6,1),age=weighted.mean(age,py))
  L$PSY=popsa%>%group_by(sex,year)%>%summarize(PY=round(sum(py)/1e6,1))
  D=canc%>%group_by(cancer,db)%>%summarize(n=n())
  A=canc%>%filter(agedx>99.5)%>%group_by(cancer)%>%summarize(n=n())
  canc$seqnum=ifelse(canc$seqnum>3,4,canc$seqnum)
  canc$seqnum=ifelse(canc$seqnum==0,1,canc$seqnum)
  Sq=canc%>%group_by(cancer,seqnum)%>%summarize(n=n())
  C=dcast(D,cancer~db,value.var="n",fun.aggregate = sum,margins=c("db"))
  #   C=dcast(D,cancer~db,value.var="n",fun.aggregate = sum,margins=TRUE)
  O=dcast(A,cancer~.,value.var="n")
  S=dcast(Sq,cancer~seqnum,value.var="n",fun.aggregate = sum,margins=c("seqnum")) #checks fine with total
  #    S=dcast(Sq,cancer~seqnum,value.var="n")
  d=left_join(C,O,by="cancer")
  names(d)[5:6]=c("total","gte100")
  d=d%>%mutate(under100=total-gte100)
  d=left_join(d,S,by="cancer")
  names(d)[8:12]=c("firsts","seconds","thirds","highers","all")
  d[is.na(d)]=0
  d$cancer=as.character(d$cancer)
  d[dim(d)[1]+1,]=data.frame("total",t(sapply(d[-1],sum)))
  d[dim(d)[1],1]="total"
  L$d=d
  cat(paste("\nTable 1. Cases per database, by >= or < 100, and by cancer sequence, through",
            max(popsa$year),"(sexes and races pooled).\n"))
  print(d) 
  cat(paste("\nTable 2. Population PY in millions and PY-weighted mean ages, through",
            max(popsa$year),"(sexes and races pooled).\n"))
  cat(paste(""))
  print(as.data.frame(L$PDR)) 
  N=sum(L$PSY$PY)
  tit=paste(N,"Million Total Person-Years")
  p=qplot(year,PY,col=sex,ylab="Person-Years (Millions)",main=tit,data=L$PSY)+
    theme(legend.position = c(.18, .75),legend.title = element_blank(),legend.text=element_text(size=rel(1.5)),
          axis.text=element_text(size=rel(2)),axis.title=element_text(size=rel(2)))
  print(p)
  
  invisible(L)
} 
