summary.seerSet<-function(object, ...) {
  age=trt=race=surv=year=py=popsa=cancer=NULL 
  cf=function (x) comma_format()(x)
  D=object$canc%>%filter(trt!="unk")%>%group_by(cancer,trt)%>%
    summarize(n=n(),age=round(mean(age),1),surv=round(median(surv),1)) #,
#               seq=mean(ifelse(seqnum==0,1,seqnum)) ) #%>%filter(n>9)
  P=object$popsa%>%group_by(year)%>%summarize(PY=round(sum(py)/1e6,1))
  A=dcast(D,cancer~trt,value.var="age")
  S=dcast(D,cancer~trt,value.var="surv")
  #   Sq=dcast(D,cancer~trt,value.var="seq")
  N=dcast(D,cancer~trt,value.var="n")
  #   N=dcast(D,cancer~trt,value.var="n",margins=c("cancer"),fun.agg=sum)
  d=left_join(N,A,by="cancer")
  d=left_join(d,S,by="cancer")
  names(d)=c("Cancer",paste(rep(c("Count","Age","Survival"),each=2),c("rad","noRad"),sep="."))
  seerSetSum=NULL
  seerSetSum$title=paste0("              Counts, Means of Ages, and Median Survivals in Years\n               Sex: ",
                          object$sex,"    Race: ",object$race,
                          "   Years: ",min(object$popsa$year),"-",max(object$popsa$year) ,"\n")
  Cnts=c(total=dim(object$canc)[1],
         unkTrt=dim(object$canc%>%filter(trt=="unk"))[1],
         unkTrtNsurv=dim(object$canc%>%filter(trt=="unk",surv>100))[1],
         unkSurv=dim(object$canc%>%filter(surv>100))[1])
  seerSetSum$cnts=Cnts
  seerSetSum$sex=object$sex
  
  seerSetSum$notes=c(paste("Of",cf(Cnts["total"]),"total",object$sex,"cases of",object$race,"race,",cf(Cnts["unkTrt"]),"with unknown treatment were not included."),
                     paste("Of",cf(Cnts["unkSurv"]),"cases with unknown survival,",cf(Cnts["unkTrtNsurv"]),
                           "were excluded due to also having unknown treatment."),
                     "In 2005, due to hurricane Katrina some PY (and cases) are kept in a separate database not used here.")
  
  seerSetSum$P=P
  seerSetSum$d=d
  class(seerSetSum)="seerSet.summary"
  seerSetSum
} 
