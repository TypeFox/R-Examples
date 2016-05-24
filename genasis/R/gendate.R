gendate<-function (x) {
  if ((length(unique(class(x))>1))|(unique(class(x))[1]!="Date")) {
    x<-as.character(substr(x,1,10))
    if (max(grepl("\\.",x),na.rm=TRUE)==1) {
      for (i in 1:length(x)) {
        y<-unlist(strsplit(x[i],"\\."))
        z<-rep(NA,3)
        z[1]<-paste0(y[3],"-"); z[2]<-paste0(y[2],"-"); z[3]<-y[1]
        x[i]<-z[3]
        if (z[2]!="NA-") {x[i]<-paste0(z[2],x[i])}
        if (z[1]!="NA-") {x[i]<-paste0(z[1],x[i])}
      }
    }
    x[which(nchar(x)==4)]<-paste0(x[which(nchar(x)==4)],"-07-01")
    x[which(nchar(x)==6)]<-paste0(x[which(nchar(x)==6)],"-15")
    x[which(nchar(x)==7)]<-paste0(x[which(nchar(x)==7)],"-15")
    for (j in 1:length(x)) {
      y<-as.numeric(unlist(strsplit(as.character(x[j]),"-")))
      y[which(nchar(y)==1)]<-paste0("0",y[which(nchar(y)==1)])
      x[j]<-paste0(y[1],"-",y[2],"-",y[3])
    }
    x<-as.Date(x)
  }
  res<-x
  return(res)
}