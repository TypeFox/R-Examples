# additional function
weatherSeverity <-
function(weather,severity,dates,EmergDate,EndEpidDate,NoReadingsH,RHthreshold)
{
nday<- as.numeric(dates-EmergDate)
#........
MeanSeverity<-apply(severity[,4:8],2,function(x) mean(x,na.rm=TRUE))
StDevSeverity<-apply(severity[,4:8],2,function(x)sd(x,na.rm=TRUE))
Sfile<-data.frame(severity[,1:2],dates,nday,MeanSeverity,StDevSeverity)
# ----------------------------
Wfile<-subset(weather,((weather[,1] >= EmergDate) &  (weather[,1] <= EndEpidDate)))
Wfile <- data.frame(Wfile,"","","")
nc<-ncol(Wfile)
for (j in 2:nc) Wfile[,j]<- as.character(Wfile[,j])
for(j in 2:nc) Wfile[Wfile[,j]==".",j]<-""
nr <-nrow(Wfile)
Wfile[,3]<-as.numeric(Wfile[,3])
Wfile[,4]<-as.numeric(Wfile[,4])
Wfile[,5]<-as.numeric(Wfile[,5])
Wfile[,6]<-as.numeric(Wfile[,6])
Wfile[,7]<-as.numeric(Wfile[,7])
Wfile[,8]<-as.numeric(Wfile[,8])

for (i in 1:nr) {
hours<-strsplit(Wfile[i,2], ":", fixed = TRUE)[[1]][1]
if (Wfile[i,2]!="") Wfile[i,6] <- as.numeric(hours)
if (!is.na(Wfile[i,4])){
if (Wfile[i,4] > 100) Wfile[i,4] <- 100
if (Wfile[i,4] > RHthreshold) Wfile[i,7] <- (1 / NoReadingsH)
}
if (!is.na(Wfile[i,7])) {
if (Wfile[i,7] == (1 / NoReadingsH)) Wfile[i,8] <- Wfile[i,3]
}
}
Rainfall <-round(as.matrix(by( Wfile[,5],Wfile[,1],function(x) sum(x,na.rm=TRUE))),1)
tmp <- as.character(row.names(Rainfall))
Date<-as.Date(tmp)
Tmp<- as.matrix(by( Wfile[,3],Wfile[,1],function(x) mean(x,na.rm=TRUE)))
HumidHrs <-as.matrix(by( Wfile[,7],Wfile[,1],function(x) sum(x,na.rm=TRUE)))
humidtmp <- as.matrix(by( Wfile[,8],Wfile[,1],function(x) mean(x,na.rm=TRUE)))
Wfile <- data.frame(Date,Rainfall,Tmp,HumidHrs,humidtmp)
return(list(Wfile=Wfile,Sfile=Sfile,EmergDate=EmergDate,EndEpidDate=EndEpidDate))
}
