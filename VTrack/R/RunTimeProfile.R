RunTimeProfile <-
function(sInputFile,sDATETIMEField,sProfilePeriod)
{
  if(sProfilePeriod == "hour"){DT <- as.character(strptime(sInputFile[,c(sDATETIMEField)],format="%Y-%m-%d %H"))}
  if(sProfilePeriod == "day"){DT <- as.character(strptime(sInputFile[,c(sDATETIMEField)],format="%Y-%m-%d"))}
  if(sProfilePeriod == "circadian"){DT <- substring(sInputFile[,sDATETIMEField],12,13)}
  if(sProfilePeriod == "month"){DT <- substring(sInputFile[,sDATETIMEField],6,7)}
  if(sProfilePeriod == "week"){DT <- ceiling(strptime(sInputFile[,c(sDATETIMEField)],format="%Y-%m-%d")$yday/7)}

  A1<-aggregate(sInputFile[,c("DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),length)
  DATETIME<-A1[1]
  TRANSMITTERID<-A1[2]
  FREQ<-A1[3]
  newtable<-data.frame(DATETIME,TRANSMITTERID,FREQ)
  
  if(names(sInputFile)[3]!="SENSOREVENT")
  {
    A1<-aggregate(sInputFile[,c("DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),sum)
    A2<-aggregate(sInputFile[,c("DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),max)
    A3<-aggregate(sInputFile[,c("DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),mean)
    A4<-aggregate(sInputFile[,c("DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),sd)
    TIMESUM<-A1[3]
    TIMEMAX<-A2[3]
    TIMEAV<-A3[3]
    TIMESTDEV<-A4[3]
    newtable<-cbind(newtable,TIMESUM,TIMEMAX,TIMEAV,TIMESTDEV)
    names(newtable)<-c("DATETIME","TRANSMITTERID","FREQ","TIMESUM","TIMEMAX","TIMEAV","TIMESTDEV")
  }
  if(names(sInputFile)[3]=="SENSOREVENT")
  {
    iminmax1 <- ifelse(names(sInputFile)[9]=="MAXSENSOR","MAXSENSOR","MINSENSOR")
    iminmax2 <- ifelse(names(sInputFile)[9]=="MAXSENSOR","SENSORMAX","SENSORMIN")
    A1<-aggregate(sInputFile[,c(iminmax1,"DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),sum)
    A2<-aggregate(sInputFile[,c(iminmax1,"DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),max)
    A2b<-aggregate(sInputFile[,c(iminmax1,"DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),min)
    A3<-aggregate(sInputFile[,c(iminmax1,"DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),mean)
    A4<-aggregate(sInputFile[,c(iminmax1,"DURATION")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),sd)
    TIMESUM<-A1[4]
    if(iminmax1=="MAXSENSOR") SENSORMINMAX<-A2[3]
    if(iminmax1=="MINSENSOR") SENSORMINMAX<-A2b[3]
    TIMEMAX<-A2[4]
    SENSORAV<-A3[3]
    TIMEAV<-A3[4]
    SENSORSTDEV<-A4[3]
    TIMESTDEV<-A4[4]
    newtable<-cbind(newtable,SENSORMINMAX,SENSORAV,SENSORSTDEV,TIMESUM,TIMEMAX,TIMEAV,TIMESTDEV)
    names(newtable)<-c("DATETIME","TRANSMITTERID","FREQ",iminmax2,"SENSORAV","SENSORSTDEV","TIMESUM","TIMEMAX","TIMEAV","TIMESTDEV")
  }
  if(length(which(names(sInputFile)=="NUMRECS"))>0)
  {
    A5 <- aggregate(sInputFile[,c("NUMRECS")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),sum)
    DETECTIONS <- A5[3]
    newtable <- cbind(newtable,DETECTIONS)
    names(newtable)[ncol(newtable)]<-"DETECTIONS"
    return(newtable)
  }
  if(length(which(names(sInputFile)=="DISTANCE"))>0)
  {
    A5 <- aggregate(sInputFile[,c("DISTANCE")],by=list(Date=DT,TRANSMITTERID=sInputFile$TRANSMITTERID),sum)
    DISTANCE <- A5[3]
    newtable <- cbind(newtable,DISTANCE)
    names(newtable)[ncol(newtable)]<-"DISTANCE"
    return(newtable)
  }else{
    return(newtable)
  }
}
