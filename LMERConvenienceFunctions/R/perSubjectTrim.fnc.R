##############################
#### PER SUBJECT TRIMMING ####
##############################
perSubjectTrim.fnc<-function(
      data,   #unquoted
      response, #quoted
      subject,  #quoted
      trim=2.5
){
  data0<-data
  sdev=as.data.frame(tapply(data[,response],data[,subject],sd))
  cat(".")
  sdev$Subject=row.names(sdev)
  colnames(sdev)[1]="SD"
  row.names(sdev)=1:length(sdev$SD)
  sdev$SD=as.numeric(sdev$SD)
  sdev$Subject=as.factor(sdev$Subject)
  row.names(sdev)=1:length(sdev$SD)
  sdev=na.omit(sdev)
  sdev=sdev[,c(2,1)]
  data=merge(data,sdev,by.x=subject,by.y="Subject")
  cat(".")
  m=as.data.frame(tapply(data[,response],data[,subject],mean))
  m$Subject=row.names(m)
  colnames(m)[1]="Mean"
  row.names(m)=1:length(m$Mean)
  m=m[,c(2,1)]
  m=na.omit(m)
  data=merge(data,m,by.x=subject,by.y="Subject")
  cat(".")
  data$SD=as.numeric(data$SD)
  data$Mean=as.numeric(data$Mean)
  Scaled=apply(data[,c(response,"SD","Mean")],1,function(row)(row[1]-row[3])/row[2])
  cat(".")
  data=cbind(data,Scaled)
  cat(".")
  outliers = as.numeric(rownames(data[abs(data$Scaled) > trim, ]))
  data=data[-outliers,,drop=TRUE]
  cat(".\n")
  cat("n.removed =",(nrow(data0)-nrow(data)),"\n")
  cat("percent.removed =",(nrow(data0)-nrow(data))/nrow(data0)*100,"\n")

  return(list(data=data,data0=data0,n.removed=nrow(data0)-nrow(data),percent.removed=(nrow(data0)-nrow(data))/nrow(data0)*100))

}
