plottest.idturtle <-
function (test) 
{
  measurements<-as.factor(test[,c("NMeasurements")])
  measurements<-data.frame(summary(measurements,maxsum=length(measurements)))
  measurements[2]<-rownames(measurements)
  names(measurements)<-list("measurements","N")
  turtles<-as.factor(test[,c("NTurtles")])
  turtles<-data.frame(summary(turtles,maxsum=length(turtles)))
  turtles[2]<-rownames(turtles)
  names(turtles)<-list("Turtles","N")
  maximo<-max(as.numeric(measurements[,c("N")]))
  figure<-data.frame()
  for (i in 1:5){
    figure[i]<-vector()
  }
  names(figure)<-list("N","measurements","Turtles","measurements%","Turtles%")
  for (i in 0:maximo){
    figure[i+1,1]<-i
    if(length(measurements[measurements$N==i,c("measurements")])==0){
      figure[i+1,2]<-0
    }
    else{
      figure[i+1,2]<-as.numeric(measurements[measurements$N==i,c("measurements")])
    }  
    if(length(turtles[turtles$N==i,c("Turtles")])==0){
      figure[i+1,3]<-0
    }
    else{
      figure[i+1,3]<-as.numeric(turtles[turtles$N==i,c("Turtles")])
    }
  }
  figure[4]<-100*figure[2]/nrow(test)
  figure[5]<-100*figure[3]/nrow(test)
  print("Plotting")
  barplot(figure[,4],main="Number of measurement sets",xlab="Measurements",ylab="Percentage",col='gray30',ylim = c(0,100),names.arg=figure[,1])
  barplot(figure[,5],main="Number of turtles",xlab="Turtles",ylab="Percentage",col='gray30',ylim = c(0,100),names.arg=figure[,1])
  
  print("Writing file in hard disc")
  write.table(figure, file="./plottest.csv", sep=",", row.names=FALSE, quote=FALSE)
  
  print("Done")
  return(figure)
}
