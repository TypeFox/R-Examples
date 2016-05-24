features.plot <-
function(R,first=FALSE,para=c(5,40,50,10,10,10,10,5,1,5,5,5,10),skip=FALSE,column=4,row=3){
  
  if(skip==TRUE) par(mfrow=c(row+1,column))
  if(skip==FALSE) par(mfrow=c(row,column))
  
  Syllabic<-vector()
  Place<-vector()
  Stop<-vector()
  Voice<-vector()
  Nasal<-vector()
  Retroflex<-vector()
  Lateral<-vector()
  Aspirated<-vector()
  Long<-vector()
  High<-vector()
  Back<-vector()
  Round<-vector()
  SkipCost<-vector()
  
  if(first==TRUE){
    for(i in 1:length(R)){
      Syllabic[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,1]),decreasing=TRUE))[1])
      Place[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,2]),decreasing=TRUE))[1])
      Stop[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,3]),decreasing=TRUE))[1])
      Voice[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,4]),decreasing=TRUE))[1])
      Nasal[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,5]),decreasing=TRUE))[1])
      Retroflex[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,6]),decreasing=TRUE))[1])
      Lateral[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,7]),decreasing=TRUE))[1])
      Aspirated[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,8]),decreasing=TRUE))[1])
      Long[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,9]),decreasing=TRUE))[1])
      High[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,10]),decreasing=TRUE))[1])
      Back[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,11]),decreasing=TRUE))[1])
      Round[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,12]),decreasing=TRUE))[1])
      SkipCost[i]<-as.integer(names(sort(table(R[[i]]$optimized_parameters[,13]),decreasing=TRUE))[1]) 
    }
  }
  
  if(first==FALSE){
    for(i in 1:length(R)){
      Syllabic<-c(Syllabic,R[[i]]$optimized_parameters[,1])
      Place<-c(Place,R[[i]]$optimized_parameters[,2])
      Stop<-c(Stop,R[[i]]$optimized_parameters[,3])
      Voice<-c(Voice,R[[i]]$optimized_parameters[,4])
      Nasal<-c(Nasal,R[[i]]$optimized_parameters[,5])
      Retroflex<-c(Retroflex,R[[i]]$optimized_parameters[,6])
      Lateral<-c(Lateral,R[[i]]$optimized_parameters[,7])
      Aspirated<-c(Aspirated,R[[i]]$optimized_parameters[,8])
      Long<-c(Long,R[[i]]$optimized_parameters[,9])
      High<-c(High,R[[i]]$optimized_parameters[,10])
      Back<-c(Back,R[[i]]$optimized_parameters[,11])
      Round<-c(Round,R[[i]]$optimized_parameters[,12])
      SkipCost<-c(SkipCost,R[[i]]$optimized_parameters[,13])
    }
  }
  
  m<-selection(R)
  
  hist(Syllabic,xlim=c(0,100))
  abline(v=c(para[1],m[1]),lwd=1.5,col=c("red","purple"))
  hist(Place,xlim=c(0,100))
  abline(v=c(para[2],m[2]),lwd=1.5,col=c("red","purple"))
  hist(Stop,xlim=c(0,100))
  abline(v=c(para[3],m[3]),lwd=1.5,col=c("red","purple"))
  hist(Voice,xlim=c(0,100))
  abline(v=c(para[4],m[4]),lwd=1.5,col=c("red","purple"))
  hist(Nasal,xlim=c(0,100))
  abline(v=c(para[5],m[5]),lwd=1.5,col=c("red","purple"))
  hist(Retroflex,xlim=c(0,100))
  abline(v=c(para[6],m[6]),lwd=1.5,col=c("red","purple"))
  hist(Lateral,xlim=c(0,100))
  abline(v=c(para[7],m[7]),lwd=1.5,col=c("red","purple"))
  hist(Aspirated,xlim=c(0,100))
  abline(v=c(para[8],m[8]),lwd=1.5,col=c("red","purple"))
  hist(Long,xlim=c(0,100))
  abline(v=c(para[9],m[9]),lwd=1.5,col=c("red","purple"))
  hist(High,xlim=c(0,100))
  abline(v=c(para[10],m[10]),lwd=1.5,col=c("red","purple"))
  hist(Back,xlim=c(0,100))
  abline(v=c(para[11],m[11]),lwd=1.5,col=c("red","purple"))
  hist(Round,xlim=c(0,100))
  abline(v=c(para[12],m[12]),lwd=1.5,col=c("red","purple"))
  if(skip==TRUE){
    hist(SkipCost,xlim=c(0,100))
    abline(v=c(para[13],m[13]),lwd=1.5,col=c("red","purple"))
  }
}
