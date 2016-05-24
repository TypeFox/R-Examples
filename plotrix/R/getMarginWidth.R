getMarginWidth<-function(side=4,labels,is.legend=FALSE) {
 currentmar<-par("mar")
 currentpin<-par("pin")[1]
 currentfin<-par("fin")[1]
 currentusr<-par("usr")
 if(is.legend) marwidth<-1.2*legend(0,0,labels,fill=NA,plot=FALSE)$rect$w
 else marwidth<-1.1*max(strwidth(labels))
 marprop<-ifelse(side==2,currentmar[2]/(currentmar[2]+currentmar[4]),
  currentmar[4]/(currentmar[2]+currentmar[4]))*(currentfin-currentpin)/currentfin
 plotprop<-currentpin/currentfin
 plotwidth<-currentusr[2]-currentusr[1]
 marusr<-(marprop/plotprop)*plotwidth
 newmar<-ifelse(side==2,currentmar[2],currentmar[4])*marwidth/marusr
 cat("plotprop",plotprop,"marprop",marprop,"plotwidth",
  plotwidth,"marwidth",marwidth,"\n")
 marcenter<-currentusr[2]+marwidth/2
 if(marwidth>plotwidth) {
  warning("figure size too small for new margin!")
  newmar<-ifelse(side==2,currentmar[2],currentmar[4])
 }
 return(list(newmar=newmar,marcenter=marcenter))
}
