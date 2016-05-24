#rm(list=ls())
#source("/home/el/lavori/Rdevel/overlapping1.1/R/cutnumeric.R")
#load("~/lavori/Rdevel/overlappingTest/provaoverlap.rda")
final.plot <- function(DD,OV) {
  
  xx <- yy <- list()
  JJ <- unique(DD$k)
  for (i in JJ) {
    xxx <- DD$x[(DD$k==i)&(DD$w==1)]
    yyy <- DD$y[(DD$k==i)&(DD$w==1)]
    xx <- c(xx,list(sort(xxx)))
    yy <- c(yy,list(yyy[order(xxx)]))
  }
  
  xyplot(y~x|factor(k),data=DD,groups=DD$j,panel=function(...){
    panel.xyplot(...,type="l")
    panel.polygon(xx[[packet.number()]],yy[[packet.number()]],col="#74d600",border="#74d600")
    panel.text(min(DD$x),max(DD$y),
               paste("overlap = ",round(OV[packet.number()]*100,2),"%",sep=""),pos=4)
  },xlab="",ylab="")
}
