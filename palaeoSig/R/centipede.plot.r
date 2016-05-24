#library(rioja)
#data(SWAP)
#mod<-WA(SWAP$spec, SWAP$pH, tolDW=TRUE)
#coef(mod)
#centipede.plot(mod, keep=colSums(SWAP$spec>0)>40, pch=20, cex.axis=.7)


centipede.plot<-function(x,keep = TRUE, xlab="", xlim, ...){
   co<-coef(x)[keep,]
   co<-co[order(co[,"Optima"]),]
   if(missing(xlim))xlim<-range(c(co[,"Optima"]-co[,"Tolerances"], co[,"Optima"]+co[,"Tolerances"]))
   plot(co[,"Optima"], 1:nrow(co), xlab=xlab, yaxt="n",ylab="",xlim=xlim,...)
   arrows(x0=co[,"Optima"]-co[,"Tolerances"], y0=1:nrow(co), x1=co[,"Optima"]+co[,"Tolerances"], y1=1:nrow(co), length=0)
   axis(2, at=1:nrow(co), labels=rownames(co),las=2, ...)
}
