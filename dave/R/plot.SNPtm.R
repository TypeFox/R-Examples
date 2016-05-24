plot.SNPtm <-
function(x,...) {
   o.SNPtm<- x
   par(omi=c(3,0,0,0))
   line.wd<- seq(0.5,3.0,0.5) 
   symbols<- c(1,4,8,17,22,18)
   colors<- c("darkolivegreen4","lightgreen","gold","darkorange","red1","darkred")
#
   x<- o.SNPtm$sim.data
   tv<- o.SNPtm$time.vector
   vegtypes<- o.SNPtm$veg.types

# Plot in color
   par(omi=c(3,0,0,0))
   plot(c(min(tv),max(tv)),c(0,100),type="n",xlab="Time (yr)",ylab="Cover %")
       for(v in 1:6){
           points(tv,x[,v],pch=symbols[7-v],cex=0.6,col=colors[7-v])
       }
  legend("topleft",vegtypes,pch=symbols,col=colors,cex=0.8,ncol=2,bty="n")


# Plot, bw
   x<- o.SNPtm$sim.data
   par(omi=c(3,0,0,0))
   colors<- rgb(0,0,0,seq(100,250,20),maxColorValue=255)
   line.wd<- seq(0.5,3.0,0.5) 
   colors<- rgb(0,0,0,seq(100,250,20),maxColorValue=255)
   plot(c(min(tv),max(tv)),c(0,100),type="n",xlab="Time (yr)",ylab="Cover %")
       for(v in 1:6){
           lines(tv,x[,v],lwd=line.wd[7-v],col=colors[7-v])
       }
   legend("topleft",vegtypes,lwd=line.wd,col=colors,cex=0.8,ncol=2,bty="n")
  }
