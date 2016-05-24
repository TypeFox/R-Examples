plot.ellipsefit<-function(x,putNumber=FALSE,values=NULL,xlim=NULL,ylim=NULL,main=NULL,newPred=TRUE,show=NULL,split.line=FALSE,...)
  {
  a <- x
  if (newPred==TRUE)
  {
ti <- (1:101)*pi/50
newX <- a$values["b.x"]*cos(ti)+a$values["cx"]
newY <- a$values["b.y"]*cos(ti)+a$values["retention"]*sin(ti)+a$values["cy"]
}
else {
newY <- a$pred.y
newY[length(newY)+1] <- newY[1]
newX <- a$pred.x
newX[length(newX)+1] <- newX[1] }
  if (is.null(xlim)) xlim <-c(min(c(a$x,newX)),max(c(a$x,newX)))
  if (is.null(ylim)) ylim <- c(min(c(a$y,newY)),max(c(a$y,newY)))                           
if (!is.null(values)) {
    if (values=="inherent") {
      plot(newY~newX,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("cx=",format(a$values["cx"],digits=3)," cy=",format(a$values["cy"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Retention=",format(a$values["retention"],digits=3)),side=3,line=0.0,cex=0.75)
    }
    
   else if (values=="hysteresis") {
      plot(newY~newX,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)," cx=",format(a$values["cx"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("cy=",format(a$values["cy"],digits=3)," Area=",format(a$values["area"],digits=3)," Lag=",format(a$values["lag"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Retention=",format(a$values["retention"],digits=3)," Coercion=",format(a$values["coercion"],digits=3)),side=3,line=0.0,cex=0.75)
    }
    
   else if (values=="hysteresis.all") {
      plot(newY~newX,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)," cx=",format(a$values["cx"],digits=3)," cy=",format(a$values["cy"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("Area=",format(a$values["area"],digits=3)," Lag=",format(a$values["lag"],digits=3)," Retention=",format(a$values["retention"],digits=3)," Split Angle=",format(a$values["beta.split.angle"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Coercion=",format(a$values["coercion"],digits=3)," Hysteresis x=",format(a$values["hysteresis.x"],digits=3)," Hysteresis y=",format(a$values["hysteresis.y"],digits=3)),side=3,line=0.0,cex=0.75)
    }
  
 else if (values=="derived") {
    plot(newY~newX,type="l",ylim=ylim,xlim=xlim,...)
    title(line=3, paste(main),cex=1.2)
    mtext(paste(
      "Coercion=",format(a$values["coercion"],digits=3)," Area=",format(a$values["area"],digits=3)),side=3,line=1.85,cex=0.75)
    mtext(paste("Lag=",format(a$values["lag"],digits=3)," Split Angle=",format(a$values["beta.split.angle"],digits=3)),side=3,line=0.95,cex=0.75)
    mtext(paste("Hysteresis x=",format(a$values["hysteresis.x"],digits=3)," Hysteresis y=",format(a$values["hysteresis.y"],digits=3)),side=3,line=0.0,cex=0.75)
  }
  
 else if (values=="ellipse") {
    plot(newY~newX,type="l",ylim=ylim,xlim=xlim,...)
    title(line=3, paste(main),cex=1.2)
    mtext(paste(
      "Ampx=",format(a$values["ampx"],digits=3)," Ampy=",format(a$values["ampy"],digits=3)),side=3,line=1.85,cex=0.75)
    mtext(paste("rote.deg=",format(a$values["rote.deg"],digits=3)," Eccentricity=",format(a$values["eccentricity"],digits=3)),side=3,line=0.95,cex=0.75)
    mtext(paste("S-major Axis=",format(a$values["semi.major"],digits=3)," S-minor Axis=",format(a$values["semi.minor"],digits=3)),side=3,line=0.0,cex=0.75)
  }

  else  if (values=="ellipse.all") {
      plot(newY~newX,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste("Cx=",format(a$values["cx"],digits=3)," Cy=",format(a$values["cy"],digits=3),
        " Ampx=",format(a$values["ampx"],digits=3)," Ampy=",format(a$values["ampy"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("rote.deg=",format(a$values["rote.deg"],digits=3)," focus.x=",format(a$values["focus.x"],digits=3)," focus.y=",format(a$values["focus.y"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("S-major Axis=",format(a$values["semi.major"],digits=3)," S-minor Axis=",format(a$values["semi.minor"],digits=3)," Eccentricity=",format(a$values["eccentricity"],digits=3)),side=3,line=0.0,cex=0.75)
    }
 else plot(newY~newX,type="l",ylim=ylim,xlim=xlim,main=main,...)
 }
 else plot(newY~newX,type="l",ylim=ylim,xlim=xlim,main=main,...) 
 
points(a$y~a$x,pch=1,cex=0.85)

if (any(show=="semi.major")) segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["semi.major"]*cos(a$values["rote.deg"]/180*pi),a$values["cy"]+a$values["semi.major"]*sin(a$values["rote.deg"]/180*pi),col="red")
if (any(show=="semi.minor")) segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["semi.minor"]*cos(a$values["rote.deg"]/180*pi+pi/2),a$values["cy"]+a$values["semi.minor"]*sin(a$values["rote.deg"]/180*pi+pi/2),col="red")

if (any(show %in% c("b.x","b.y"))) segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["b.x"],a$values["cy"]+a$values["b.y"],col="blue")
  if (any(show %in% c("focus.x","focus.y"))) points(c(a$values["cx"]+a$values["focus.x"],a$values["cx"]-a$values["focus.x"]),c(a$values["cy"]+a$values["focus.y"],a$values["cy"]-a$values["focus.y"]),col="gold",cex=2,pch=19)
  if (any(show=="rote.deg")) {arrows(a$values["cx"]+a$values["coercion"],a$values["cy"],a$values["cx"]+a$values["focus.x"],a$values["cy"]+a$values["focus.y"],lty=2)
                              segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["coercion"],a$values["cy"],lty=2)
}
  if (any(show=="retention")) segments(a$values["cx"],a$values["cy"],a$values["cx"],a$values["cy"]+a$values["retention"],col="purple")
  
  if (any(show=="coercion")) segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["coercion"],a$values["cy"],col="green")
  
if(putNumber==TRUE) text(a$x,a$y,as.character(format(1:length(a$y),digits=4)))

if (split.line==TRUE) {
tiS <- (1:101)*pi/50
newXS <- a$values["b.x"]*cos(tiS)+a$values["cx"]
splitLine <- a$values["b.y"]*cos(tiS)+a$values["cy"]
lines(newXS,splitLine,lty=2)
}

}
