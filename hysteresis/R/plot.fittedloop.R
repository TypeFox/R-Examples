plot.fittedloop <- function (x,split.line=TRUE,xlim=NULL,ylim=NULL,putNumber=FALSE,values=NULL,main=NULL,show=NULL,...) {
  a <- x
  ti <- (1:101)*pi/50
  
  if (a$method=="harmonic2") {
  Input <- a$values["b.x"]*cos(ti)+a$values["cx"]
  if (a$extended.classical==FALSE) Output <- a$values["b.y"]*cos(ti)^a$values["n"]+a$values["retention"]*sin(ti)^a$values["m"]+a$values["cy"]
  else Output <- sign(cos(ti))*a$values["b.y"]*abs(cos(ti))^a$values["n"]+a$values["retention"]*sin(ti)^a$values["m"]+a$values["cy"]
  } 
  else {
    costp <- cos(ti) 
    sintp <- sin(ti) 
    direc <- sign(costp)
    direcsin <- sign(sintp)
    Input <- a$values["cx"]+a$values["b.x"]*costp 
    Output <- a$values["cy"]+direcsin*a$values["retention"]*abs(sintp)^a$values["m"]+direc*a$values["b.y"]*abs(costp)^a$values["n"]
    
  }
  if (a$extended.classical==FALSE & a$method=="harmonic2") splitLine <- a$values["b.y"]*cos(ti)^a$values["n"]+a$values["cy"]
     else splitLine <- sign(cos(ti))*a$values["b.y"]*abs(cos(ti))^a$values["n"]+a$values["cy"]
  
  if (is.null(xlim)) xlim <-range(c(a$x,Input))
  if (is.null(ylim)) ylim <- range(c(a$y,Output,splitLine))                           
 if (is.null(values)) plot(Output~Input,type="l",ylim=ylim,xlim=xlim,main=main,...)
  else {
    if (values=="inherent") {
      plot(Output~Input,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("cx=",format(a$values["cx"],digits=3)," cy=",format(a$values["cy"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Retention=",format(a$values["retention"],digits=3)),side=3,line=0.0,cex=0.75)
    }
    
   else if (values=="hysteresis") {
      plot(Output~Input,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)," cx=",format(a$values["cx"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("cy=",format(a$values["cy"],digits=3)," Area=",format(a$values["area"],digits=3)," Lag=",format(a$values["lag"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Retention=",format(a$values["retention"],digits=3)," Coercion=",format(a$values["coercion"],digits=3)),side=3,line=0.0,cex=0.75)
    }
    
  else  if (values=="hysteresis.all") {
      plot(Output~Input,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)," cx=",format(a$values["cx"],digits=3)," cy=",format(a$values["cy"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("Area=",format(a$values["area"],digits=3)," Lag=",format(a$values["lag"],digits=3)," Retention=",format(a$values["retention"],digits=3)," Split Angle=",format(a$values["beta.split.angle"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Coercion=",format(a$values["coercion"],digits=3)," Hysteresis x=",format(a$values["hysteresis.x"],digits=3)," Hysteresis y=",format(a$values["hysteresis.y"],digits=3)),side=3,line=0.0,cex=0.75)
    }
  
  else if (values=="derived") {
    plot(Output~Input,type="l",ylim=ylim,xlim=xlim,...)
    title(line=3, paste(main),cex=1.2)
    mtext(paste(
      "Coercion=",format(a$values["coercion"],digits=3)," Area=",format(a$values["area"],digits=3)),side=3,line=1.85,cex=0.75)
    mtext(paste("Lag=",format(a$values["lag"],digits=3)," Split Angle=",format(a$values["beta.split.angle"],digits=3)),side=3,line=0.95,cex=0.75)
    mtext(paste("Hysteresis x=",format(a$values["hysteresis.x"],digits=3)," Hysteresis y=",format(a$values["hysteresis.y"],digits=3)),side=3,line=0.0,cex=0.75)
  }
  else plot(Output~Input,type="l",ylim=ylim,xlim=xlim,main=main,...)
  }
  
  if (any(show %in% c("b.x","b.y"))) segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["b.x"],a$values["cy"]+a$values["b.y"],col="blue")
  if (any(show=="retention")) segments(a$values["cx"],a$values["cy"],a$values["cx"],a$values["cy"]+a$values["retention"],col="purple") 
  if (any(show=="coercion")) segments(a$values["cx"],a$values["cy"],a$values["cx"]+a$values["coercion"],a$values["cy"],col="green")
  
  points(a$y~a$x,pch=1,cex=0.85)
  if (split.line==TRUE) {
  lines(Input,splitLine,lty=2)}
  if(putNumber==TRUE) text(a$x,a$y,as.character(format(1:length(a$y),digits=4)))
}
