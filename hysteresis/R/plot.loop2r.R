plot.loop2r <- function (x,split.line=TRUE,xlim=NULL,ylim=NULL,putNumber=FALSE,main=NULL,values=NULL,...) {
  a <- x
  ti <- (1:101)*pi/50
  if (a$method=="harmonic2") {
  Ind <- (ti < pi) & (ti > 0)
  Input <- a$values["b.x"]*cos(ti)+a$values["cx"]
  if (a$extended.classical==FALSE) Output <- a$values["b.y"]*cos(ti)^a$values["n"]+Ind*a$values["retention.above"]*sin(ti)^a$values["m"]+(1-Ind)*a$values["retention.below"]*sin(ti)^a$values["m"]+a$values["cy"]
  else Output <- sign(cos(ti))*a$values["b.y"]*abs(cos(ti))^a$values["n"]+Ind*a$values["retention.above"]*sin(ti)^a$values["m"]+(1-Ind)*a$values["retention.below"]*sin(ti)^a$values["m"]+a$values["cy"]
  }
  else {
    costp <- cos(ti) 
    sintp <- sin(ti) 
    direc <- sign(costp)
    direcsin <- sign(sintp)
    Input <- a$values["cx"]+a$values["b.x"]*costp 
    Output <- a$values["cy"]+(direcsin < 0)*direcsin*a$values["retention.below"]*abs(sintp)^a$values["m"]+(direcsin > 0)*direcsin*a$values["retention.above"]*abs(sintp)^a$values["m"]+direc*a$values["b.y"]*abs(costp)^a$values["n"]
    
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
      mtext(paste("Retention.above=",format(a$values["retention.above"],digits=3)," Retention.below=",format(a$values["retention.below"],digits=3)),side=3,line=0.0,cex=0.75)
    }
  else  
    if (values=="hysteresis" | values=="hysteresis.all") {
      plot(Output~Input,type="l",ylim=ylim,xlim=xlim,...)
      title(line=3, paste(main),cex=1.2)
      mtext(paste(
        "b.x=",format(a$values["b.x"],digits=3)," b.y=",format(a$values["b.y"],digits=3)," cx=",format(a$values["cx"],digits=3)," cy=",format(a$values["cy"],digits=3)),side=3,line=1.85,cex=0.75)
      mtext(paste("Area=",format(a$values["area"],digits=3)," Lag above=",format(a$values["lag.above"],digits=3)," Lag below=",format(a$values["lag.below"],digits=3)),side=3,line=0.95,cex=0.75)
      mtext(paste("Retention.above=",format(a$values["retention.above"],digits=3)," Retention.below=",format(a$values["retention.below"],digits=3)," Coercion above=",format(a$values["coercion.above"],digits=3)," Coercion below=",format(a$values["coercion.below"],digits=3)),side=3,line=0.0,cex=0.75)
    }
  else
  if (values=="derived") {
    plot(Output~Input,type="l",ylim=ylim,xlim=xlim,...)
    title(line=3, paste(main),cex=1.2)
    mtext(paste(
      "Coercion above=",format(a$values["coercion.above"],digits=3)," Coercion below=",format(a$values["coercion.below"],digits=3)," Area=",format(a$values["area"],digits=3)),side=3,line=1.85,cex=0.75)
    mtext(paste(" Lag above=",format(a$values["lag.above"],digits=3)," Lag below=",format(a$values["lag.below"],digits=3)," Split Angle=",format(a$values["beta.split.angle"],digits=3)),side=3,line=0.95,cex=0.75)
    mtext(paste("Hysteresis x above=",format(a$values["hysteresis.x.above"],digits=3),"Hysteresis x below=",format(a$values["hysteresis.x.below"],digits=3)," Hysteresis y above=",format(a$values["hysteresis.y.above"],digits=3),"Hysteresis y below=",format(a$values["hysteresis.y.below"],digits=3)),side=3,line=0.0,cex=0.75)
  }
  else plot(Output~Input,type="l",ylim=ylim,xlim=xlim,main=main,...)
  }
 
  points(a$y~a$x,pch=1,cex=0.85)
  if (split.line==TRUE) {

  lines(Input,splitLine,lty=2)}
  if(putNumber==TRUE) text(a$x,a$y,as.character(format(1:length(a$y),digits=4)))
}
