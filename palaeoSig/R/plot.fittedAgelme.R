#'Plots fitted agelme model and dates 
#'@param x Fitted agelme model.
#'@param main Title of the plot.
#'@param xlab xlabel of the plot.
#'@param ylab ylabel of the plot.
#'@param \dots Other arguments to plot.
#'@examples 
#'data(STOR)

#'fit.mod <- with(STOR,agelme(depthup,depthdo,cageup,cagedo))

#'#Predicting using the constant variance model,
#'#for each cm between 70 and 400 cm.
#'fit.pre <- predict(fit.mod,1,70:400)
#'plot(fit.pre)
#'

plot.fittedAgelme<-function(x,main,xlab="Depth", ylab="Calibrated Age",...){
  if(missing(main))main=x$v  
  xlim=range(c(x$data$depup, x$data$depthdo, x$fit$depth))
  ylim=range(c(x$data$cageup, x$data$cagedo,x$fit[,2:4])) 
  with(x$data,plot((depthup+depthdo)/2,(cageup+cagedo)/2,xlim=xlim, ylim=ylim,ylab=ylab,xlab=xlab,pch=20,...))
  with(x$data,arrows(x0=depthup, y0=(cageup+cagedo)/2, x1=depthdo, y1=(cageup+cagedo)/2, length=0, col=ifelse(use, "black", "grey70")))
  with(x$data,arrows(x0=(depthup+depthdo)/2, y0=cageup, x1=(depthup+depthdo)/2, y1=cagedo, length=0, col=ifelse(use,"black","grey70")))
  title(main=main)
  matlines(x$fit$depth,x$fit[,2:4],lty=c(1,2,2),col=c(2,4,4))
}
