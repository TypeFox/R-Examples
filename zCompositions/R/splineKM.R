splineKM <- function(x,label=NULL,dl=NULL,n.knots=NULL,
                     legend.pos="bottomright",
                     ylab="ECDF",
                     xlab="Value",
                     col.km="black",lty.km=1,lwd.km=1,
                     col.sm="red",lty.sm=2,lwd.sm=2,...){
  
  if (is.character(dl)) stop("The argument dl must be a numeric vector")
  if (length(dl)!=length(x)) stop("x and dl must be two vectors of the same length")
  
  if (is.null(label)) stop("A value for label must be given")
  if (!is.na(label)){
    if (label!=0 & any(x==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data")
    if (any(is.na(x))) stop(paste("NA values not labelled as censored values were found in the data"))
  }
  if (is.na(label)){
    if (any(x==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data")
    if (!any(is.na(x),na.rm=T)) stop(paste("Label",label,"was not found in the data"))
  }
  
  if ((!is.null(n.knots)) & (length(n.knots)!=1)) stop("n.knots must contain a single value")
  
  x[x==label] <- NA
  
  who <- is.na(x); w <- which(who)
  
  xcen <- ifelse(who,TRUE,FALSE)
  x[who] <- dl[who]
  
  dat <- data.frame(x,xcen)
  km.ecdf <- cenfit(dat$x,dat$xcen)
  
  x <- rev(km.ecdf@survfit$time) 
  y <- rev(km.ecdf@survfit$surv)
  
  if (is.null(n.knots)) {scdf <- smooth.spline(x,y)}
  if (!is.null(n.knots)) {scdf <- smooth.spline(x,y,nknots=n.knots)}
  scdf <- approxfun(scdf$x,scdf$y)
  
  plot(km.ecdf,conf.int=FALSE,ylab=ylab,xlab=xlab,
       col=col.km,lty=lty.km,lwd=lwd.km, ...)
  lines(x,scdf(x),type="l",
       col=col.sm,lty=lty.sm,lwd=lwd.sm)
  abline(h=1,col="white",lwd=4)
  legend(legend.pos,bty="n",
         legend=c("KM estimate","KMSS estimate"),
         lty=c(lty.km,lty.sm),col=c(col.km,col.sm),lwd=c(lwd.km,lwd.sm)) 

}