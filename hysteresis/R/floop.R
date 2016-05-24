floop <- function(x,y=NULL,n=1,m=1,times="equal",period=NULL,subjects=NULL, subset=NULL,na.action=getOption("na.action"),extended.classical=FALSE,boot=FALSE,method="harmonic2",...) {
 if (boot==TRUE) return(summary(floop(x,y,n,m,times,period,subjects,subset,na.action,extended.classical,method=method),...))
  if (m==1 & n==1) return(fel(x,y,times=times,period=period,subjects=subjects,subset=subset,na.action=na.action,method=method))
  floopcall <- match.call()
  if (method!="geometric" & method!="harmonic2") stop("method should be either geometric or harmonic2")
  if (ncol(matrix(x)) > 2)
    times <- x[,3]
  dat <- xy.coords(x,y)
  if (!is.null(subset)) {
    dat$x<-dat$x[subset]; dat$y<-dat$y[subset];
    if (!is.null(subjects)) {
      if (!is.list(subjects)) {
        subjects<-factor(subjects[subset])
      }
    else subjects <- lapply(subjects,function (x) factor(x[subset])) }   
    if (is.numeric(times))
      times<-times[subset]}
  if (!is.null(subjects)) {
    dat <- cbind("x"=dat$x,"y"=dat$y)
    if (is.numeric(times))
      ans <- by(cbind(dat,times),subjects,floop,m=m,n=n,period=period,na.action=na.action,extended.classical=extended.classical)
    else
      ans <- by(dat,subjects,floop,m=m,n=n,period=period,times=times,na.action=na.action,extended.classical=extended.classical)
    if (!is.list(subjects)) names(ans) <- levels(subjects)
    
    values <- t(sapply(ans,function (x) x["values"]$values))
    Std.Errors <- t(sapply(ans,function (x) x["Std.Errors"]$Std.Errors))
    
    if (is.list(subjects)){
      subjectmat <- matrix(NA,nrow=prod(dim(ans)),ncol=length(subjects))
      if (length(subjects) > 2)  {
        for (i in 2:(length(subjects)-1)){
          subjectmat[,i] <- rep(dimnames(ans)[[i]],times=prod(dim(ans)[(i+1):length(subjects)]),each=prod(dim(ans)[1:(i-1)]))
        }}
      subjectmat[,1] <- rep(dimnames(ans)[[1]],times=prod(dim(ans)[2:length(subjects)]))
      subjectmat[,length(subjects)] <- rep(dimnames(ans)[[length(subjects)]],each=prod(dim(ans)[1:(length(subjects)-1)]))
      
      colnames(subjectmat) <- names(subjects)
      values <- data.frame(subjectmat,values)
      Std.Errors <- data.frame(subjectmat,Std.Errors)
    }
    ans <- list("models"=ans,"Estimates"=values,"Std.Errors"=Std.Errors)
    class(ans) <- "fittedlooplist" 
    attr(ans,"call") <- floopcall
    return(ans)
  }
 

  if (is.null(period))
    period <- length(dat$x)
 suppressWarnings(if (times=="equal")
  t <- (1:length(dat$x))/period*pi*2
 else if (is.numeric(times)) t <- 2*times/period*pi)
  if (method!="geometric") {
 matx <- cbind(rep(1,length(dat$x)),sin(t),cos(t))
 
 xfit <- lm.fit(matx,dat$x)
 cx <- as.vector(coef(xfit)[1])
 b.x <- as.vector(sqrt(coef(xfit)[2]^2+coef(xfit)[3]^2))
 phase.angle<- atan2(coef(xfit)[3],coef(xfit)[2])-pi/2
 costp <- cos(t+phase.angle)
 if (extended.classical==FALSE) maty <- cbind(rep(1,length(dat$x)),sin(t+phase.angle)^m,costp^n)
 if (extended.classical==TRUE) {
   direc <- sign(costp)
   maty <- cbind(rep(1,length(dat$x)),sin(t+phase.angle)^m,direc*abs(costp)^n)
 }
  yfit <- lm.fit(maty,dat$y)
 cy <- as.vector(coef(yfit)[1])
 retention <- as.vector(coef(yfit)[2])
 b.y <- as.vector(coef(yfit)[3])
 pred.x<-cx+b.x*cos(t+phase.angle)
  if (extended.classical==FALSE) pred.y<-cy+retention*sin(t+phase.angle)^m+b.y*costp^n
 if (extended.classical==TRUE)  pred.y<-cy+retention*sin(t+phase.angle)^m+direc*(b.y*abs(costp)^n)
fit <- list(xfit,yfit)
  } else {
   start <- direct(dat$x,dat$y) 
#Starting values taken as mean of those from direct ellipse fit, straight line from max to min.
#m and n start chosen by user.
if (times=="unknown") {
ti<-numeric(length(dat$x))
for (i in 1:length(dat$x)) {
  x0<-dat$x[i]
  y0<-dat$y[i]
  zmin1<-optimize(ellipsespot,c(0,pi),"x0"=x0,"y0"=y0,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"semi.major"=start$vals["semi.major"],"semi.minor"=start$vals["semi.minor"],"rote.rad"=start$vals["theta"])
  zmin2<-optimize(ellipsespot,c(pi,2*pi),"x0"=x0,"y0"=y0,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"semi.major"=start$vals["semi.major"],"semi.minor"=start$vals["semi.minor"],"rote.rad"=start$vals["theta"])
  ti[i]<-ifelse(zmin1$objective < zmin2$objective, zmin1, zmin2)[[1]]
}
ti<-c(ti[1],diff(ti))
}
else ti <- c(t[1],diff(t))
inti <- internal.1(start$vals["semi.major"],start$vals["semi.minor"],start$vals["theta"])
   mod=optim(par=c("t"=ti,"cx"=(start$vals["cx"]+mean(dat$x))/2,"cy"=(start$vals["cy"]+mean(dat$y))/2,"b.x"=(inti[1]+diff(range(dat$x))/2)/2,"b.y"=(inti[2]+diff(range(dat$y))/2)/2,"logm"=log(m),
                   "logn"=log(n),"retention"=inti[3]/2),fn=floopCauchyLoss,x=dat$x,y=dat$y,
             method="BFGS",hessian=TRUE)
   par = as.vector(mod$par)
   times <- par[1:length(x)]
   cx <- par[length(x)+1]
   cy <- par[length(x)+2]
   b.x <- par[length(x)+3]
   b.y <- par[length(x)+4]
   logm <- par[length(x)+5]
   logn <- par[length(x)+6]
   m <- exp(logm)
   n <- exp(logn)
   retention <- par[length(x)+7]
   t <- cumsum(times)
   phase.angle <- t[1]
   costp <- cos(t) 
   sintp <- sin(t) 
   t <- t - phase.angle
   direc <- sign(costp)
   direcsin <- sign(sintp)
   pred.x <- cx+b.x*costp
   pred.y <- cy+direcsin*retention*abs(sintp)^exp(logm)+direc*(b.y*abs(costp)^exp(logn))
   fit <- mod
 } 
  residuals <- sqrt((dat$x-pred.x)^2+(dat$y-pred.y)^2)
  if (n==1) beta.split.angle<-atan2(b.y,b.x)*180/pi 
    else if (n >= 2) beta.split.angle <- 0
    else beta.split.angle<-NA
    if (m==n) {
  hysteresis.x <- 1/sqrt(1+(b.y/abs(retention))^(2/m))
  coercion <- hysteresis.x*b.x
  } else {
  warning("hysteresis.x and coercion only available if m=n")
  hysteresis.x <- NA
  coercion <- NA
  }
  hysteresis.y <- retention/b.y
  area <- (0.5/(beta((m+3)/2,(m+3)/2)*(m+2))+1/beta((m+1)/2,(m+1)/2)-1/beta((m+3)/2,(m-1)/2))/(2^m)*pi*abs(retention*b.x)
if (method=="harmonic2" & m %% 2==0) area <- 0
 lag<-abs(atan2(retention,b.y))*period/(pi*2)
  ans <- list("values"=c("n"=n, "m"=m,"b.x"=b.x,"b.y"=b.y,"phase.angle"=as.vector(phase.angle),"cx"=cx,"cy"=cy,"retention"=retention,
               "coercion"=coercion,"area"=area, "lag"=lag,"beta.split.angle"=beta.split.angle,"hysteresis.x"=hysteresis.x, "hysteresis.y"=hysteresis.y),"fit"=fit,
              "x"=dat$x,"y"=dat$y,"pred.x"=pred.x,"pred.y"=pred.y,"period"=period, "period.time"=t+phase.angle,"residuals"=residuals,"call"=floopcall, "extended.classical"=extended.classical,"method"=method)
ans$call <- floopcall
 ans$Std.Errors <- try(unlist(delta.error.loop(ans)))
 ans$Estimates <- ans$values
 ans$Std.Errors <- ifelse(is.na(ans$Estimates),NA,ans$Std.Errors)
  class(ans) <- "fittedloop"
  ans
}
