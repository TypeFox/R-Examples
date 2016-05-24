floop2r <- function(x,y=NULL,n=1,m=1,times="equal",period=NULL,subjects=NULL, subset=NULL,na.action=getOption("na.action"),extended.classical=FALSE,boot=FALSE,method="harmonic2",...) {
 if (boot==TRUE) return(summary(floop2r(x,y,n,m,times,period,subjects,subset,na.action,extended.classical,"method"=method),...))
  if (m==1 & n==1 & method=="harmonic2") return(fel(x,y,times=times,period=period,subjects=subjects,subset=subset,na.action=na.action,method="harmonic2"))
  floopcall <- match.call()
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
    class(ans) <- "fittedlooplist2r" 
    attr(ans,"call") <- floopcall
    return(ans)
  }

 
  if (is.null(period))
    period <- length(dat$x)
 suppressWarnings(if (times=="equal")
  t <- (0:(length(dat$x)-1))/period*pi*2
 else if (is.numeric(times)) t <- 2*times/period*pi)
 if (method=="harmonic2") {
 matx <- cbind(rep(1,length(dat$x)),sin(t),cos(t))
 
 xfit <- lm.fit(matx,dat$x)
 cx <- as.vector(coef(xfit)[1])
 b.x <- as.vector(sqrt(coef(xfit)[2]^2+coef(xfit)[3]^2))
 phase.angle <- atan2(coef(xfit)[3],coef(xfit)[2])-pi/2
 t2 <- t + phase.angle
 costp <- cos(t2)
 Ind <- (t2 < pi) & (t2 > 0) 
 if (extended.classical==FALSE) maty <- cbind(rep(1,length(dat$x)),sin(t2)^m,costp^n,Ind*sin(t2)^m)
 if (extended.classical==TRUE) {
   direc <- sign(costp)
   maty <- cbind(rep(1,length(dat$x)),sin(t2)^m,direc*abs(costp)^n,Ind*sin(t2)^m)
 }
  yfit <- lm.fit(maty,dat$y)
 cy <- as.vector(coef(yfit)[1])
 retention.below <- abs(as.vector(coef(yfit)[2]))
 retention.above <- abs(as.vector(coef(yfit)[2])+as.vector(coef(yfit)[4]))
 b.y <- as.vector(coef(yfit)[3])
 pred.x<-cx+b.x*cos(t2)
  if (extended.classical==FALSE) pred.y<-cy+(1-Ind)*retention.below*sin(t2)^m+Ind*retention.above*sin(t2)^m+b.y*costp^n
 if (extended.classical==TRUE)  pred.y<-cy+(1-Ind)*retention.below*sin(t2)^m+Ind*retention.above*sin(t2)^m+direc*(b.y*abs(costp)^n)
  fit <- list(xfit,yfit)
 }
 else {
 start <- direct(dat$x,dat$y) 
 if (times=="unknown") {
ti<-numeric(length(x))
for (i in 1:length(x)) {
  x0<-x[i]
  y0<-y[i]
  zmin1<-optimize(ellipsespot,c(0,pi),"x0"=x0,"y0"=y0,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"semi.major"=start$vals["semi.major"],"semi.minor"=start$vals["semi.minor"],"rote.rad"=start$vals["theta"])
  zmin2<-optimize(ellipsespot,c(pi,2*pi),"x0"=x0,"y0"=y0,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"semi.major"=start$vals["semi.major"],"semi.minor"=start$vals["semi.minor"],"rote.rad"=start$vals["theta"])
  ti[i]<-ifelse(zmin1$objective < zmin2$objective, zmin1, zmin2)[[1]]
}
ti<-c(ti[1],diff(ti))
}
else ti<-c(t[1],diff(t))
inti <- internal.1(start$vals["semi.major"],start$vals["semi.minor"],start$vals["theta"])
   mod=optim(par=c("t"=ti,"cx"=(start$vals["cx"]+mean(dat$x))/2,"cy"=(start$vals["cy"]+mean(dat$y))/2,"b.x"=(inti[1]+diff(range(dat$x))/2)/2,"b.y"=(inti[2]+diff(range(dat$y))/2)/2,"logm"=log(m),
                   "logn"=log(n),"retention.above"=inti[3]/2,"retention.below"=inti[3]/2),fn=floopCauchyLoss,x=dat$x,y=dat$y,
             method="BFGS",hessian=TRUE)
   par <- as.vector(mod$par)
   times <- par[1:length(x)]
   cx <- par[length(x)+1]
   cy <- par[length(x)+2]
   b.x <- par[length(x)+3]
   b.y <- par[length(x)+4]
   logm <- par[length(x)+5]
   logn <- par[length(x)+6]
   retention.above <- par[length(x)+7]
   retention.below <- par[length(x)+8]
   m <- exp(logm)
   n <- exp(logn)
   t <- cumsum(times)
   phase.angle <- t[1]
   t2 <- t - phase.angle
   costp <- cos(t) 
   sintp <- sin(t) 
   direc <- sign(costp)
   direcsin <- sign(sintp)
   pred.x <- cx+b.x*costp 
   pred.y <- cy+(direcsin < 0)*direcsin*retention.below*abs(sintp)^exp(logm)+(direcsin > 0)*direcsin*retention.above*abs(sintp)^exp(logm)+direc*(b.y*abs(costp)^exp(logn))
  fit <- mod 
 } 
  residuals <- sqrt((dat$x-pred.x)^2+(dat$y-pred.y)^2)
  if (n==1) beta.split.angle<-atan2(b.y,b.x)*180/pi 
    else if (n >= 2) beta.split.angle <- 0
    else beta.split.angle<-NA
       if (m==n) {
  hysteresis.x.above <- 1/sqrt(1+(b.y/retention.above)^(2/m))
 hysteresis.x.below <- 1/sqrt(1+(b.y/retention.below)^(2/m))
  coercion.above <- hysteresis.x.above*b.x
 coercion.below <- hysteresis.x.below*b.x
  } else {
  warning("hysteresis.x and coercion only available if m=n")
  hysteresis.x.above <- NA
 hysteresis.x.below <- NA
  coercion.above <- NA
 coercion.below <- NA
  }
  hysteresis.x.above <- 1/sqrt(1+(b.y/retention.above)^(2/m))
 hysteresis.x.below <- 1/sqrt(1+(b.y/retention.below)^(2/m))
  coercion.above <- hysteresis.x.above*b.x
 coercion.below <- hysteresis.x.below*b.x
  hysteresis.y.above <- retention.above/b.y
 hysteresis.y.below <- retention.below/b.y
  area <- (0.5/(beta((m+3)/2,(m+3)/2)*(m+2))+1/beta((m+1)/2,(m+1)/2)-1/beta((m+3)/2,(m-1)/2))/(2^m)*pi*abs((retention.above+retention.below)*b.x)/2
 if (method=="harmonic2" & m %% 2==0) area <- 0
 lag.above<-abs(atan2(retention.above,b.y))*period/(pi*2)
 lag.below<-abs(atan2(retention.below,b.y))*period/(pi*2)
  ans <- list("values"=c("n"=n, "m"=m,"b.x"=b.x,"b.y"=b.y,"phase.angle"=as.vector(phase.angle),"cx"=cx,"cy"=cy,"retention.above"=retention.above,
              "retention.below"=retention.below, "coercion.above"=coercion.above,"coercion.below"=coercion.below,"area"=area, "lag.above"=lag.above,"lag.below"=lag.below,"beta.split.angle"=beta.split.angle,
              "hysteresis.x.above"=hysteresis.x.above,"hysteresis.x.below"=hysteresis.x.below, "hysteresis.y.above"=hysteresis.y.above,"hysteresis.y.below"=hysteresis.y.below),"fit"=fit,
              "x"=dat$x,"y"=dat$y,"pred.x"=pred.x,"pred.y"=pred.y,"period"=period, "period.time"=t2,"residuals"=residuals,"call"=floopcall, "extended.classical"=extended.classical,"method"=method)
ans$call <- floopcall
 ans$Estimates <- ans$values
 ans$Std.Errors <- try(unlist(delta.error.loop2r(ans)))
  ans$Std.Errors <- ifelse(is.na(ans$Estimates),NA,ans$Std.Errors)
  class(ans) <- "loop2r"
  ans
}
