loopbootgeom <-
  function(j=NULL,pred.x,pred.y,xresid,yresid,obs,extended.classical,cbb,joint,period){
    if (is.numeric(cbb)==TRUE) {
    index2 <- sample(1:(obs+3),5,replace=FALSE)
      index <- index2[1:3]
      xresid2 <- c(xresid,xresid)
      yresid2 <- c(yresid,yresid)
      k <- round(obs/cbb)
      k2 <- round((obs-2)/cbb) #Round so not exact if obs/cbb and obs-2/cbb aren't integers.
      xblocks <- sample(1:(obs+3),k,replace=TRUE)
      if (joint==FALSE) yblocks <- sample(1:(obs+3),k2,replace=TRUE)
      else yblocks <- xblocks[1:k2]
      xressamp <- c(t(outer(xblocks,0:(cbb-1),FUN="+")))
      yressamp <- c(t(outer(yblocks,0:(cbb-1),FUN="+")))
      y<-yresid2[yressamp]+pred.y[-index2]
      x<-xresid2[xressamp]+pred.x[-index]
    }
    else {  
    index2 <- sample(1:(obs+3),5,replace=FALSE)
      index <- index2[1:3]
      if (joint==FALSE) {
    y<-sample(yresid,obs-2,replace=TRUE)+pred.y[-index2]
    x<-sample(xresid,obs,replace=TRUE)+pred.x[-index]
      }
      else {
        resid.sampler <- sample(1:(obs+3),obs,replace=TRUE)
        y<-yresid[resid.sampler[1:(obs-2)]]+pred.y[-index2]
        x<-xresid[resid.sampler]+pred.x[-index]
      }
    }
    
    
 start <- direct(x[1:length(y)],y) 

ti<-rep(2*pi/length(x),length(x))
inti <- internal.1(start$vals["semi.major"],start$vals["semi.minor"],start$vals["theta"])
   mod=optim(par=c("t"=ti,"cx"=start$vals["cx"],"cy"=start$vals["cy"],"b.x"=inti[1],"b.y"=inti[2],"logm"=0,
                   "logn"=0,"retention"=inti[3]/2),fn=floopCauchyLoss,x=x,y=y,
             method="BFGS",hessian=TRUE)
    par = mod$par
    times <- par[1:length(x)]
    cx <- par[length(x)+1]
    cy <- par[length(x)+2]
    b.x <- par[length(x)+3]
    b.y <- par[length(x)+4]
    logm <- par[length(x)+5]
    m <- exp(logm)
    logn <- par[length(x)+6]
    n <- exp(logn)
    retention <- par[length(x)+7]
    t <- cumsum(times)
    phase.angle <- t[1]
    if (n==1) beta.split.angle<-atan2(b.y,b.x)*180/pi 
    else if (n >= 2) beta.split.angle <- 0
    else beta.split.angle<-NA
    hysteresis.x <- 1/sqrt(1+(b.y/retention)^(2/m))
    coercion <- hysteresis.x*b.x
    hysteresis.y <- retention/b.y
    lag<-abs(atan2(retention,b.y))*period/(pi*2)
    area <- (0.5/(beta((m+3)/2,(m+3)/2)*(m+2))+1/beta((m+1)/2,(m+1)/2)-1/beta((m+3)/2,(m-1)/2))/(2^m)*pi*abs(retention*b.x)
    
    z <- c("n"=n, "m"=m,"b.x"=b.x,"b.y"=b.y,"phase.angle"=as.vector(phase.angle),"cx"=cx,"cy"=cy,"retention"=retention,
           "coercion"=coercion,"area"=area,"lag"=lag, "beta.split.angle"=beta.split.angle,"hysteresis.x"=hysteresis.x, "hysteresis.y"=hysteresis.y)
z
  }
