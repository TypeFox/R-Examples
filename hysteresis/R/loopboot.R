loopboot <-
  function(j=NULL,pred.x,pred.y,xresid,yresid,ti,obs,n,m,extended.classical,cbb,joint,period){
    if (is.numeric(cbb)==TRUE) {
      index <- sample(1:(obs+3),3,replace=FALSE)
      xresid2 <- c(xresid,xresid)
      yresid2 <- c(xresid,xresid)
      k <- obs/cbb
      xblocks <- sample(1:(obs+3),k,replace=TRUE)
      if (joint==FALSE) yblocks <- sample(1:(obs+3),k,replace=TRUE)
      else yblocks <- xblocks
      xressamp <- c(t(outer(xblocks,0:(cbb-1),FUN="+")))
      yressamp <- c(t(outer(yblocks,0:(cbb-1),FUN="+")))
      y<-yresid2[yressamp]+pred.y[-index]
      x<-xresid2[xressamp]+pred.x[-index]
    }
    else {  
      index <- sample(1:(obs+3),3,replace=FALSE)
      if (joint==FALSE) {
    y<-sample(yresid,obs,replace=T)+pred.y[-index]
    x<-sample(xresid,obs,replace=T)+pred.x[-index]
      }
      else {
        resid.sampler <- sample(1:(obs+3),obs,replace=TRUE)
        y<-yresid[resid.sampler]+pred.y[-index]
        x<-xresid[resid.sampler]+pred.x[-index]
      }
    }
    
    Ta.lm<-lm.fit(cbind(rep(1,obs),sin(ti[-index]),cos(ti[-index])),x)            
    b.x<-sqrt(coef(Ta.lm)[[2]]^2+coef(Ta.lm)[[3]]^2)
    phase.angle<- atan2(coef(Ta.lm)[[3]],coef(Ta.lm)[[2]])-pi/2
    rad<-ti[-index]+phase.angle
    cx<-coef(Ta.lm)[[1]]
    if (extended.classical==FALSE) {
    Tb.lm<-lm.fit(cbind(rep(1,obs),sin(rad)^m,cos(rad)^n),y)}
    else {
      direc <- sign(cos(rad))
      Tb.lm<-lm.fit(cbind(rep(1,obs),sin(rad)^m,direc*abs(cos(rad))^n),y)}
    
    b.y<-coef(Tb.lm)[[3]]
    retention<- coef(Tb.lm)[[2]]
    cy<-coef(Tb.lm)[[1]]
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
