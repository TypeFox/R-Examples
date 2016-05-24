splitloopboot <-
  function(j=NULL,pred.x,pred.y,xresid,yresid,ti,obs,n,m,extended.classical,cbb,joint,period){
    if (is.numeric(cbb)==TRUE) {
      index <- sample(1:(obs+3),3,replace=FALSE)
      xresid2 <- c(xresid,xresid)
      yresid2 <- c(yresid,yresid)
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
    t2 <- ti[-index] + phase.angle
    cx<-coef(Ta.lm)[[1]]
    costp <- cos(t2)
    Ind <- (t2 < pi) & (t2 > 0) 
    if (extended.classical==FALSE) maty <- cbind(rep(1,length(t2)),sin(t2)^m,costp^n,Ind*sin(t2)^m)
    if (extended.classical==TRUE) {
      direc <- sign(costp)
      maty <- cbind(rep(1,length(t2)),sin(t2)^m,direc*abs(costp)^n,Ind*sin(t2)^m)
    }
    yfit <- lm.fit(maty,y)
    cy <- as.vector(coef(yfit)[1])
    retention.below <- abs(as.vector(coef(yfit)[2]))
    retention.above <- abs(as.vector(coef(yfit)[2])+as.vector(coef(yfit)[4]))
    b.y <- as.vector(coef(yfit)[3])
    if (n==1) beta.split.angle<-atan2(b.y,b.x)*180/pi 
    else if (n >= 2) beta.split.angle <- 0
    else beta.split.angle<-NA
    hysteresis.x.above <- 1/sqrt(1+(b.y/retention.above)^(2/m))
    hysteresis.x.below <- 1/sqrt(1+(b.y/retention.below)^(2/m))
    coercion.above <- hysteresis.x.above*b.x
    coercion.below <- hysteresis.x.below*b.x
    hysteresis.y.above <- retention.above/b.y
    hysteresis.y.below <- retention.below/b.y
    area <- (0.5/(beta((m+3)/2,(m+3)/2)*(m+2))+1/beta((m+1)/2,(m+1)/2)-1/beta((m+3)/2,(m-1)/2))/(2^m)*pi*abs((retention.above+retention.below)*b.x)/2
    lag.above<-abs(atan2(retention.above,b.y))*period/(pi*2)
    lag.below<-abs(atan2(retention.below,b.y))*period/(pi*2)
    z <- c("n"=n, "m"=m,"b.x"=b.x,"b.y"=b.y,"phase.angle"=as.vector(phase.angle),"cx"=cx,"cy"=cy,"retention.above"=retention.above,
                           "retention.below"=retention.below, "coercion.above"=coercion.above,"coercion.below"=coercion.below,"area"=area, "lag.above"=lag.above,"lag.below"=lag.below,"beta.split.angle"=beta.split.angle,
                           "hysteresis.x.above"=hysteresis.x.above,"hysteresis.x.below"=hysteresis.x.below, "hysteresis.y.above"=hysteresis.y.above,"hysteresis.y.below"=hysteresis.y.below)
    z
  }
