curvetest.raw <-
function (fits1, fits2, equal.var,   
    conf.level , plotit) {   
    if(missing(plotit)) plotit=F
    if(missing(conf.level)) conf.level=0.05
    if(missing(equal.var)) equal.var=TRUE
    ww=fits1$kernel 
    myx=fits1$myx
    n1 <- length(fits1$data.model$x)    
    if (!missing(fits2)) {
        n2 <- length(fits2$data.model$x) 
        myx=sort(unique(c(myx, fits2$myx)))
    }
    nn<-length(myx)
    if (missing(fits2)) { 
       sx<-fits1$sx; y=fits1$data.model$y;res<-fits1$res; delta<-fits1$delta
       sigma.square <- sum(res^2)/delta; k0=fits1$k0; vv<-fits1$vv;
       CC <- max(abs(sx%*%y)/(sqrt(sigma.square) * sqrt(rowSums(sx^2))))
       p <- k0/3.1415926 * (1 + (CC * CC)/vv)^(-vv/2) + pt(-CC, vv) * 2
       out<-list(Statistic=CC, p=ifelse(p>1, 1, p),  eDF=vv, sigma.square= sigma.square, k0=k0, fits=fits1)
   } else  if (!missing(fits2))
     if(equal.var) { ##Fits2 and equal.var        
        delta11<-fits1$delta; delta12=fits1$delta2
        delta21<-fits2$delta; delta22<-fits2$delta2
        resid1 <- fits1$res
        resid2 <- fits2$res
        sx<-fits1$sx; y1<-fits1$data.model$y
        tx<-fits2$sx; y2<-fits2$data.model$y
        stx<-sqrt(rowSums(sx^2)+rowSums(tx^2))         
        T1X <- sweep(sx, 1, stx, FUN = "/")
        T2X <- sweep(tx, 1, stx, FUN = "/")
        k0 <- 0
        for (ii in 2:nn) k0 <- k0 + sum(sqrt(distance(T1X[ii, 
            ], T1X[ii - 1, ])^2 + distance(T2X[ii, ], T2X[ii-1, ])^2))
        vv<-(delta11+delta21)^2/(delta12+delta22)
        sigma.square <- (sum(resid1^2) + sum(resid2^2))/(delta11 +delta21)
        CC <- max(abs(sx %*% y1 - tx %*% y2)/(sqrt(sigma.square)*stx))
        p <- k0/3.1415926 * (1 + (CC * CC)/vv)^(-vv/2) + pt(-CC,  vv) * 2
        out<-list(Statistic=CC, p=ifelse(p>1, 1, p), eDF=vv, equal.var=equal.var, 
          sigma.square= sigma.square, k0=k0, fits1=fits1, fits2=fits2)
     }else if (!equal.var) {
        delta11<-fits1$delta; delta12=fits1$delta2; v1 <- delta11^2/delta12
        delta21<-fits2$delta; delta22<-fits2$delta2; v2 <- delta21^2/delta22
        resid1 <- fits1$res
        resid2 <- fits2$res
        esigma1 <- sqrt((sum(resid1^2))/delta11)
        esigma2 <- sqrt(sum(resid2^2)/delta21)
        sx<-fits1$sx; y1<-fits1$data.model$y
        tx<-fits2$sx; y2<-fits2$data.model$y
        vv <- (n1 + n2)^2 * v1 * v2/(n2^2 * v1 + n1^2 * v2)
        stx <- sqrt(esigma1^2 * rowSums(sx^2) + esigma2^2 * rowSums(tx^2))
        T1X <- esigma1 * sweep(sx, 1, stx, FUN = "/")
        T2X <- esigma2 * sweep(tx, 1, stx, FUN = "/")
        k0 <- 0
        for (ii in 2:nn) k0 <- k0 + sum(sqrt(distance(T1X[ii, 
            ], T1X[ii - 1, ])^2 + distance(T2X[ii, ], T2X[ii - 
            1, ])^2))
        CC <- max(abs(sx %*% y1 - tx %*% y2)/stx)
        p <- k0/3.1415926 * exp(-CC^2/2) + pnorm(-CC) * 2
        out<-list(Statistic=CC, p=ifelse(p>1, 1, p), eDF=vv, equal.var=equal.var, 
        esigma1=esigma1,esigma2=esigma2, k0=k0, fits1=fits1, fits2=fits2)
    }
    class(out)<-"curvetest"    
    invisible(out)
}
