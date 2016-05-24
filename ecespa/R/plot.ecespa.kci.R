`plot.ecespa.kci` <-
function(x, type=1, q=0.025, kmean=TRUE, add=FALSE,
                             maine=NULL, xlabe=NULL, ylabe=NULL, xlime=NULL, ylime=NULL,    
     lty=c(1,2,3), col=c(1,2,3), lwd=c(1,1,1), ...)
{
  
   if (type == 1) {                      ## univariado t<ed>pico L1
       cosa <- cbind(x$kia, 
                  t(apply(x$kia.s,1,quantile, c(q,1-q),na.rm=TRUE)))
       cosamean <- sqrt(apply(x$kia.s,1,mean,na.rm=TRUE)/pi)- x$r
       cosa <- sqrt(cosa/pi)-x$r
       if(is.null(ylabe)) ylabe <- expression(hat(L)[1])
       if(is.null(maine)) maine <- x$datanamea
   }
   if (type == 2) {                      ## univariado t<ed>pico L2
       cosa <- cbind(x$kib, 
                  t(apply(x$kib.s,1,quantile, c(q,1-q),na.rm=TRUE)))
       cosamean <- sqrt(apply(x$kib.s,1,mean,na.rm=TRUE)/pi)- x$r
       cosa <- sqrt(cosa/pi)-x$r
       if(is.null(ylabe)) ylabe <- expression(hat(L)[2])
       if(is.null(maine)) maine <- x$datanameb
   }
   if (type == 12) {                      ## bivariado t<ed>pico L12
       cosa <- cbind(x$kci.ab.o, 
                  t(apply(x$kci.ab.s,1,quantile, c(q,1-q),na.rm=TRUE)))
       cosa <- sqrt(cosa/pi)-x$r
       cosamean <- apply(x$kci.ab.s,1,mean,na.rm=TRUE)
       cosamean <- sqrt(cosamean/pi)-x$r
       if(is.null(ylabe)) ylabe <- expression(hat(L)[12])
       if(is.null(maine)) maine <- paste(x$datanamea, " vs. ", x$datanameb) 
   }
   
   if (type == 21) {                       ## bivariado t<ed>pico L21
       cosa <- cbind(x$kci.ba.o, 
                  t(apply(x$kci.ba.s,1,quantile, c(q,1-q),na.rm=TRUE)))
       cosa=sqrt(cosa/pi)-x$r
       cosamean <- apply(x$kci.ba.s,1,mean,na.rm=TRUE)
       cosamean <- sqrt(cosamean/pi)-x$r

       if(is.null(ylabe)) ylabe <- expression(hat(L)[21])
       if(is.null(maine)) maine <- paste(x$datanamea, " vs. ", x$datanameb) 
  }

   if (type == 112) {                       ## segregaci<f3>n primer ppp (K1-K12)
       d1.12.o <- x$kia - x$kci.ab.o
       d1.12.s <- x$kia - x$kci.ab.s
       cosa <- cbind(d1.12.o, t(apply(d1.12.s,1,quantile, c(q,1-q),na.rm=TRUE)))
       cosamean <- apply(d1.12.s,1,mean,na.rm=TRUE)
       if(is.null(ylabe)) ylabe <- expression(hat(K)[1]- hat(K)[12])
       if(is.null(maine)) maine <- paste(x$datanamea, " vs. ", x$datanameb) 
   }

   if (type == 221) {                       ## segregaci<f3>n segundo ppp (K2-K21)
       d2.21.o = x$kib - x$kci.ba.o
       d2.21.s = x$kib.s - x$kci.ba.s
       cosa <- cbind(d2.21.o, t(apply(d2.21.s,1,quantile, c(q,1-q),na.rm=TRUE)))
       cosamean <- apply(d2.21.s,1,mean,na.rm=TRUE)
       if(is.null(ylabe)) ylabe <- expression(hat(K)[2]- hat(K)[21])
       if(is.null(maine)) maine <- paste(x$datanamea, " vs. ", x$datanameb) 

   }
   if(is.null(xlabe)) xlabe <-"distance"
   matplot(x$r,cosa, type="l", lty=c(lty[1],lty[2],lty[2]), 
           col=c(col[1],col[2],col[2]), lwd=c(lwd[1],lwd[2],lwd[2]),
           main= maine,xlab= xlabe, ylab= ylabe, ylim=ylime, xlim=xlime,
           mgp=c(2.5,1,0), add=add)
  if (kmean==TRUE) lines(x$r, cosamean, lty=lty[3], col=col[3],lwd=lwd[3])

}

