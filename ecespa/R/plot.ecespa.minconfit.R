`plot.ecespa.minconfit` <-
function(x, type="L", add=FALSE, xlim=NULL, ylim=NULL,
                            lwd=c(1,1),lty=c(1,2), col=c(1,2), main=NULL, ...){
   if(type!="L"){
       if(add==FALSE){
           plot(x$r, x$Kobs, xlab="r", ylab="K(r)", type="l", xlim=xlim, ylim=ylim,
                     lwd=lwd[1], lty =lty[1], col=col[1], ...)
           if(!is.null(main)) title(main=main) else title(main=x$dataname)
       }
       if(add==TRUE) lines(x$r, x$Kobs, lwd=lwd[1], lty =lty[1], col=col[1])
       lines(x$r, x$Kfit, lwd=lwd[2], lty =lty[2], col=col[2])
       print("dashed line is K fited")
    }
    else if(type =="L"){
            if(add==FALSE) {
            if(is.null(ylim)) ylim= range(c(range(sqrt(x$Kobs/pi)- x$r), range(sqrt(x$Kfit/pi)- x$r)))
            plot(x$r, sqrt(x$Kobs/pi)- x$r, xlab="r", ylab="L(r)" , 
                 type="l", xlim=xlim, ylim=ylim,
                 lwd=lwd[1], lty =lty[1], col=col[1], ...)
            if(!is.null(main)) title(main=main) else title(main=x$dataname)
         }
         if(add==TRUE) lines(x$r, sqrt(x$Kobs/pi)- x$r, lwd=lwd)
         lines(x$r, sqrt(x$Kfit/pi)- x$r, lwd=lwd[2], lty =lty[2], col=col[2])
         print("dashed line is L fited")
   }
    invisible(NULL)
}

