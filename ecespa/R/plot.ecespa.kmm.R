`plot.ecespa.kmm` <-
function(x,type="Kmm.n", q=0.025, 
            xlime=NULL, ylime=NULL,  maine=NULL, add=FALSE, kmean=TRUE,
            ylabe=NULL, xlabe=NULL, lty=c(1,2,3), col=c(1,2,3), lwd=c(1,1,1),
             ...)
{
   if(!is.null(x$nsim)){
if (type == "Kmm.n") {        
    cosa <- cbind(x$kmm.n, 
                  t(apply(x$kmmsim.n, 1, quantile, c(q,1-q), na.rm=TRUE)))
            cosamean <- t(apply(x$kmmsim.n, 1, mean, na.rm=TRUE))
    if(is.null(ylabe)) ylabe <- expression(normalized.K[mm](r)) 
            
}
        if (type == "Kmm") {        
    cosa <- cbind(x$kmm, 
                  t(apply(x$kmmsim, 1, quantile, c(q,1-q), na.rm=TRUE)))
            cosamean <- t(apply(x$kmmsim, 1, mean, na.rm=TRUE))
            if(is.null(ylabe)) ylabe <- expression(K[mm](r))
            

}
       if(is.null(maine)) maine <- x$dataname
       if(is.null(xlabe)) xlabe <- "distance"
       matplot(x$r, cosa, type="l", lty=c(lty[1],lty[2],lty[2]), 
           col=c(col[1],col[2],col[2]), lwd=c(lwd[1],lwd[2],lwd[2]),
           main= maine,xlab= xlabe, ylab= ylabe, ylim=ylime, xlim=xlime,
           mgp=c(2.5,1,0), add=add)
        if (kmean==TRUE) lines(x$r, cosamean, lty=lty[3], col=col[3],lwd=lwd[3])
   }
   
   if(is.null(x$nsim)){
     if (type == "Kmm.n") {        
    cosa <- x$kmm.n 
    if(is.null(ylabe)) ylabe <- expression(normalized.K[mm](r)) 
            
}
        if (type == "Kmm") {        
    cosa <- x$kmm
     if(is.null(ylabe)) ylabe <- expression(K[mm](r))
}
       if(is.null(maine)) maine <- x$dataname
       if(is.null(xlabe)) xlabe <- "distance"
       matplot(x$r, cosa, type="l", lty=lty[1], 
           col=col[1], lwd=lwd[1],
           main= maine,xlab= xlabe, ylab= ylabe, ylim=ylime, xlim=xlime,
           mgp=c(2.5,1,0), add=add)
        
   }
   
}

