plot.aCGHsegmented <-
function(x, add=FALSE, y=TRUE, psi.lines=TRUE, ...){
#se y=TRUE disegna anche le osservazioni e add è ignorato
        if(y) add<-FALSE
        arg<-list(...)
        if(is.null(arg$col)) arg$col=1
        if(is.null(arg$lwd)) arg$lwd=1
        if(is.null(arg$lty)) arg$lty=1
        if(is.null(arg$xlab)) arg$xlab<-"Genome Location"
        if(is.null(arg$ylab)) arg$ylab<-"LogRatio"
        
        yy<-x$y
        n<-length(yy)
        if(x$n.psi<=0) {
            if(y){
              plot(yy, ylab=arg$ylab, xlab=arg$xlab, pch=20, col=grey(.7))
              abline(h=x$est.means, ...)
              } else {
                if(add) {abline(h=x$est.means, ...)} else {
                  plot(yy, ylab=arg$ylab, xlab=arg$xlab, type="n")
                  abline(h=x$est.means, ...) }
                  }
          return(invisible(NULL))
          }
        #coll=arg$col
        psi<-x$psi
        est.means<-x$est.means
        if(y) {
            plot(yy, ylab=arg$ylab, xlab=arg$xlab, pch=20, col=grey(.7))
            points(c(1,psi,n), c(est.means,est.means[length(est.means)]),
              type="s",...)
              } else {
        if(!add) {
          plot(c(1,psi,n), c(est.means,est.means[length(est.means)]),
              type="s",ylab=arg$ylab, xlab=arg$xlab, col=arg$col, lwd=arg$lwd, lty=arg$lty) 
              } else {
          points(c(1,psi,n), c(est.means,est.means[length(est.means)]),
              type="s",col=arg$col, lwd=arg$lwd, lty=arg$lty)
              }
          }
          if(psi.lines){
            segments(x0=psi,y0=par()$usr[3],x1=psi,y1=est.means[-1],lty=3 ,col=arg$col)
            points(psi,rep(par()$usr[3],length(psi)),pch=19 ,col=arg$col)
            }
          invisible(NULL)
            }

