p.synth.peak.error <- function(peaks, y.max=(max(peaks,na.rm=TRUE)),
peak.cluster = NULL, peak.palette=grey(c(0,0.6,0.8)), use.layout=TRUE,
show.errors=1:n.errors, peak.lty=rep(1,n.errors),
mfrow=c(2, ceiling(length(show.errors)/2)), plot.legend=TRUE, print.error.nr=TRUE){
    n.errors=dim(peaks)[2]
    n.levels=dim(peaks)[3]
    old.pal <- palette(peak.palette)
    if(use.layout){
        nf <- layout(matrix(1:10,2,5,byrow=TRUE), c(1.5, 1,1,1,1), c(1,1.55), TRUE)
        par(oma=c(1,1,0,0)+0.1)
    } else if(!is.null(mfrow)){
        par(mfrow = mfrow)
        par(oma=c(1,1,0,0)+0.1)
    }
    if(use.outer <- !is.null(mfrow)){
        xlab <- ""
        ylab <- ""
    } else {
        ylab <- "specific discharge/mm/h"
        xlab <- "time/h"
    }
    #layout.show(nf)
    for(error in show.errors){
        if(use.layout){
           if(error > 5) {
               b.mar <-4
               x.axt = "s"
           } else {
               b.mar <-0
               x.axt = "n"
           }
           if(error%%5==1){
               l.mar <-4
               y.axt = "s"
           } else {
               l.mar <-0
               y.axt = "n"
           }
           par(mar=c(b.mar,l.mar,0,0)+0.1 )
        } else {
           par(mar=c(4,4,0,0)+0.1 )
           x.axt = "s"
           y.axt = "s"
        }
       plot(peaks[2,error,1,],type="n", xlab=xlab,ylab=ylab, xaxt=x.axt, yaxt=y.axt, ylim=c(0,y.max))
       if(print.error.nr){
	       text(x=150, y=0.8*y.max, error, cex=2)
       }
       for(level in 1:n.levels){
           if(is.null(peak.cluster)){
                peak.col <- ((level-1)%/%3)+2
                t.peak.lty <- peak.lty[((level-1)%/%3)+1]
           } else {
                peak.col <- peak.cluster[error,level]
                t.peak.lty <- peak.lty[peak.cluster[error,level]]
           }
           lines(peaks[1,error,level,], col=peak.col, lty=t.peak.lty)
           lines(peaks[2,error,level,])
       }
    }
    if(use.outer){
        mtext(outer=TRUE, side=2, line=-1, text="specific discharge/mm/h")
        mtext(outer=TRUE, side=1, line=-1, text="time/h")
    }

    if(!is.null(peak.cluster) & plot.legend){
	    plot.new()
            legend("left", legend=c("reference", paste("Cluster", LETTERS[1:max(peak.cluster)])), lty=c(1,peak.lty), col=c("black", peak.palette))
    }
    
    palette(old.pal)
}

