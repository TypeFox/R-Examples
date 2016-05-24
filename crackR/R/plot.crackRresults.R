plot.crackRresults <-
function(x, boot=FALSE, thresh=NA, sfpof.int=FALSE, ...)
{

    obj <- x
    rm(x)

    ## for sfpof.int==TRUE, use pof.int to approx SFPOF in each interval
    if( sfpof.int ) obj$sfpof <- calcSfpofFromPofInt(obj$pof.int)

    ## if sfpof has zero values, we need to provide a small positive number for plotting
    if(min(obj$sfpof$sfpof)<=0) obj$sfpof$sfpof[obj$sfpof$sfpof <= 0] <- min(obj$sfpof$sfpof[obj$sfpof$sfpof>0])
    
    ## using base R graphics

    temp.mar <- par("mar")
    par(mar=c(5.1,6.1,4.1,2.1))
    plot(x=obj$sfpof$flight, y=obj$sfpof$sfpof, log="y", xlab=NA, ylab=NA,
         axes=FALSE, type="l", ...)
    box()
    axis(side=1)
    axis(side=2, las=1, at=10^(seq(from=-40, to=0, by=2)),
         labels=expression(10^40,10^-38,10^-36,10^-34,10^-32,
             10^-30,10^-28,10^-26,10^-24,10^-22,
             10^-20,10^-18,10^-16,10^-14,10^-12,
             10^-10,10^-08,10^-06,10^-04,10^-02,10^00))
    
    mtext(side = 1, "Flight Number", line = 2.5)
    mtext(side = 2, "SFPOF Estimate", line = 3.5)

    ## plot bootstrap lines as well
    if( boot )
        {
            columns.in.sfpof <- dim(obj$sfpof)[2]
            if( columns.in.sfpof > 2 )
                {
                    for(bbb in 3:columns.in.sfpof)
                        lines(x=obj$sfpof$flight, y=obj$sfpof[,bbb], col=4, lty=2, ...)
                }
        }

    if( !is.na(thresh) )
        {
            abline(h=thresh, lty=3, lwd=2, col="red")
            mtext(paste("SFPOF threshold =", thresh),
                  side=3, adj=1, cex=0.7, line=0.5, col="red")
        }

    par(mar=temp.mar)
    
}
