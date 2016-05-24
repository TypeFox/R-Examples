"plot.hdrconf" <- function(x, den, ...)
{
    ## Plots output from hdr.conf
    maxden <- max(den$y)
    xrange <- range(den$x)
    plot(den,ylim=c(-0.25*maxden,maxden),type="l",ylab="",xlab="",yaxt="n", ...)
    axis(2,at=pretty(c(0,maxden)),srt=90)
    abline(0,0)
    lines(xrange,rep(x$falpha,2),lty=4)
    lines(xrange,rep(x$falpha.ci[1],2),lty=2)
    lines(xrange,rep(x$falpha.ci[2],2),lty=2)

    add.hdr(x$hdr.hi, -0.16*maxden, 0.04*maxden, col=gray(0.9), horiz=TRUE)
    add.hdr(x$hdr,    -0.20*maxden, 0.04*maxden, col=gray(0.8), horiz=TRUE)
    add.hdr(x$hdr.lo, -0.24*maxden, 0.04*maxden, col=gray(0.9), horiz=TRUE)
}

"add.hdr" <- function(hdr, pos, width, col=2, horiz=FALSE, border=TRUE)
{
    ## Routine to add a single HDR on a graph
    ## Called by plot.hdr.conf().

    nint <- length(hdr[!is.na(hdr)])/2
    if(nint == 0)
        return(invisible())

    for(i in 1:nint)
    {
        l <- i*2-1  # lower
        tempx <- pos + c(-0.5,-0.5,0.5,0.5) * width
        tempy <- c(hdr[l],hdr[l+1],hdr[l+1],hdr[l])
        if(horiz==TRUE)
            polygon(tempy, tempx, col=col, border=border)
        else
            polygon(tempx, tempy, col = col, border=border)
    }
}
