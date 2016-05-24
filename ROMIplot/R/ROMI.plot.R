ROMI.plot <-
function(Dx=NULL, Nx=NULL, mx=NULL, smooth=TRUE) {
    if (all(is.null(Dx), is.null(Nx), is.null(mx)))
        stop("You need to specify either 'Dx' and 'Nx' or, alternatively, 'mx'!")
    if (!smooth) {
        if (is.null(mx)) {
            mx <- Dx/Nx
        }
        mx.2 <- mx[,-1]
        mx.1 <- mx[,-(ncol(mx))]            
        aai <- 100 * -log(mx.2/mx.1)
    } else {
        if (any(is.null(Dx), is.null(Nx))) 
            stop("If you want smooth your data, You need to specify 'Dx' AND 'Nx'!")
        if (any(Nx==0))
            stop("The current version does not support zero exposures.")
        the.agesDx <- as.numeric(rownames(Dx))
        the.yearsDx <- as.numeric(colnames(Dx))
        the.agesNx <- as.numeric(rownames(Nx))
        the.yearsNx <- as.numeric(colnames(Nx))
        if ((!identical(the.agesDx, the.agesNx)) || 
            (!identical(the.yearsDx, the.yearsNx)))
            stop("It seems you use different dimensions (=> Years or Ages).")

        fitDx <- Mort2Dsmooth(x=the.agesDx, y=the.yearsDx,
                              Z=Dx, offset=log(Nx))
        mx <- fitDx$fitted.values / Nx
        mx.2 <- mx[,-1]
        mx.1 <- mx[,-(ncol(mx))]            
        aai <- 100 * -log(mx.2/mx.1)
    }

    mylevels <- c(-10000, -5, -3, -1, seq(-0.5,5,by=0.5), 6, 10000)
    mycolors <- c(grey(seq(0.1,0.7, length=4)),  rep("white",2),
                  hsv(h=seq(from=240, to=180, length=3)/360, s=1, v=1),
                  hsv(h=120/360, s=1, v=seq(0.5,1,length=4)),
                  hsv(h=seq(from=0, to=60, length=4)/360, s=1, v=1))


    layout.matrix <- matrix(c(rep(1,8),2,2), nrow=1)
    the.layout <- layout(layout.matrix)
    par(mar=c(5,4,4,0)+ 0.1)
    image(x=the.yearsDx, y=the.agesDx, z=t(aai),
          xlim=range(the.yearsDx), ylim=range(the.agesDx),
          col=mycolors, breaks=mylevels, xlab="Year", ylab="Age")
    grid(lty=1, col="white")

    par(mar=c(5,0,4,2)+ 0.1)    
    plot(x=1,y=1, type="n", xlab="", ylab="", axes=FALSE,
         xlim=c(0,10), ylim=range(the.agesDx))
    yb <- min(the.agesDx)
    yt <- min(the.agesDx) + 0.9 * (max(the.agesDx) - min(the.agesDx))
    ycoords <- seq(from=yb, to=yt, length=length(mycolors)+1)
    rect(xleft=0, xright=5,
         ybottom=ycoords[-(length(ycoords))], ytop=ycoords[-1],
         col=mycolors)
    keylabels <- c("", "-5.0", "-3.0", "-1.0", "-0.5", "0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0",
                   "3.5", "4.0", "4.5", "5.0", "6.0", "")
    text(x=8, y=ycoords, labels=keylabels, pos=2)
    segments(x0=0, x1=5.5,
             y0=ycoords, y1=ycoords,
             col="black", lty=1)
    text(x=5, y=yt, labels=expression(paste(rho, " (in %)")), pos=3, cex=1.25)
    return(aai)
}
