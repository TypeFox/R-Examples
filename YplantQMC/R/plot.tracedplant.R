#'@method plot tracedplant
#'@S3method plot tracedplant
plot.tracedplant <- function (x, ...) 
{
	# 'Greens' from RColorBrewer
	PAL <- c("white", "#A1D99B", "#74C476", "#41AB5D", "#238B45", 
        "#005A32")
	r <- x$rays
	
    colIndex <- pmax(1, r$p + 1)
    colIndex <- pmin(colIndex, 6)
    r$plotcol <- PAL[colIndex]
    Xlim <- c(min(r$X), max(r$X))
    Ylim <- c(min(r$Y), max(r$Y))
    Ywid <- Ylim[2] - Ylim[1]
    Xwid <- Xlim[2] - Xlim[1]
    if (Ywid > Xwid) 
        Xlim <- Xlim * Ywid/Xwid
    else Ylim <- Ylim * Xwid/Ywid
    # with(r, plot(X, Y, pch = 15, cex = 0.3, col = plotcol, pty = "s", 
        # ylim = Ylim, xlim = Xlim))
    with(r, plot(X, Y, type='n', pty = "s", 
        ylim = Ylim, xlim = Xlim))	
	w2 <- x$gridpointwidth / 2
	for(i in 1:nrow(r)){
		X <- r$X[i]
		Y <- r$Y[i]
		rect(X-w2,Y-w2,X+w2,Y+w2,col=r$plotcol[i],border=NA)
	}		
}