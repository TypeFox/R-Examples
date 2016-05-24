##"shift" factor to shift away from the rect legend default=1
metre <- function(xl, yb, xr, yt, lab, cols = risk.colors(length(lab) - 1),
                  shift = 0, cex = 1)
{
    n <- length(lab)-1
    dx <- xr-xl
    dy <- yt-yb
    dxy <- max(dx, dy)/n ##increasing step
    drift <- min(xr-xl, yt-yb)*shift
    if(dx>dy) {
        rect(xl+(1:n-1)*dxy, yb, xl+(1:n)*dxy, yt, col=cols)
        text(xl+(0:n)*dxy, yb-drift, lab, cex=cex, pos=1)
    } else {
        rect(xl, yb+(1:n-1)*dxy, xr, yb+(1:n)*dxy, col=cols)
        text(xr+drift, yb+(0:n)*dxy, lab, cex=cex, pos=4)
    }
	return()
}
