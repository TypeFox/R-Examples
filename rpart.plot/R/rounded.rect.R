# rounded.rect.R

# similar to the standard function rect(), but round the rectangle corners
rounded.rect <- function(x1, y1, x2, y2, xlim, ylim, r, col, border.col, lty, lwd)
{
    if(all(r == 0)) { # all boxes have square corners? for efficiency treat specially
        # start at bottom-left corner, move clockwise, NA lifts pen
        x <- rbind(x1, x1, x2, x2, x1, NA)
        y <- rbind(y1, y2, y2, y1, y1, NA)
    } else {
        r <- recycle(r, x1)
        # adjustment xr to compensate for the aspect ratio of the figure
        xr <- r * (xlim[2] - xlim[1]) / (ylim[2] - ylim[1])
        # use pmin to prevent us drawing a quarter-circle bigger than avail space
        xr <- pmin(xr, abs(x2 - x1), abs(y2 - y1), na.rm=TRUE) / 2
        yr <- pmin(r, abs(x2 - x1), abs(y2 - y1), na.rm=TRUE) / 2
        theta <- seq(0, pi/2, length.out=7)
        xmat <- ymat <- matrix(nrow=length(x1), ncol=7)
        theta <- seq(0, pi/2, length.out=7)
        for(i in 1:7) {
            xmat[, i] <- xr * cos(theta[i])
            ymat[, i] <- yr * sin(theta[i])
        }
        x <- cbind( # TODO functionify these
            x1+xr-xmat[,1], x1+xr-xmat[,2], x1+xr-xmat[,3], x1+xr-xmat[,4],
                x1+xr-xmat[,5], x1+xr-xmat[,6], x1+xr-xmat[,7],
            x2-xr+xmat[,7], x2-xr+xmat[,6], x2-xr+xmat[,5], x2-xr+xmat[,4],
                x2-xr+xmat[,3], x2-xr+xmat[,2], x2-xr+xmat[,1],
            x2-xr+xmat[,1], x2-xr+xmat[,2], x2-xr+xmat[,3], x2-xr+xmat[,4],
                x2-xr+xmat[,5], x2-xr+xmat[,6], x2-xr+xmat[,7],
            x1+xr-xmat[,7], x1+xr-xmat[,6], x1+xr-xmat[,5], x1+xr-xmat[,4],
                x1+xr-xmat[,3], x1+xr-xmat[,2], x1+xr-xmat[,1],
            NA)
        x <- t(x)
        y <- rbind(
            y2-yr+ymat[,1], y2-yr+ymat[,2], y2-yr+ymat[,3], y2-yr+ymat[,4],
                y2-yr+ymat[,5], y2-yr+ymat[,6], y2-yr+ymat[,7],
            y2-yr+ymat[,7], y2-yr+ymat[,6], y2-yr+ymat[,5], y2-yr+ymat[,4],
                y2-yr+ymat[,3], y2-yr+ymat[,2], y2-yr+ymat[,1],
            y1+yr-ymat[,1], y1+yr-ymat[,2], y1+yr-ymat[,3], y1+yr-ymat[,4],
                y1+yr-ymat[,5], y1+yr-ymat[,6], y1+yr-ymat[,7],
            y1+yr-ymat[,7], y1+yr-ymat[,6], y1+yr-ymat[,5], y1+yr-ymat[,4],
                y1+yr-ymat[,3], y1+yr-ymat[,2], y1+yr-ymat[,1],
            NA)
    }
    # polygon doesn't recycle lwd the way we want, so use loop
    col        <- recycle(col, x1)
    border.col <- recycle(border.col, x1)
    lty        <- recycle(lty, x1)
    lwd        <- recycle(lwd, x1)
    for(i in 1:length(x1))
        polygon(x[,i], y[,i], col=col[i], border=border.col[i], lty=lty[i], lwd=lwd[i])
}
# similar to rounded.rect but blurs the edges of the box
draw.shadow <- function(x1, y1, x2, y2, xlim, ylim, r,
                        shadow.col, shadow.offset)
{
    # lighten color by alpha 0 ... 1 where 0 is white.
    lighten <- function(col, alpha)
    {
        if(is.character(col)) {
            # needed when elems of col were set to 0 although col has type character
            col[col=="0"] = par("bg")
        }
        rgb <- col2rgb(col) / 255
        if(device.supports.alpha)
            rgb(rgb[1,], rgb[2,], rgb[3,], alpha=alpha)
        else {
            rgb <- rgb + (1 - alpha) * (c(1,1,1) - rgb) # move each r,g,b towards 1
            rgb(rgb[1,], rgb[2,], rgb[3,])
        }
    }
    #--- draw.shadow begins here ---

    # Guess if device supports an alpha channel.  We want to use the alpha
    # channel where possible because overlapping shadows look better with it.
    # TODO make this more accurate i.e. check for devices other than postscript
    #      is there no way of querying a device's capabilities?

    device.supports.alpha <- names(dev.cur())[1] != "postscript"

    if(device.supports.alpha)
        alphas <- c(1, 0.64, 0.3, 0.2, 0.1)
    else
        alphas <- c(1, 0.8, 0.5, 0.3, 0.1)

    strwidth1 <- my.strwidth("M")
    width <- shadow.offset * strwidth1
    blur.width <- width / 8
    for(i in 5:1)
        rounded.rect(x1 + width - i * blur.width,
                     y1 - width - i * blur.width,
                     x2 + width + i * blur.width,
                     y2 - width + i * blur.width,
                     xlim, ylim, r,
                     lighten(shadow.col, alphas[i]),
                     0, 0, 1)
}
