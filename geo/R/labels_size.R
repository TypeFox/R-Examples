#' Label symbols of given size
#' 
#' Labels of given size.
#' 
#' 
#' @param cont Contours
#' @param digits Number of digits to use in labels
#' @param sizes Sizes (of what?)?
#' @param xlim,ylim Limits
#' @param fill Fill? Default FALSE
#' @param n Number of ??
#' @param rat Ratio of ??
#' @param minsym Minimum symbol for lowest category in labels
#' @param label.resolution Label resolution ?
#' @param open Open legend/label, default FALSE
#' @param lwd Line width
#' @param col Color
#' @return No value, labels added to current plot.
#' @note Needs further elaboration, document with other labelling functions??
#' @seealso Called by \code{\link{colsymbol}}.
#' @keywords aplot
#' @export labels_size
labels_size <-
function(cont, digits, sizes, xlim = c(0, 1), ylim = c(0, 1), fill = F, n,
        rat, minsym = "<", label.resolution = 0, open = F, lwd = 1, col = 1)
{
        xlim <- sort(xlim)
        ylim <- sort(ylim)
        ncont <- length(cont)
        lbox <- ncont + 1
        if(fill)
                lbox <- max(lbox, 20)
        boxy <- c(1:lbox)
        boxy <-  - boxy/lbox + 1
        boxy1 <- boxy + 1/(1.2 * lbox)
        if(fill) {
                boxy <- boxy[1:(ncont + 1)]
                boxy1 <- boxy1[1:(ncont + 1)]
        }
        yloc <- (boxy + boxy1)/2
        xloc <- matrix(0.85, length(yloc))
        theta <- (c(0:n) * 2 * pi)/n
        theta <- c(theta, NA)
        theta <- c(matrix(theta, n + 2, length(yloc)))
        par(adj = 0)
        cont <- round(cont, digits = digits)
        textx <- c(1:(length(cont) - 1))
        textx1 <- textx
        textx <- format(round(cont[1:(length(cont) - 1)] + label.resolution,
                digits = digits))
        textx1 <- format(round(cont[2:length(cont)], digits = digits))
        textx <- paste(textx, "-", textx1)
        minsym <- paste(minsym, " ", sep = "")
        tmp1 <- paste(minsym, format(round(cont[1], digits = digits)))
        tmp2 <- paste("> ", format(round(cont[ncont], digits = digits)))
        textx <- c(tmp1, textx, tmp2)
        boxx <- c(matrix(0.1, 1, length(boxy)))
        boxx <- xlim[1] + abs((xlim[2] - xlim[1])) * boxx
        xloc <- xlim[1] + abs((xlim[2] - xlim[1])) * xloc
        yloc <- ylim[1] + abs((ylim[2] - ylim[1])) * yloc
        boxy <- ylim[1] + (ylim[2] - ylim[1]) * boxy
        ll <- (ylim[2] - ylim[1]) * 0.05
        # put the labels.
        if(fill) text(boxx, boxy + ll/2, textx) else text(boxx, boxy + ll,
                        textx)
        # put the labels.
        theta <- (c(0:n) * 2 * pi)/n
        theta <- c(theta, NA)
        theta <- c(matrix(theta, n + 2, length(boxy)))
        y <- c(t(matrix(yloc, length(yloc), n + 2)))
        x <- c(t(matrix(xloc, length(xloc), n + 2)))
        sizes <- c(t(matrix(sizes, length(boxx), n + 2)))
        y <- y + sizes * rat * sin(theta)
        x <- x + sizes * rat * cos(theta)
        if(!open)
                polygon(x, y, col = col, border = T)
        else lines(x, y, col = col, lwd = lwd)
}

