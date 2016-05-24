#' Compute letter value summary table.
#'
#' @param x numeric vector
#' @param qu quantiles to compute
#' @param out binary vector of outliers (\code{TRUE} for outlier,
#'   \code{FALSE} otherwise)
#' @inheritParams determineDepth
#' @keywords internal
#' @return
#'  \item{letter.val}{letter value statistic, distinguishes between upper and
#'    lower LV statistic for all statistics but the median}
#'  \item{conf.int}{confidence interval of corresponding letter value
#'    statistic}
#'  \item{out}{list of defined outliers}
#' @noRd
outputLVplot <- function(x,qu,k,out,alpha) {
  n <- length(x)

  depth <- getDepth(k,n)

  LV <- cbind(depth,lower=qu[k:1],upper=qu[k-1+1:k])
  conf <- confintLV(x, k, alpha=alpha)

  dimnames <- nameLV(k)
  row.names(LV) <- dimnames[[1]]
  row.names(conf) <- dimnames[[2]]

  result <- list(letter.val = LV, conf.int= conf,outliers = which(out))
  return(result)
}

#' Draw an LV plot.
#'
#' @param x x positions
#' @param y y positions
#' @param k number of letter value statistics used
#' @param out indices of outliers
#' @param qu quantiles
#' @param horizontal display horizontally (TRUE) or vertically (FALSE)
#' @param col vector of colours to use
#' @param border vector of colours to use
#' @param width maximum height/width of box
#' @param width.method one of 'linear', 'height' or 'area'. Methods 'height' and 'area' ensure that these dimension are proportional to the number of observations within each box.
#' @keywords internal
drawLVplot <- function(x,y,k,out,qu,horizontal,col, border="grey50", width=0.9, width.method = "linear", median.col, ...) {
  i <- seq_len(k)
  y <- rep(y, length(x))

  lower <- i[-k]
  upper <- rev(seq_len(k-1) + k)
  col <- rev(col[i])

  if (width.method=="linear") {
    offset <- width/2* c(lower / k, 1, rev(lower) / k)
  } else {
    if (width.method=="area") {
      areas <- .5*c(2^(lower-k), rev(2^(lower-k)))
      height <- diff(qu)
      # offset is half the width
      offset <- width/2*c(areas[lower]/height[lower], .9*min(1, 1/width),
                          rev(areas[upper-1])/rev(height[upper-1]))
    } else {
      if (width.method=="height") {
        height <- 2*c(2^(lower-k), rev(2^(lower-k)))
        offset <- width/2*c(height[lower], 1, rev(height[upper-1]))
      } else stop("parameter width.method is not specified. Use 'linear', 'area', or 'height'")
    }
  }


  if (horizontal) {
    points(x[out], y[out], pch = 1, cex=0.7)

    # bottom rectangles:
    rect(qu[lower], y[lower] + offset[lower],
         qu[lower+1], y[lower] - offset[lower], col = col, border=border)
    # top rectangles:
    rect(qu[upper-1], y[upper] + offset[upper],
         qu[upper], y[upper] - offset[upper], col = col, border=border)

    # draw the median as a line
    med = length(lower) + 1
    oldpar <- par()
    par(lwd=2*oldpar$lwd)
    lines(x=c(qu[med],qu[med]),
          y=c(y[med]-offset[med], y[med]+offset[med]), col=median.col)
    par(lwd=oldpar$lwd)
  } else { # draw vertical plot
    points(y[out], x[out], pch = 1, cex=0.7)
    # bottom rectangles:
    rect(y[lower] + offset[lower], qu[lower],
         y[lower] - offset[lower], qu[lower+1], col = col, border=border)
    # top rectangles:
    rect(y[upper] + offset[upper], qu[upper-1],
         y[upper] - offset[upper], qu[upper], col = col, border=border)

    med = length(lower)+1
    oldpar <- par()
    par(lwd=2*oldpar$lwd)
    lines(x=c(y[med]-offset[med], y[med]+offset[med]),
          y=c(qu[med],qu[med]), col=median.col)
    par(lwd=oldpar$lwd)
  }
}
