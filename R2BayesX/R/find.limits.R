find.limits <- function(map, mar.min = 2, ...)
{
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  n <- length(map)
  myrange <- function(x, c.select = 1L, ...) {
    return(na.omit(x[, c.select], ...))
  }
  xlim <- range(unlist(lapply(map, myrange, c.select = 1L, ...)))
  ylim <- range(unlist(lapply(map, myrange, c.select = 2L, ...)))
  mar <- NULL
  asp <- attr(map, "asp")
  if(is.null(asp))
    asp <- (diff(ylim) / diff(xlim)) / cos((mean(ylim) * pi) / 180)
  if(!is.null(height2width <- attr(map, "height2width"))) {
    height2width <- height2width * 0.8
    if(!is.null(mar.min)) {
      if(height2width > 1) {
        side <- 17.5 * (1 - 1/height2width) + mar.min / height2width
        mar <- c(mar.min, side, mar.min, side)
      }
      else {
        top <- 17.5  * (1 - height2width) + mar.min * height2width
        mar <- c(top, mar.min, top, mar.min)
      }
    }
  }

  return(list(ylim = ylim, xlim = xlim, mar = mar, asp = asp))
}

