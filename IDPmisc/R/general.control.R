## general.control.R

## private functions for plot.rose
general.control <- function(stacked = FALSE,
                            rose.rad = NULL,
                            rose.x = NULL,
                            rose.y = NULL,
                            mar = rep(0.3,4),
                            rev.col = FALSE,
                            shift = 0,
                            cex = 1,
                            col = NULL,
                            lty = 1:3,
                            lwd = 1,
                            type = "s")
  ## Author: Rene Locher
  ## Version: 2009-03-16
  ## helper function for plot.rose

  {
    if (is.null(rose.rad)) auto <- TRUE else auto <- FALSE

    if (!is.null(rose.rad) && !is.unit(rose.rad))
      rose.rad <- unit(rose.rad,"mm")

    if (!auto && (!is.null(rose.x) || !is.null(rose.y)))
      warning("When 'rad.rose' is NULL 'rad.x' and 'rad.y' are ignored!")

    if (!is.null(rose.x) && !is.unit(rose.x))
      rose.x <- unit(rose.x,"mm")

    if (!is.null(rose.y) && !is.unit(rose.y))
      rose.y <- unit(rose.y,"mm")

    if (!is.element(type,c("s","l"))) {
        type <- "s"
        warning(paste("type",type,"is not a valid option. Type 's' is used instead"))}

    return(list(rose =
                list(rad = rose.rad,
                     x = rose.x,
                     y = rose.y,
                     auto = auto),
                mar = mar,
                stacked = stacked,
                rev.col = rev.col,
                shift = shift,
                cex = cex,
                col = col,
                lty = lty,
                lwd = lwd,
                ncp = 1000,
                type = type))
  } ## general.control
