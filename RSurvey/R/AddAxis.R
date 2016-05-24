# Adds an axis to the current plot.

AddAxis <- function(side, lim, ticks.inside=FALSE, minor.ticks=FALSE, ...) {

  ## Additional functions

  # Locate minor tickmark locations

  LocateMinorTicks <- function() {
    at.major <- axTicks(side)
    x1 <- min(at.major)
    x2 <- max(at.major)
    n <- length(at.major) - 1

    delta <- (x2 - x1) / n
    x1 <- x1 - delta
    x2 <- x2 + delta
    n <- n + 2

    if (inherits(lim, "POSIXt")) {
      len.minor <- 6 * n + 1
      at.minor <- pretty(lim, n=len.minor)
    } else {
      typ <- round(delta / 10^floor(log10(delta) + .Machine$double.eps))
      mult <- if (typ == 1 || typ == 2) 5 else 10

      i <- 0
      no.match <- TRUE
      while (no.match) {
        i <- i + 1
        len.minor <- mult * n + i
        at.minor <- seq(x1, x2, length=len.minor)

        if (i > 100) {
          msg <- "Unexpected interval; minor tickmarks not plotted."
          warning(call.=FALSE, msg)
          return(NULL)
        } else if (all(at.major %in% at.minor)) {
          no.match <- FALSE
        }
      }
    }
    at.minor[!at.minor %in% at.major]
  }

  ## Main program

  add.labels <- if (side %in% c(1, 2)) TRUE else FALSE
  las <- if (side %in% c(1, 3)) 1 else 0
  tcl.dir <- if (ticks.inside) 1 else -1

  tcl.major <- tcl.dir * (0.50 / (6 * par("csi")))
  tcl.minor <- tcl.dir * (0.25 / (6 * par("csi")))

  lwd.ticks <- 0.5 * (96 / (6 * 12))

  if (inherits(lim, c("integer", "numeric"))) {
    axis(side, tcl=tcl.major, cex.axis=0.8, las=las,
         labels=add.labels, lwd=-1, lwd.ticks=lwd.ticks, ...)
    if (minor.ticks) {
      at <- LocateMinorTicks()
      if (!is.null(at))
        axis(side, at=at, tcl=tcl.minor, labels=FALSE,
             lwd=-1, lwd.ticks=lwd.ticks)
    }
  } else if (inherits(lim, "POSIXt")) {
    axis.POSIXct(side, lim, tcl=tcl.major, cex.axis=0.8, las=las,
                 labels=add.labels, lwd=-1, lwd.ticks=lwd.ticks, ...)
    if (minor.ticks) {
      at <- LocateMinorTicks()
      axis.POSIXct(side, at=at, tcl=tcl.minor, labels=FALSE,
                   lwd=-1, lwd.ticks=lwd.ticks)
    }
  }
}
