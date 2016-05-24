PlotSeasonalEffect <- function(bsts.object,
                               nseasons = 7,
                               season.duration = 1,
                               same.scale = TRUE,
                               ylim = NULL,
                               get.season.name = NULL,
                               burn = SuggestBurn(.1, bsts.object),
                               ...) {
  ## Creates a set of plots similar to a 'month plot' showing how the
  ## effect of each season has changed over time.  This function uses
  ## mfrow to create a page of plots, so it cannot be used on the same
  ## page as other plotting functions.
  ##
  ## Args:
  ##   bsts.object:  A bsts model containing a seasonal component.
  ##   ylim:  The limits of the vertical axis.
  ##   same.scale: Used only if ylim is NULL.  If TRUE then all
  ##     figures are plotted on the same scale.  If FALSE then each
  ##     figure is independently scaled.
  ##   nseasons:  The number of seasons in the seasonal component to be plotted.
  ##   season.duration: The duration of each season in the seasonal
  ##     component to be plotted.
  ##   get.season.name: A function taking a Date, POSIXt, or other
  ##     time object used as the index of the original data series
  ##     used to fit 'bsts.object,' and returning a character string
  ##     that can be used as a title for each panel of the plot.  If
  ##     this is NULL and nseasons is one of the following time units
  ##     then the associated function will be used. (see ?weekdays)
  ##     - 4  quarters
  ##     - 7  weekdays
  ##     - 12 months
  ##   burn: The number of MCMC iterations to be discarded as burn-in.
  ##   ...:  Extra arguments passed to PlotDynamicDistribution.
  ##
  ## Returns:
  ##   Invisible NULL.
  ##
  effect.names <- dimnames(bsts.object$state.contributions)$component
  position <- grep("seasonal", effect.names)
  if (length(position) == 1) {
    name.components <- strsplit(effect.names[position], ".", fixed = TRUE)[[1]]
    nseasons <- as.numeric(name.components[2])
    season.duration <- as.numeric(name.components[3])
  } else {
    effect.name <- paste("seasonal", nseasons, season.duration, sep = ".")
    position <- grep(effect.name, effect.names)
  }
  if (length(position) != 1) {
    stop("The desired seasonal effect could not be located. ",
         "Did you specify 'nseasons' and 'season.duration' correctly?")
  }
  effect <- bsts.object$state.contributions[, position, ]
  if (burn > 0) {
    effect <- effect[-(1:burn), , drop = FALSE]
  }
  if (is.null(ylim) && same.scale == TRUE) {
    ylim <- range(effect)
  }
  vary.ylim <- is.null(ylim)
  time <- index(bsts.object$original.series)
  nr <- floor(sqrt(nseasons))
  nc <- ceiling(nseasons / nr)
  if (is.null(get.season.name) && inherits(time, c("Date", "POSIXt"))) {
    if (nseasons == 7) {
      get.season.name <- weekdays
    } else if (nseasons == 12) {
      get.season.name <- months
    } else if (nseasons == 4) {
      get.season.name <- quarters
    }
  }

  opar <- par(mfrow = c(nr, nc))
  on.exit(par(opar))
  for (season in 1:nseasons) {
    time.index <- seq(from = 1 + (season - 1) * season.duration,
                      to = length(time),
                      by = nseasons * season.duration)
    season.effect <- effect[, time.index]
    if (vary.ylim) {
      ylim <- range(season.effect)
    }
    dates <- time[time.index]
    PlotDynamicDistribution(season.effect,
                            dates,
                            ylim = ylim,
                            xlim = range(time),
                            ...)
    lines(dates, apply(season.effect, 2, median), col = "green")
    if (inherits(dates, c("Date", "POSIXt")) && !is.null(get.season.name)) {
      season.name <- get.season.name(dates[1])
      title(main = season.name)
    } else {
      title(main = paste("Season", season))
    }
  }
  return(invisible(NULL))
}
