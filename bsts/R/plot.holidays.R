PlotHolidays <- function(model, ylim = NULL, same.scale = TRUE, ...) {
  ## Makes a series of side-by-side boxplots for all the holiday state
  ## components in 'model.'
  ##
  ## Args:
  ##   model: A model fit by 'bsts' containing one or more holiday
  ##     state components.
  ##   ylim: limits on the vertical axis of the plots.  If ylim is
  ##     specified, all plots will have the same vertical axis.
  ##   same.scale: If ylim is NULL, this flag determines whether all
  ##     plots share the same scale for there vertical axis
  ##     (same.scale == TRUE), or each plot is independently scaled
  ##     (same.scale == FALSE).
  ##   ...:  Extra arguments passed to boxplot
  ## Returns:
  ##   Invisible NULL.
  number.of.components <- dim(model$state.contributions)[2]

  holiday.list <- list()

  times <- index(model$original.series)
  component.names <- dimnames(model$state.contributions)[[2]]

  for (i in 1:number.of.components) {
    component <- model$state.contributions[, i, ]
    identically.zero <- apply(component, 2, function(x) all(x == 0))
    zero.fraction <- sum(identically.zero) / length(identically.zero)
    if (zero.fraction > .9) {
      holiday <- component[, !identically.zero]
      attr(holiday, "times") <- times[!identically.zero]
      holiday.list[[length(holiday.list) + 1]] <- holiday
      names(holiday.list)[length(holiday.list)] <- component.names[i]
    }
  }

  number.of.holidays <- length(holiday.list)

  if (same.scale && is.null(ylim)) {
    ylim <- range(unlist(holiday.list))
  }

  nr <- floor(sqrt(number.of.holidays))
  nc <- ceiling(number.of.holidays / nr)
  opar <- par(mfrow = c(nr, nc))
  on.exit(par(opar))

  for (i in 1:number.of.holidays) {
    times <- attr(holiday.list[[i]], "times")
    boxplot(holiday.list[[i]],
            names = as.character(attr(holiday.list[[i]], "times")),
            main = names(holiday.list)[i],
            ylim = ylim,
            las = 2,
            ...)
    dt <- diff(times, unit = "days")
    sep = (1:length(times))[dt > 300]
    abline(v = sep + .5, lty = 2, col = "lightgray")
  }
  return(invisible(NULL))
}
