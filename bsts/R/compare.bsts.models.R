CompareBstsModels <- function(model.list,
                              burn = SuggestBurn(.1, model.list[[1]]),
                              filename = "",
                              colors = NULL,
                              lwd = 2,
                              xlab = "Time",
                              main = "",
                              grid = TRUE) {
  ## Produce a set of line plots showing the cumulative absolute one
  ## step ahead prediction errors for different models.  This plot not
  ## only shows which model is doing the best job predicting the data,
  ## it highlights regions of the data where the predictions are
  ## particularly good or bad.
  ##
  ## Args:
  ##   model.list:  A list of bsts models.
  ##   burn: The number of initial MCMC iterations to remove from each
  ##     model as burn-in.
  ##   filename: A string.  If non-empty string then a pdf of the
  ##     plot will be saved in the specified file.
  ##   colors: A vector of colors to use for the different lines in
  ##     the plot.  If NULL then the rainbow pallette will be used.
  ##   lwd: The width of the lines to be drawn.
  ##   xlab: Labels for the horizontal axis.
  ##   main: Main title for the plot.
  ##   grid: Logical.  Should gridlines be drawn in the background?
  ##
  ## Returns:
  ##   Nothing (invisible NULL).
  stopifnot(is.list(model.list))
  stopifnot(length(model.list) > 1)
  stopifnot(all(sapply(model.list, inherits, "bsts")))
  time.dimension <-
    sapply(model.list, function(m) {dim(m$state.contributions)[3]})
  stopifnot(all(time.dimension == time.dimension[1]))

  model.names <- names(model.list)
  if (is.null(model.names)) {
    model.names <- paste("Model", 1:length(model.list))
  }
  number.of.models <- length(model.list)
  if (filename != "") pdf(filename)
  opar <- par(mfrow=c(2, 1))
  original.margins <- c(5.1, 4.1, 4.1, 2.1)
  margins <- original.margins
  opar$mar <- original.margins
  margins[1] <- 0
  par(mar = margins)
  cumulative.errors <-
    matrix(nrow = number.of.models,
           ncol = dim(model.list[[1]]$state.contributions)[3])
  for (i in 1:number.of.models) {
    prediction.errors <-
        bsts.prediction.errors(model.list[[i]])[-(1:burn), , drop = FALSE]
    cumulative.errors[i, ] <- cumsum(abs(colMeans(prediction.errors)))
  }

  if (is.null(colors)) colors <- c("black", rainbow(number.of.models-1))

  idx <- index(as.zoo(model.list[[1]]$original.series))
  plot(zoo(cumulative.errors[1, ], order.by = idx),
       ylim = range(cumulative.errors),
       ylab = "cumulative absolute error",
       xaxt = "n",
       lwd = lwd,
       col = colors[1],
       yaxs = "i",
       main = main)
  axis(2)
  for (i in 2:number.of.models) {
    lines(zoo(cumulative.errors[i, ], order.by = idx),
          lty = i,
          col = colors[i],
          lwd = lwd)
  }

  if (grid) {
    grid()
  }

  legend("topleft",
         model.names,
         lty = 1:number.of.models,
         col = colors,
         bg  = "white",
         lwd = lwd)

  margins <- original.margins
  margins[3] <- 0
  par(mar = margins)
  plot(model.list[[1]]$original.series,
       main = "",
       ylab = "scaled values",
       xlab = xlab,
       yaxs = "i")
  if (grid) {
    grid()
  }
  par(opar)
  if (filename != "") dev.off()
  return(invisible(NULL))
}
