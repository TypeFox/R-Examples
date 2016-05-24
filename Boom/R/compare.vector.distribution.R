# Copyright 2012 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

CompareVectorBoxplots <- function(draws, main = NULL, colors = NULL,
                                  burn = 0) {
  ## Creates a boxplot comparing the distributions of several vector
  ## valued parameters.
  ## Args:
  ##   draws: A list of MCMC draws.  Each list element is a matrix
  ##     with rows corresponding to MCMC iterations and columns to
  ##     variables.
  ##   main:  Main title of the plot.
  ##   colors:  Colors to use for the boxplots.  This should either be
  ##   burn: The number of initial MCMC iterations to discard before
  ##     making the plot.
  stopifnot(is.list(draws))
  stopifnot(all(sapply(draws, is.matrix)))
  dimension <- ncol(draws[[1]])
  stopifnot(all(sapply(draws, ncol) == dimension))
  number.of.components <- length(draws)

  list.position <- 0
  y <- list()
  for (column in 1:dimension) {
    for (component in 1:number.of.components) {
      list.position <- list.position + 1
      y[[list.position]] <- draws[[component]][, column]
      if (burn > 0) y[[list.position]] <- y[[list.position]][-(1:burn)]
    }
  }

  if (length(colors) == 1) {
    colors <- rep(colors, number.of.components)
  }
  if(!is.null(colors)) {
    stopifnot(length(colors) == number.of.components)
  }
  boxplot(y, col = colors, main = main, axes = FALSE)
  axis(2)
  box()
  split.positions <- number.of.components * (1:(dimension-1)) + .5
  tick.positions <- number.of.components * (0:(dimension-1)) +
    .5*(number.of.components + 1)

  vnames <- colnames(draws[[1]])
  if (is.null(vnames)) {
    tick.labels <- paste(1:dimension)
  } else {
    tick.labels <- vnames
  }
  axis(side = 1, at = tick.positions, labels = tick.labels)
  abline(v = split.positions)

  component.names <- names(draws)
  if (is.null(component.names)) {
    component.names <- paste(0:(number.of.components - 1))
  }

  if(!is.null(colors)) {
    legend("topright", fill = colors,
           legend = component.names,
           title = "Component",
           bg = "white")
  }

  return(invisible(NULL))
}
