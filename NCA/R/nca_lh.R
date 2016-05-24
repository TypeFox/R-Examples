nca_lh <-
function (loop.data, mpy, cutoff, bottleneck.x) {
  x <- loop.data$x
  y <- loop.data$y

  # ceiling line through point with lowest empirical X and highest empirical Y
  # y coordinate for x.low
  y.xlow <- max(y[which(x==loop.data$x.low)])
  # x coordinate for y.high
  x.yhigh <- min(x[which(y==loop.data$y.high)])
  # coordinates for Low and High points
  x.lh <- c(loop.data$x.low,  x.yhigh)
  y.lh <- c(y.xlow, loop.data$y.high)
  
  ceiling   <- 0.5 * (loop.data$y.high - y.xlow) * (x.yhigh - loop.data$x.low)
  effect    <- ceiling / loop.data$scope
  if (effect > 0) {
    slope     <- (loop.data$y.high - y.xlow) / (x.yhigh - loop.data$x.low)
    intercept <- loop.data$y.high - (slope * x.yhigh)
    ineffs    <- p_ineffs(loop.data, intercept, slope)
    above     <- p_above(loop.data, slope, intercept)
  } else {
    slope     <- NA
    intercept <- NA
    ineffs    <- list(x=NA, y=NA, abs=NA, rel=NA)
    above     <- NA    
  }
  
  bottleneck  <- p_bottleneck(loop.data, mpy, slope, intercept, cutoff, bottleneck.x)

  return(list(line=list(x.lh, y.lh),
              ceiling=ceiling, slope=slope, effect=effect,
              intercept=intercept, above=above, ineffs=ineffs,
              bottleneck=bottleneck))
}