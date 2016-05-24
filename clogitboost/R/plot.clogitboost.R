
plot.clogitboost<-function(x, d, grid = NULL, ...){
  fit <- x
  gridstep <- (fit$xmax[d] - fit$xmin[d]) / 100
  if (is.null(grid)) {
    grid <- seq(fit$xmin[d], fit$xmax[d], by = gridstep)
  }
  plot(grid, marginal(fit, grid, d), ...)
}
