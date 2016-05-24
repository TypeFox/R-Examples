.makehnp <-
function(obj, conf, halfnormal, how.many.out, paint.out, col.paint.out, print.on, plot.sim, ...) {
  # A few checks
  if(print.on) how.many.out <- T
  if(paint.out) {
    how.many.out <- T
    if(missing(col.paint.out)) col.paint.out <- 2
  }
  # Residuals and simulated envelope
  res.original <- obj[,1]
  res <- obj[,-1]
  env <- apply(res, 1, quantile, c((1-conf)/2, .5, (1-conf)/2+conf))
  # Saving / plotting
  n <- nrow(res)
  i <- 1:n
  if(halfnormal) q.x <- qnorm((i+n-1/8)/(2*n+1/2)) else q.x <- qnorm((i-3/8)/(n+1/4))
  simdata <- list(q.x, t(env)[,1], t(env)[,2], t(env)[,3], res.original)
  class(simdata) <- "hnp"
  names(simdata) <- c("x", "lower", "median", "upper", "residuals")
  if(how.many.out) {
    mat <- cbind(t(env), res.original, q.x)
    out <- sum(mat[,4] > mat[,3] | mat[,4] < mat[,1])
    if(paint.out) {
      simdata$out.index <- matrix(mat[mat[,4] > mat[,3] | mat[,4] < mat[,1], 4:5], ncol=2)
      simdata$col.paint.out <- col.paint.out
    }
    simdata$how.many.out <- TRUE
    simdata$total <- nrow(mat)
    simdata$out <- out
    simdata$print.on <- print.on
    simdata$paint.out <- paint.out
  } else {
    simdata$how.many.out <- FALSE
  }
  if(plot.sim) {
    plot(simdata, ...)
    return(invisible(simdata))
  } else {
    return(simdata)
  }
}
