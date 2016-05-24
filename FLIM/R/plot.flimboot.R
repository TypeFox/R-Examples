plot.flimboot <-
function(x, response, grouping=NULL, col=1:20, ylab="Response",
                          xlab="Times", ylim=NULL, main=NULL, ...) {
  stds <- flimSD(x, response, grouping)
  ofo <- x$org
  if(is.null(main)) main <- "Hypothetical mean with bootstrap confidence bands"
  if(is.null(grouping)) {
    mr <- flimMean(ofo, response)[, 1]
    lb <- mr - 2*stds
    ub <- mr + 2*stds
    if(is.null(ylim)) ylim <- c(min(lb), max(ub))
    plot(ofo$times, mr, type="l", col=col[1], ylab=ylab, xlab=xlab, ylim=ylim,
         main=main, ...)
    points(ofo$times, mr, lty=1, col=col[1], ...)
    lines(ofo$times, lb, lty=3, col=col[1])
    lines(ofo$times, ub, lty=3, col=col[1])
  } else {
    mr.mat <- flimMean(ofo, response, grouping)
    nr.levels <- length(unique(ofo$df[, grouping]))
    if(is.null(ylim)) ylim <- c(min(mr.mat[, 1:nr.levels] - 2*stds), 
                                max(mr.mat[, 1:nr.levels] + 2*stds)) 
    plot(ofo$times, flimMean(ofo, response)[, 1], type="n", axes=T, ylab=ylab,
         xlab=xlab, main=main, ylim=ylim, ...)
    for(k in 1:nr.levels){
      mr <- mr.mat[, k]
      std <- stds[, k]
      lb <- mr - 2*std
      ub <- mr + 2*std
      lines(ofo$times, mr, lty=1, col=col[k], ...)
      points(ofo$times, mr, lty=1, col=col[k], ...)
      lines(ofo$times, lb, lty=3, col=col[k])
      lines(ofo$times, ub, lty=3, col=col[k])
    }
  }
}
