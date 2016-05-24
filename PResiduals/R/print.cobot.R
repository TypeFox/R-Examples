#' cobot class print method
#' @param x cobot object
#' @param ... arguments passed to print.default
#' @keywords print
#' @export
#' @method print cobot

print.cobot <- function(x, ...) {
  y <- matrix(nrow=length(x$TS),ncol=3)
  dims <- character(length(x$TS))
  for (i in 1:length(x$TS)){
    y[i,] <- c(x$TS[[i]]$ts, sqrt(x$TS[[i]]$var),x$TS[[i]]$pval)
    dims[i] <- x$TS[[i]]$label
  }
  dimnames(y) <- list(dims,c('est','stderr','p'))
  invisible(print(y,...))
  cat('Fisher Transform:',x$fisher,'\n')
  cat(format(x$conf.int*100,digits=3),'% Confidence Interval: (',
      x$TS[[2]]$lower,', ',x$TS[[2]]$upper,')\n',sep='')
  cat('Number of Observations:',x$data.points,'\n')
  invisible(y)
}
