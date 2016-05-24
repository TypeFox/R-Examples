#' cocobot class print method
#' @method print cocobot
#' @param x cocobot object
#' @param ... arguments passed to print.default
#' @keywords print
#' @export
print.cocobot <- function(x, ...) {
  y <- matrix(nrow=length(x$TS),ncol=5)
  dims <- character(length(x$TS))
  for (i in 1:length(x$TS)){
    y[i,] <- c(x$TS[[i]]$ts, sqrt(x$TS[[i]]$var),x$TS[[i]]$pval,x$TS[[i]]$lower,x$TS[[i]]$upper)
    dims[i] <- x$TS[[i]]$label
  }
  dimnames(y) <- list(dims,c('est','stderr','p','lower CI','upper CI'))
  invisible(print(y,...))
  cat('Fisher Transform:',x$fisher,'\n')
  cat('Confidence Interval: ', format(x$conf.int*100,digits=3),'%\n', sep='')
  cat('Number of Observations:',x$data.points,'\n')
  invisible(y)
}
