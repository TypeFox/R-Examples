plot.learnPattern <- function(x, type="l", main=deparse(substitute(x)), ...) {
  if(is.null(x$ooberrors))
    stop("No out-of-bag (OOB) error returned by learnPattern")

  matplot(1:x$ntree, x$ooberrors, type = type, xlab="trees", ylab="OOB Error",
            main=main, ...)
            
  invisible(x$ooberrors)
}

  
