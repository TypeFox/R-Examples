##' show both head and tail in summary with ellipsis in between
##'
##' @param x object to summarize
##' @param k number of elements taken from head and tail
##' @return calls cat on value
##' @export
headtail <- function(x, k=3) {
  out <- capture.output(x); n <- length(out)
  val <- paste(c(out[1:(k+1)], "   ...", out[(n-k):n]), collapse="\n")
  cat(val, "\n")
  invisible()
}


## headtail.matrix <- function(x, k=3) {
##   headtail.data.frame(x, k)
## }



## confint for htest object
##' confint for htest objects
##'
##' @param object an htest object
##' @param parm ignored
##' @param level ignored
##' @return a confidence interval
##'
##' @export
confint.htest <- function(object, parm, level, ...) {
  x <- object
  tmp <- x$conf.int
  transform <- list(...)$transform
  if(!is.null(transform))
    tmp <- transform(tmp)
  cat(sprintf("(%.2f, %.2f) with %d percent confidence",
              tmp[1], tmp[2], 
              floor(100*attr(tmp, "conf.level"))))
  invisible()
}
##
