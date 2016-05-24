plot.snpRF <- function(x, type="l", main=deparse(substitute(x)), ...) {
  if(x$type == "unsupervised")
    stop("No plot for unsupervised snpRF.")
  test <- !is.null(x$test$err.rate)

  err <- x$err.rate
  if(test) err <- cbind(err, x$test$err.rate)

  if(test) {
    colnames(err) <- c("OOB", "Test")
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  } else {
    matplot(1:x$ntree, err, type = type, xlab="trees", ylab="Error",
            main=main, ...)
  }
  invisible(err)
}

  
