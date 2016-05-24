print.transiogram <-
function(x, ...) {
  mylen <- length(x$lags)
  for(i in 1:mylen) {
    cat(x$type, "transition probabilities at lag:", x$lags[i], "\n\n")
    print(x$Tmat[,,i], ...)
    if(i != mylen) cat("\n")
  }
  invisible(x)
}
