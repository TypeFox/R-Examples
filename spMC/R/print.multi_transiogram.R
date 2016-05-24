print.multi_transiogram <-
function(x, ...) {
  for(i in 1:dim(x$Tmat)[3]) {
    cat(x$type, " transition probabilities at lag (", sep = "")
    cat(unlist(x$lags[i, ]), sep = ", ")
    cat(")\n\n")
    print(x$Tmat[, , i], ...)
    if(i != length(x)) cat("\n")
  }
  invisible(x)
}

