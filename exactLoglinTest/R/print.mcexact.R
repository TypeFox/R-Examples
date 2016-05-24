print.bab <- function(x, ...){
  rval <- as.data.frame(rbind(x$dobs, x$phat, x$mcse))
  rownames(rval) <- c("observed.stat", "pvalue", "mcse")
  colnames(rval) <- x$stat(rowlabels = TRUE)
  print(rval)
  invisible(rval)
}

print.cab <- function(x, ...){
  rval <- as.data.frame(rbind(x$dobs, x$phat, x$mcse))
  rownames(rval) <- c("observed.stat", "pvalue", "mcse")
  colnames(rval) <- x$stat(rowlabels = TRUE)
  print(rval)
  invisible(rval)
}

print.babSummary <- function(x, ...){
  cat("Number of iterations       = ", x$startiter - 1, "\n")
  cat("t degrees of freedom       = ", x$tdf, "\n")
  cat("Number of counts           = ", x$n, "\n")
  cat("df                         = ", x$n1, "\n")
  cat("Next update has nosim      = ", x$nosim, "\n")
  cat("Next update has maxiter    = ", x$maxiter, "\n")
  cat("Proportion of valid tables = ", x$perpos, "\n")
  cat("\n")

  rval <- as.data.frame(rbind(x$dobs, x$phat, x$mcse))
  rownames(rval) <- c("observed.stat", "pvalue", "mcse")
  colnames(rval) <- x$stat(rowlabels = TRUE)
  print(rval)
  invisible(rval)
}

print.cabSummary <- function(x, ...){
  cat("Number of iterations       = ", x$startiter - 1, "\n")
  cat("t degrees of freedom       = ", x$tdf, "\n")
  cat("Number of counts           = ", x$n, "\n")
  cat("df                         = ", x$n1, "\n")
  cat("Number of batches          = ", x$nobatches, "\n")
  cat("Batchsize                  = ", x$batchsize, "\n")
  cat("Next update has nosim      = ", x$nosim, "\n")
  cat("Proportion of valid tables = ", x$perpos, "\n")
  cat("\n")
  
  rval <- as.data.frame(rbind(x$dobs, x$phat, x$mcse))
  rownames(rval) <- c("observed.stat", "pvalue", "mcse")
  colnames(rval) <- x$stat(rowlabels = TRUE)
  print(rval)
  invisible(rval)
}
