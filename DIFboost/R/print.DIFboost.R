print.DIFboost <-
function(x, ...){
    I <- x$I
    dif.mat <- x$dif.mat
    dif.items <- x$DIF.Items
    cat("Number of (valid) persons: P =",x$P,"\n")
    cat("Number of items: I =",I,"\n")
    cat("DIF Items:",dif.items,"\n")
    cat("\n")
    cat("Matrix of estimated item-specific coefficients:\n")
    print.default(dif.mat, ...)
    invisible(x)
  }
