


rvdf <- function (..., row.names = NULL, check.rows = FALSE,
                  check.names = TRUE) {
  L <- list(...)
  L <- lapply(L, as.rvobj)
  dummy <- lapply(L, rvmedian)
  X <- as.data.frame(dummy, row.names=row.names, optional=TRUE)
  X <- structure(as.list(X), row.names=row.names(X))
  for (name in names(L)) {
    X[[name]] <- L[[name]]
  }
  class(X) <- c("rvdf", class(X))
  return(X)  
}

print.rvdf <- function (x, ..., digits=NULL, quote=FALSE, right=TRUE, row.names=TRUE) {
  constants <- sapply(x, is.constant)
  df <- as.data.frame(lapply(x, rvmedian))
  n.cols <- ncol(df)
  n.rows <- nrow(df)
  cat(gettextf("rv data frame with %d columns and %d rows\n", n.cols, n.rows))
  rows.to.print <- 10
  if (n.rows > 2 * rows.to.print) {
    h <- head(df, n=rows.to.print)
    base::print.data.frame(h, ..., digits=digits, quote=quote, right=right, row.names=row.names)
    cat(gettextf("... (%d rows omitted) ...\n", n.cols - 2 * rows.to.print))
    t <- tail(df, n=rows.to.print)
    base::print.data.frame(t, ..., digits=digits, quote=quote, right=right, row.names=row.names)
  } else {
    base::print.data.frame(df, ..., digits=digits, quote=quote, right=right, row.names=row.names)    
  }
    invisible(x)
}

"[.rvdf" <- function (x, i, j, ...) {
  
}
