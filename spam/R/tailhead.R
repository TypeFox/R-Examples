# This is file ../spam/R/tailhead.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


head.spam <- function(x, n = 6L, m = n, ...) {
  stopifnot(length(n) == 1L, length(m) == 1L)
   n <- if (n < 0L) 
        max(nrow(x) + n, 0L)
    else min(n, nrow(x))
   m <- if (m < 0L) 
        max(ncol(x) + m, 0L)
    else min(m, ncol(x))
    as.matrix(x[seq_len(n), seq_len(m), drop = FALSE])
}

tail.spam <- function (x, n = 6L, m = n, addrownums = TRUE,  ...) 
{
    stopifnot(length(n) == 1L, length(m) == 1L)
    nrx <- nrow(x)
    ncx <- ncol(x)
    n <- if (n < 0L) 
        max(nrx + n, 0L)
    else min(n, nrx)
    m <- if (m < 0L) 
        max(ncx + m, 0L)
    else min(m, ncx)
    selr <- seq.int(to = nrx, length.out = n)
    selc <- seq.int(to = ncx, length.out = n)
    ans <- as.matrix( x[selr, selc, drop = FALSE])
    if (addrownums) {
#    if (addrownums && is.null(rownames(x))) { # rownames is null by default
#        rownames(ans) <- paste0("[", selr, ",]") # can be used from R2.15
#        colnames(ans) <- paste0("[,", selc, "]") 
        rownames(ans) <- paste("[", selr, ",]", sep = "")
        colnames(ans) <- paste("[,", selc, "]", sep = "")
      }
    ans
}


setMethod("head","spam",head.spam)
setMethod("tail","spam",tail.spam)

setMethod("head","spam.chol.NgPeyton", function(x, n = 6L, m = n, ...) head.spam(as.spam(x), n = 6L, m = n, ...))
setMethod("tail","spam.chol.NgPeyton", function(x, n = 6L, m = n, addrownums = TRUE,  ...) tail.spam(as.spam(x), n = 6L, m = n, addrownums = TRUE,  ...))
