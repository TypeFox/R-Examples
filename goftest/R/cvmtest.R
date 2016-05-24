##
## cvmtest.R
##
## Cramer-von Mises test
##
## $Revision: 1.4 $ $Date: 2014/06/24 02:13:27 $
##

cvm.test <- function(x, null="punif", ..., nullname) {
  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if(is.character(null)) nulltext <- null
  if(missing(nullname) || is.null(nullname)) {
    reco <- recogniseCdf(nulltext)
    nullname <- if(!is.null(reco)) reco else 
                paste("distribution", sQuote(nulltext))
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  n <- length(x)
  F0 <- if(is.function(null)) null else
        if(is.character(null)) get(null, mode="function") else
        stop("Argument 'null' should be a function, or the name of a function")
  U <- F0(x, ...)
  if(any(U < 0 | U > 1))
    stop("null distribution function returned values outside [0,1]")
  U <- sort(U)
  k <- seq_len(n)
  omega2 <- 1/(12 * n) + sum((U - (2*k - 1)/(2*n))^2)
  PVAL <- pCvM(omega2, n=n, lower.tail=FALSE)
  names(omega2) <- "omega2"
  METHOD <- c("Cramer-von Mises test of goodness-of-fit",
              paste("Null hypothesis:", nullname))
  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if(length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for(i in seq_along(parnames))
      pard[i] <- paste(parnames[i], "=", paste(pars[[i]], collapse=" "))
    pard <- paste("with",
                  ngettext(length(pard), "parameter", "parameters"),
                  "  ", 
                  paste(pard, collapse=", "))
    METHOD <- c(METHOD, pard)
  }
  out <- list(statistic = omega2,
               p.value = PVAL,
               method = METHOD,
               data.name = xname)
  class(out) <- "htest"
  return(out)
}

