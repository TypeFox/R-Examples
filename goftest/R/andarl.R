##
## andarl.R
##
##  Anderson-Darling test and null distribution
##
## $Revision: 1.6 $ $Date: 2014/06/24 02:12:20 $
##

ad.test <- function(x, null="punif", ..., nullname) {
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
  ## call Marsaglia C code
  z <- .C("ADtestR",
          x = as.double(U),
          n = as.integer(n),
          adstat = as.double(numeric(1)),
          pvalue = as.double(numeric(1))
          )
  STATISTIC <- z$adstat
  names(STATISTIC) <- "An"
  PVAL <- z$pvalue
  METHOD <- c("Anderson-Darling test of goodness-of-fit",
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
  out <- list(statistic = STATISTIC,
               p.value = PVAL,
               method = METHOD,
               data.name = xname)
  class(out) <- "htest"
  return(out)
}

pAD <- function(q, n=Inf, lower.tail=TRUE, fast=TRUE) {
  q <- as.numeric(q)
  p <- rep(NA_real_, length(q))
  if(any(ones <- is.infinite(q) & (q == Inf)))
    p[ones] <- 1
  if(any(zeroes <- (is.finite(q) & q <= 0) | (is.infinite(q) & (q == -Inf))))
    p[zeroes] <- 0
  ok <- is.finite(q) & (q > 0)
  nok <- sum(ok)
  if(nok > 0) {
    if(is.finite(n)) {
      z <- .C("ADprobN",
              a       = as.double(q[ok]),
              na      = as.integer(nok),
              nsample = as.integer(n),
              prob    = as.double(numeric(nok))
              )
      p[ok] <- z$prob
    } else if(fast) {
      ## fast version adinf()
      z <- .C("ADprobApproxInf",
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok))
              )
      p[ok] <- z$prob
    } else {
      ## slow, accurate version ADinf()
      z <- .C("ADprobExactInf",
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok))
              )
      p[ok] <- z$prob
    }
      
  }
  if(!lower.tail)
    p <- 1 - p
  return(p)
}

qAD <- local({

  f <- function(x, N, P, Fast) {
    pAD(x, N, fast=Fast) - P
  }
    
  qAD <- function(p, n=Inf, lower.tail=TRUE, fast=TRUE) {
    ## quantiles of null distribution of Anderson-Darling test statistic
    stopifnot(all(p >= 0))
    stopifnot(all(p <= 1))
    if(!lower.tail) p <- 1-p
    ans <- rep(NA_real_, length(p))
    for(i in which(p >= 0 & p < 1)) 
      ans[i] <- uniroot(f, c(0, 1), N=n, P=p[i], Fast=fast, extendInt="up")$root
    return(ans)
  }

  qAD
})


  

