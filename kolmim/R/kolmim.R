# Original cdf by Wang et al
pkolm <- function (d, n) {
  if (d <= 0 || n <= 0) stop("non-positive argument")
  .Call(K_pKolmogorov2x, as.double(d), as.integer(n), PACKAGE="kolmim")
}

pkolmim <- function (d, n) {
  if (any(d <= 0) || any(n <= 0)) stop("non-positive argument")
  l <- max(length(d), length(n))
  .Call(K_pkolmim, as.double(rep(d, length=l)),
     as.integer(rep(n, length=l)), as.integer(l),
     PACKAGE="kolmim")
}

# Performs one-sample two-sided exact KS test using improved routine
ks.test.imp <- function(x, y, ...) {
  if(is.character(y)) # avoid matching anything in this function
    y <- get(y, mode = "function", envir = parent.frame())
  if(!is.function(y))
    stop("'y' must be a function or a string naming a valid function")

  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if(n < 1L) stop("not enough data")
  if(length(unique(x)) < n)
    warning("ties should not be present for the Kolmogorov-Smirnov test")

  METHOD <- "One-sample two-sided exact Kolmogorov-Smirnov test"
  x <- y(sort(x), ...) - (0 : (n-1)) / n
  STATISTIC <- max(c(x, 1/n - x))
  names(STATISTIC) <- "D"
  PVAL <- 1 - .Call(K_pkolmim2x, as.double(STATISTIC), as.integer(n),
                    PACKAGE = "kolmim")
  PVAL <- min(1.0, max(0.0, PVAL)) # fix possible overshoot
  RVAL <- list(statistic = STATISTIC,
               p.value = PVAL,
               alternative = "two-sided",
               method = METHOD,
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

