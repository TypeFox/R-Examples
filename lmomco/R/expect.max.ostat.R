"expect.max.ostat" <-
function(n, para=NULL, cdf=NULL, pdf=NULL, qua=NULL,
         j=NULL, lower=-Inf, upper=Inf, aslist=FALSE, ...) {

   if(is.null(j)) j <- n
   if(j > n)         stop("j can not be greater than n")
   if(is.null(para)) stop("parameter list not specified and it needs to be in lmomco style")

   if(is.null(qua)) {
      if(is.null(cdf)) stop("cdf function using lmomco parameter style not specified")
      if(is.null(pdf)) stop("pdf function using lmomco parameter style not specified")
   
      "fna" <- function(x) {
         F <- cdf(x, para=para, ...)
         #a <- F^(j-1); b <- (1-F)^(n-j)
         # Say F = 0, n=10, and j=1, then a --> 0*log(0) is NaN
         # but b --> 10*log(0) is -Inf, -Inf is not ok, but more likely
         # to still continue propogation of computations, but rounding
         # will likely still compound, override and set to zero
         a <- (j-1)*log(F); b <- (n-j)*log(1-F)
         a[is.nan(a)]      <- 0;      b[is.nan(b)] <- 0
         a[! is.finite(a)] <- 0; b[! is.finite(b)] <- 0
         p <- pdf(x, para=para, ...)
         p[is.nan(p)] <- 0; p[! is.finite(p)] <- 0
         v <- p*exp(log(x*x) + a + b)/x # the x*x, /x to for negative x
         return(v)
      }
      tmp <- NULL
      try(tmp <- integrate(fna, lower, upper, subdivisions = 200L))
      if(is.null(tmp)) return(NA)
      v <- tmp$value
      VAL <- exp(log(v*v) - lbeta(j,n-j+1))/v # the v*v, /v for negative v
      if(aslist) {
         zz <- list(type="bypdfcdf", value=VAL,
                    abs.error=tmp$abs.error, subdivisions=tmp$subdivisions, message=tmp$message)
         return(zz)
      } else {
         return(VAL)
      }
   } else {

      if(is.null(qua)) stop("qua function using lmomco parameter style not specified")
      if(! is.finite(lower)) lower <- 0
      if(! is.finite(upper)) upper <- 1   
      "fnb" <- function(F) {
         x <- qua(F, para=para, ...)
         a <- (j-1)*log(F); b <- (n-j)*log(1-F) # in log
         a[is.nan(a)] <- 0; b[is.nan(b)] <- 0
         a[! is.finite(a)] <- 0; b[! is.finite(b)] <- 0
         return(exp(log(x*x) + a + b)/x) # the x*x, /x to for negative x
      }
      tmp <- NULL
      try(tmp <- integrate(fnb, lower, upper, subdivisions = 200L))
      if(is.null(tmp)) return(NA)
      v <- tmp$value
      N <- log(n) + lchoose(n-1, j-1) # in log
      VAL <- exp(log(v*v) + N)/v # the v*v, /v for negative v
      if(aslist) {
         zz <- list(type="byqua", value=VAL,
                    abs.error=tmp$abs.error, subdivisions=tmp$subdivisions, message=tmp$message)
         return(zz)
      } else {
         return(VAL)
      }
   }
}

"expect.min.ostat" <-
function(n, ...) {
  return(expect.max.ostat(n, j=1, ...))
}


"eostat" <-
function(...) {
  return(expect.max.ostat(...))
}

