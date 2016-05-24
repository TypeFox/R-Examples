## from <http://mathworld.wolfram.com/LambertW-Function.html>,
## ultimately from Corless et al 1996
lWasymp <- function(z,logz) {
  if (!missing(z) && any(!is.finite(log(z))))
    stop("overflow error: please specify logz")
  L1 <- if (!missing(logz)) logz else log(z)
  L2 <- log(logz)
  L1-L2+L2/L1+L2*(L2-2)/(2*L1^2)+L2*(6-9*L2+2*L2^2)/(6*L1^3)+
    L2*(-12+36*L2-22*L2^2+3*L2^3)/(12*L1^4)+L2*(60-300*L2+350*L2^2-125*L2^3+12*L2^4)/(60*L1^5)
}


lambertW <- function(z,...) {
  bigz <- Mod(z)>1e307
  if (!any(na.omit(bigz))) return(lambertW_base(z,...))
  res <- rep(NA_real_,length(z))
  res[bigz] <- lWasymp(logz=log(z[bigz]))
  res[!bigz] <- lambertW_base(z[!bigz],...)
  res
}

lambertW_base <- function(z,b=0,maxiter=10,eps=.Machine$double.eps,min.imag=1e-9) {
  badz <- !is.finite(z)
  z.old <- z
  z <- z[!badz]
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z <- as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1))
  c = (c > 1.45 - 1.1*abs(b))
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    ## if (any(is.na(t) | is.na(w))) {
    ## FIXME: what to do here?
    ## stop("NAs encountered in LambertW")
    ## print(t)
    ## print(w)
    ## }
    ok <- is.finite(t) & is.finite(w)
    ## iterate until ALL values converged -- perhaps
    ## inefficient
    if (all(abs(Re(t[ok])) < (2.48*eps)*(1.0 + abs(Re(w[ok])))
            & abs(Im(t[ok])) < (2.48*eps)*(1.0 + abs(Im(w[ok])))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
        ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w[!is.na(w)])<min.imag)) w = as.numeric(w)
  if (sum(badz)>0) {
    w.new <- rep(NA_real_,length(z.old))
    w.new[!badz] <- w
    w <- w.new
  }
  return(w)
}
