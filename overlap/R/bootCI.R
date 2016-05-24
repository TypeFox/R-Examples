# Calculate confidence intervals based on a set of bootstrap estimates.
# In bootCIlogit, the corrections are perfomed on the logistic scale.

# t0 is point estimate of parameter of interest
# bt is a vector of bootstrap estimates
# conf is required confidence interval
# Returns: matrix with estimators in rows, lower/upper limits in columns.

bootCI <-
function(t0, bt, conf=0.95)  {
  bt <- as.vector(bt)
  bt <- bt[is.finite(bt)]
  out <- matrix(NA, 5, 2)
  dimnames(out) <- list(c("norm", "norm0", "basic", "basic0", "perc"),
                        c("lower","upper"))
  bias <- mean(bt) - t0
  merr <- sd(bt) * qnorm((1 + conf)/2)
  out[1, ] <- c(t0 - bias - merr, t0 - bias + merr) # norm, same as boot.ci
  out[2, ] <- c(t0 - merr, t0 + merr)               # norm0
  out[5, ] <- quantileInter(bt, conf)               # perc
  out[3, ] <- 2 * t0 - out[5, 2:1]                  # basic, same as boot.ci
  out[4, ] <- out[5, ] - bias                       # basic0
  return(out)
}

bootCIlogit <-
function(t0, bt, conf=0.95)  
  plogis(bootCI(qlogis(t0), qlogis(bt), conf=conf))

# Normally-interpolated quantile confidence interval.
quantileInter <-
function(bt, conf=0.95)
# Args:
#   bt : numeric vector, no missing values
#   conf : required confidence interval
# Returns: a vector of length 2, lower and upper confidence limits,
#   or NAs if the vector is not long enough.
# Not exported.
{
  R <- length(bt)
  alpha <- (1 + c(-conf, conf)) / 2
  rk <- (R+1) * alpha
  if (!all(rk>1 & rk<R) )  {
    out <- rep(NA_real_, 2)
  } else {
    k <- trunc(rk)
    ts <- sort(bt, partial = sort(c(k, k+1)))[c(k, k+1)]
    if(all(k == rk)) {
      out <- ts[1:2]
    } else {
      wanted <- qnorm(alpha)
      tooLo  <- qnorm(k/(R+1))
      tooHi  <- qnorm((k+1)/(R+1))
      tLo <- ts[1:2]
      tHi <- ts[3:4]
      out <- tLo + (wanted-tooLo)/(tooHi - tooLo)*(tHi - tLo)  # interpolation
    }
  }
  return(out)
}
