# Find parameters corresponding to mean and sd, given alpha
setMomentsFMstable <- function(mean=1, sd=1,
    alpha, oneminusalpha, twominusalpha){
  if(!missing(alpha)){
    # Case where alpha is specified
    if(length(alpha) != 1) stop(
      "setMomentsFMstable: alpha must be of length 1")
    if(!missing(oneminusalpha) || !missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    oneminusalpha <- 1 - alpha
    twominusalpha <- 2 - alpha
  } else if(!missing(oneminusalpha)){
    # Case where oneminusalpha is specified
    if(length(oneminusalpha) != 1) stop(
      "setMomentsFMstable: oneminusalpha must be of length 1")
    if(!missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    alpha <- 1 - oneminusalpha
    twominusalpha <- 1 + oneminusalpha
  } else {
    if(missing(twominusalpha))stop(
      "Specify one of alpha, oneminusalpha, twominusalpha")
    # Case where twominusalpha is specified
    if(length(twominusalpha) != 1) stop(
      "setMomentsFMstable: twominusalpha must be of length 1")
    alpha <- 2 - twominusalpha
    oneminusalpha <- twominusalpha - 1
  }
  if(length(mean) != 1) stop("setMomentsFMstable: mean must be of length 1")
  if(length(sd) != 1) stop("setMomentsFMstable: sd must be of length 1")
  if(alpha <= 0 || twominusalpha < 0) stop(
    "setMomentsFMstable: alpha must be >= 0 and <= 2")

  cv <- sd/mean
  logmr <- log1p(cv*cv)

  # Case when alpha is less than 0.5 
  if(alpha < 0.5){
    scale.to.alpha <- logmr/(2- 2^alpha)
    logscale <- log(scale.to.alpha)/alpha
    location <- -log(mean) - scale.to.alpha
  } else if(oneminusalpha == 0){
    # Case when alpha is precisely 1
    scale <- pi*logmr/(4*log(2))
    logscale <- log(scale)
    location <- logscale*2*scale/pi - log(mean)
  } else{
    # Case when alpha exceeds 0.5 but is not precisely 1
    s <- if(alpha < 1.5) sin(oneminusalpha*.5*pi) else -cos(twominusalpha*.5*pi)
    sinalpha <- if(alpha < 1) sin(alpha*.5*pi) else sin(twominusalpha*.5*pi)
    scale.to.alpha <- -.5*logmr*s/expm1(-oneminusalpha*log(2))
    logscale <- log(scale.to.alpha)/alpha
    scale <- exp(logscale)
    location <- (scale*sinalpha - scale.to.alpha)/s - log(mean)
  }

  res <- list(alpha=alpha, oneminusalpha=oneminusalpha,
    twominusalpha=twominusalpha, location=location, logscale=logscale,
    created.by="setMomentsFMstable")
  class(res) <- "stableParameters"
  return(res)
}
# Find parameters corresponding to mean and sd, given alpha
setMomentsFMstable <- function(mean=1, sd=1,
    alpha, oneminusalpha, twominusalpha){
  if(!missing(alpha)){
    # Case where alpha is specified
    if(length(alpha) != 1) stop(
      "setMomentsFMstable: alpha must be of length 1")
    if(!missing(oneminusalpha) || !missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    oneminusalpha <- 1 - alpha
    twominusalpha <- 2 - alpha
  } else if(!missing(oneminusalpha)){
    # Case where oneminusalpha is specified
    if(length(oneminusalpha) != 1) stop(
      "setMomentsFMstable: oneminusalpha must be of length 1")
    if(!missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    alpha <- 1 - oneminusalpha
    twominusalpha <- 1 + oneminusalpha
  } else {
    if(missing(twominusalpha))stop(
      "Specify one of alpha, oneminusalpha, twominusalpha")
    # Case where twominusalpha is specified
    if(length(twominusalpha) != 1) stop(
      "setMomentsFMstable: twominusalpha must be of length 1")
    alpha <- 2 - twominusalpha
    oneminusalpha <- twominusalpha - 1
  }
  if(length(mean) != 1) stop("setMomentsFMstable: mean must be of length 1")
  if(length(sd) != 1) stop("setMomentsFMstable: sd must be of length 1")
  if(alpha <= 0 || twominusalpha < 0) stop(
    "setMomentsFMstable: alpha must be >= 0 and <= 2")

  cv <- sd/mean
  logmr <- log1p(cv*cv)

  # Case when alpha is less than 0.5 
  if(alpha < 0.5){
    scale.to.alpha <- logmr/(2- 2^alpha)
    logscale <- log(scale.to.alpha)/alpha
    location <- -log(mean) - scale.to.alpha
  } else if(oneminusalpha == 0){
    # Case when alpha is precisely 1
    scale <- pi*logmr/(4*log(2))
    logscale <- log(scale)
    location <- logscale*2*scale/pi - log(mean)
  } else{
    # Case when alpha exceeds 0.5 but is not precisely 1
    s <- if(alpha < 1.5) sin(oneminusalpha*.5*pi) else -cos(twominusalpha*.5*pi)
    sinalpha <- if(alpha < 1) sin(alpha*.5*pi) else sin(twominusalpha*.5*pi)
    scale.to.alpha <- -.5*logmr*s/expm1(-oneminusalpha*log(2))
    logscale <- log(scale.to.alpha)/alpha
    scale <- exp(logscale)
    location <- (scale*sinalpha - scale.to.alpha)/s - log(mean)
  }

  res <- list(alpha=alpha, oneminusalpha=oneminusalpha,
    twominusalpha=twominusalpha, location=location, logscale=logscale,
    created.by="setMomentsFMstable")
  class(res) <- "stableParameters"
  return(res)
}
if(FALSE){
yearDsn <- setMomentsFMstable(mean=10, sd=5, alpha=1.7)
# Compute mean by integration
f <- function(x) dFMstable(x, yearDsn) * x
integrate(f, lower=-Inf, upper=Inf, rel.tol=1.e-12)$value -10
f <- function(x) dFMstable(x, yearDsn) * x*x
integrate(f, lower=-Inf, upper=Inf, rel.tol=1.e-6)$value -125

yearDsn <- setMomentsFMstable(mean=10, sd=5, alpha=1.000001)
yearDsn
yearDsn <- setMomentsFMstable(mean=10, sd=5, alpha=.999999)
yearDsn
yearDsn <- setMomentsFMstable(mean=10, sd=5, alpha=1)
yearDsn
}

print.stableParameters <- function(x, ...){
  cat(paste(" ******** Stable parameter object created by", x$created.by))
  cat(paste("\n ** Alpha =", x$alpha,"  location =", x$location,
    "  logscale =", x$logscale))
  if(abs(x$oneminusalpha) < .1) cat(paste("\n ** oneminusalpha =",
    x$oneminusalpha))
  if(abs(x$twominusalpha) < .1) cat(paste("\n ** twominusalpha =",
    x$twominusalpha))
  if(max(abs(1:2 -x$alpha - c(x$oneminusalpha, x$twominusalpha))) >
    .Machine$double.eps) cat(paste("\n Versions of alpha not consistent"))
  m <- moments(1:2, x)
  cat(paste("\n ** Logstable distribution has mean=", m[1],
    "  sd=", sqrt(m[2] - m[1]^2),"\n ********\n"))
}

# Replace numbers outside range -1.e-308 to +1.e308 with -Inf or Inf
bigAsInf <- function(x){
  x[x > 1.e308] <- Inf
  x[x < -1.e308] <- -Inf
  return(x)
}

tailsEstable <- function(x, stableParamObj){
  if(class(stableParamObj) != "stableParameters")stop(paste("tailsEstable:",
    "Parameter stableParamObj must be of class stableParameters"))
  storage.mode(x) <- "double"
  n <- as.double(length(x))
  temp <- .C("RtailsMSS", stableParamObj$alpha, stableParamObj$oneminusalpha,
    stableParamObj$twominusalpha, stableParamObj$location,
    stableParamObj$logscale, n, x, out1=double(n), out2=double(n),
    out3=double(n), out4=double(n), out5=double(n), out6=double(n),
    COPY=rep(c(FALSE, TRUE), c(7, 6)), PACKAGE="FMStable")
  list(density=bigAsInf(temp$out1), F=bigAsInf(temp$out3),
    righttail=bigAsInf(temp$out5), logdensity=bigAsInf(temp$out2),
    logF=bigAsInf(temp$out4), logrighttail=bigAsInf(temp$out6))
}

# Given a set of x values for a log
#  maximally skew stable distribution, find density,
#  distribution function, and right tail probability.
tailsFMstable <- function(x, stableParamObj){
  if(class(stableParamObj) != "stableParameters")stop(paste("tailsFMstable:",
    "Parameter stableParamObj must be of class stableParameters"))
  storage.mode(x) <- "double"
  # Send finite positive x values to C function
  forC <- is.finite(x) & x>0
  n <- as.double(sum(forC))
  if(n > 0) temp <- .C("RtailslogMSS", stableParamObj$alpha,
    stableParamObj$oneminusalpha, stableParamObj$twominusalpha,
    stableParamObj$location, stableParamObj$logscale, n, x[forC],
    out1=double(n), out2=double(n), out3=double(n), out4=double(n),
    out5=double(n), out6=double(n), COPY=rep(c(FALSE, TRUE), c(7, 6)),
    PACKAGE="FMStable")
  nx <- length(x)
  righttail <- F <- density <- double(nx)
  logrighttail <- logF <- logdensity <- rep(-Inf, nx)
  nas <- !is.finite(x)
  if(any(nas)){
    righttail[nas] <- F[nas] <- density[nas] <- NA
    logrighttail[nas] <- logF[nas] <- logdensity[nas] <- NA
  }
  le0 <- is.finite(x) & x <= 0
  if(any(le0)){
    righttail[le0] <- 1
    logrighttail[le0] <- F[le0] <- density[le0] <- 0
    logF[le0] <- logdensity[le0] <- -Inf
  }
  if(n > 0){
    density[forC] <- temp$out1
    F[forC] <- temp$out3
    righttail[forC] <- temp$out5
    logdensity[forC] <- temp$out2
    logF[forC] <- temp$out4
    logrighttail[forC] <- temp$out6
  }
  list(density=density, F=F, righttail=righttail, logdensity=logdensity,
    logF=logF, logrighttail=logrighttail)
}

dFMstable <- function(x, stableParamObj, log=FALSE){
  tls <- tailsFMstable(x, stableParamObj)
  if(log) tls$logdensity else tls$density
}

pFMstable <- function(x, stableParamObj, log=FALSE, lower.tail=TRUE){
  tls <- tailsFMstable(x, stableParamObj)
  if(log){ if(lower.tail) tls$logF else tls$logrighttail } else {
    if(lower.tail) tls$F else tls$righttail }
}

dEstable <- function(x, stableParamObj, log=FALSE){
  tls <- tailsEstable(x, stableParamObj)
  if(log) tls$logdensity else tls$density
}

pEstable <- function(x, stableParamObj, log=FALSE, lower.tail=TRUE){
  tls <- tailsEstable(x, stableParamObj)
  if(log){ if(lower.tail) tls$logF else tls$logrighttail } else {
    if(lower.tail) tls$F else tls$righttail }
}

# Find alpha for FMStable distribution by specifying mean,
#  standard deviation and probability of exceeding a value
fitGivenQuantile <- function(mean, sd, prob, value, tol=1.e-10){
  if(length(mean) != 1) stop("fitGivenQuantile: mean must be of length 1")
  if(length(sd) != 1) stop("fitGivenQuantile: sd must be of length 1")
  if(length(prob) != 1) stop("fitGivenQuantile: prob must be of length 1")
  if(length(value) != 1) stop("fitGivenQuantile: value must be of length 1")
  if(sd <= 0) stop("alpha.given.quantile: sd must be > 0")
  if(prob <= 0 || prob >= 1) stop("
    alpha.given.quantile: prob must be > 0 and < 1")
  if(value <= 0) stop("alpha.given.quantile: value must be > 0")
  if(mean <= 0) stop("alpha.given.quantile: mean must be > 0")
  # Check whether alpha=0 distribution would have sd larger than specified
  q <- 1-prob
  A <- mean/q
  minvar <- prob*mean*mean + q*(A-mean)^2
  if(sd*sd <= minvar) stop(paste("fitGivenQuantile: Quantile cannot be fitted.",
    "\nDistribution with probability", prob, "at zero and", q, "at", A,
    "\nhas standard deviation", sqrt(minvar), "which exceeds specified sd."))
  cv <- sd/mean
  Pdiscrepancy <- function(x) pFMstable(value/mean,
    setMomentsFMstable(mean=1, sd=cv, twominusalpha=x)) - prob
  xx <- uniroot(Pdiscrepancy, lower=0, upper=2-1.e-15, tol=tol)$root
  setMomentsFMstable(mean=mean, sd=mean*cv, twominusalpha=xx)
}

# Given alpha, find parameters of logstable distribution to matchQuartiles
matchQuartiles <- function(quartiles, alpha, oneminusalpha, twominusalpha,
    tol=1.e-10){
  if(!missing(alpha)){
    # Case where alpha is specified
    if(length(alpha) != 1) stop(
      "setMomentsFMstable: alpha must be of length 1")
    if(!missing(oneminusalpha) || !missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    case <- 1
  } else if(!missing(oneminusalpha)){
    # Case where oneminusalpha is specified
    if(length(oneminusalpha) != 1) stop(
      "setMomentsFMstable: oneminusalpha must be of length 1")
    if(!missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    case <- 2
  } else {
    if(missing(twominusalpha))stop(
      "Specify one of alpha, oneminusalpha, twominusalpha")
    # Case where twominusalpha is specified
    if(length(twominusalpha) != 1) stop(
      "setMomentsFMstable: twominusalpha must be of length 1")
    case <- 3
  }
  if(length(quartiles) != 2) stop(
    "matchQuartiles: quartiles must be of length 2")
  if(quartiles[1] <= 0 || quartiles[2] <= quartiles[1]) stop(
    "matchQuartiles: Require 0 < quartile[1] < quartile[2]")

  # First approximation to mean and sd
  av <- mean(quartiles)
  sd <- diff(quartiles)*.74
  switch(case, {
    disc <- function(param) sum((pFMstable(quartiles, setMomentsFMstable(
      exp(param[1]), exp(param[2]), alpha=alpha)) - c(.25,.75))^2)
  },{
    disc <- function(param) sum((pFMstable(quartiles, setMomentsFMstable(
      exp(param[1]), exp(param[2]), oneminusalpha=oneminusalpha)) -
      c(.25,.75))^2)
  },{
    disc <- function(param) sum((pFMstable(quartiles, setMomentsFMstable(
      exp(param[1]), exp(param[2]), twominusalpha=twominusalpha)) -
      c(.25,.75))^2)
  })

  # Find mean and sd that give match to specified quartiles
  opt <- optim(log(c(av, sd)), disc, control=list(reltol=tol*tol))
  if(opt$convergence > 0) print("matchQuartiles: Poor convergence")
  av <- exp(opt$par[1])
  sd <- exp(opt$par[2])
  switch(case,
    setMomentsFMstable(mean=av, sd, alpha=alpha),
    setMomentsFMstable(mean=av, sd, oneminusalpha=oneminusalpha),
    setMomentsFMstable(mean=av, sd, twominusalpha=twominusalpha))
}

# Find location and scale for convolution of n (not necessarily integral)
#  log maximally skew stable distributions.
# New moment generating function is original to power n
iidcombine <- function(n, stableParamObj){
  if(length(n) != 1) stop("iidcombine: n must be of length 1")
  if(class(stableParamObj) != "stableParameters") stop(paste("iidcombine:",
    "Parameter stableParamObj must be of class stableParameters"))
  if(n <=0) stop("iidcombine: n must be > 0")

  logscale.conv <- stableParamObj$logscale +log(n)/stableParamObj$alpha
  # For C parametrization, location parameters add and scales^(alpha) add
  if(stableParamObj$alpha < .5){
    logscale.conv <- stableParamObj$logscale +log(n)/stableParamObj$alpha
    location.conv <- n*stableParamObj$location
  } else {
    scale <- exp(stableParamObj$logscale)
    if (stableParamObj$oneminusalpha == 0){
      location.conv <- n*stableParamObj$location - 2/pi*(
         n*stableParamObj$logscale*scale -
              exp(logscale.conv)*logscale.conv)
    } else if (abs(stableParamObj$oneminusalpha) < 0.1) {
      location.conv <- n*stableParamObj$location +
        (exp(logscale.conv)-n*scale)/tan(.5*pi*stableParamObj$oneminusalpha)
    } else {
      tana <- if(stableParamObj$oneminusalpha > 0){
            tan(.5*pi*stableParamObj$alpha)
          } else { -tan(.5*pi*stableParamObj$twominusalpha)}
      location.conv <- n*stableParamObj$location +
        (exp(logscale.conv)-n*scale)*tana
    }
  }
  res <- list(alpha=stableParamObj$alpha,
    oneminusalpha=stableParamObj$oneminusalpha,
    twominusalpha=stableParamObj$twominusalpha, location=location.conv,
    logscale=logscale.conv, created.by="iidcombine")
  class(res) <- "stableParameters"
  return(res)
}

# Set parameters of stable distribution using any of three parametrizations
setParam <- function(alpha, oneminusalpha, twominusalpha,
    location, logscale, pm){
  if(!missing(alpha)){
    # Case where alpha is specified
    if(length(alpha) != 1) stop(
      "setParam: alpha must be of length 1")
    if(!missing(oneminusalpha) || !missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    oneminusalpha <- 1 - alpha
    twominusalpha <- 2 - alpha
  } else if(!missing(oneminusalpha)){
    # Case where oneminusalpha is specified
    if(length(oneminusalpha) != 1) stop(
      "setParam: oneminusalpha must be of length 1")
    if(!missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    alpha <- 1 - oneminusalpha
    twominusalpha <- 1 + oneminusalpha
  } else {
    if(missing(twominusalpha))stop(
      "Specify one of alpha, oneminusalpha, twominusalpha")
    # Case where twominusalpha is specified
    if(length(twominusalpha) != 1) stop(
      "setParam: twominusalpha must be of length 1")
    alpha <- 2 - twominusalpha
    oneminusalpha <- twominusalpha - 1
  }
  if(length(location) != 1) stop("setParam: location must be of length 1")
  if(length(logscale) != 1) stop("setParam: logscale must be of length 1")
  if(length(pm) != 1) stop("setParam: pm must be of length 1")
  pmcode <- rep(0:2, each=3)[ match(as.character(pm),
    c("0","S0","M", "1","S1","A", "2","CMS","C"))]
  if(is.na(pmcode)) stop(paste("setParam:  Parameterization must be",
    "specified as 0, 1 or 2; S0, S1 or CMS; or M, A or C"))
  if(alpha <= 0 || alpha > 2) stop(
    "setParam: Must have 0 < alpha <= 2")
  if(pmcode != 0)if(abs(alpha - 1) < .01)stop(
    "setParam: Only S0=M parametrization suitable near alpha=1")
  if(pmcode == 0)if(alpha < .1)stop(
    "setParam: S0=M parametrization not suitable for small alpha")
  angle <- 0.5*pi*alpha
  if(alpha < .5){
    if(pmcode == 0) {
      location <- location-tan(angle)*exp(logscale)
      logscale <- logscale-log(cos(angle))/alpha
    }
    if(pmcode == 1) logscale <- logscale-log(cos(angle))/alpha
  } else {
    tanangle <- if(alpha > 1) -tan(.5*pi*twominusalpha) else tan(angle)
    if(pmcode == 1) location <- location+tanangle*exp(logscale)
    if(pmcode == 2) {
      logscale <- logscale+log(abs(sin(.5*pi*oneminusalpha)))/alpha
      location <- location+tanangle*exp(logscale)
    }
  }
  res <- list(alpha=alpha, oneminusalpha=oneminusalpha,
    twominusalpha=twominusalpha, location=location, logscale=logscale,
    created.by=paste("setParam with pmcode =",pmcode))
  class(res) <- "stableParameters"
  return(res)
}

moments <- function(powers, stableParamObj, log=FALSE){
  if(class(stableParamObj) != "stableParameters") stop(paste("moments:",
    "Parameter stableParamObj must be of class stableParameters"))
  if(any(powers <= 0)) stop("moments: Must have powers > 0")
  logpowers <- log(powers)
  logprod <- logpowers + stableParamObj$logscale
  if(stableParamObj$alpha < 0.5){	# C parametrization
    # powers*location-(scale*powers)^alpha
    # If alpha is very small (e.g. .01) then scale might be 10^400?
    logmoment <- -powers * stableParamObj$location -
      exp(stableParamObj$alpha*logprod)
  } else {				# M = S0 parametrization
    prod <- exp(logprod)
    if(stableParamObj$oneminusalpha ==0){
      logmoment <- -powers * stableParamObj$location + logprod * 2*prod/pi
    } else {
      # Case when alpha > 0.5 but not precisely 1.
      sineEpsTerm <- sin(.5*pi*stableParamObj$oneminusalpha)
      oneminussineAlphaTerm <- 2 * sin(.25*pi*stableParamObj$oneminusalpha)^2
      logmoment <- -powers*stableParamObj$location - prod *
        (expm1(-logprod*stableParamObj$oneminusalpha)+oneminussineAlphaTerm)/
        sineEpsTerm
    }
  }
  if(log) logmoment else exp(logmoment)
}

# Compute distribution function of logstable distribution for alpha=0
pFMstable.alpha0 <- function(x, mean=1, sd=1, lower.tail=TRUE){
  if(length(mean) != 1) stop("pFMstable.alpha0: mean must be of length 1")
  if(length(sd) != 1) stop("pFMstable.alpha0: sd must be of length 1")
  # Probability of value A=mean+sd^2/mean is pA=mean/A.  Otherwise zero.
  A <- mean + sd^2/mean
  pA <- mean/A
  res <- double(length(x))
  if(lower.tail){
    res[x > 0] <- 1-pA
    res[x >= A] <- 1
  } else {
    res[x < A] <- pA
    res[x <= 0] <- 1
  }
  return(res)
}

# Value of put options under finite moment logstable distribution
putFMstable <- function(strike, paramObj, rel.tol=1.e-10){
  if(!all(is.finite(strike) & strike>0)) stop(
    "putFMstable: All values of strike must be finite and strictly positive")
  result <- double(length(strike))
  f <- function(x) pFMstable(x, paramObj)
  o <- order(strike)
  sum <- 0
  upto <- 0
  for (i in 1:length(strike)){
    s <- strike[o[i]]
    sum <- sum + integrate(f, lower=upto, upper=s, rel.tol=rel.tol)$value
    result[o[i]] <- sum
    upto <- s
  }
  lastput <- result[o[length(strike)]]
  #if(lastput+moments(1, paramObj)-s < .01 * lastput) cat(paste(
  #  "putFMstable: For large strikes, consider using callFMstable",
  #  "and the relationship\nput = call + strike - mean\n"))
  result
}

# Value of call options under finite moment logstable distribution
callFMstable <- function(strike, paramObj, rel.tol=1.e-10){
  if(!all(is.finite(strike) & strike>0)) stop(
    "callFMstable: All values of strike must be finite and strictly positive")
  result <- double(length(strike))
  f <- function(x) pFMstable(x, paramObj, lower.tail=FALSE)
  o <- order(strike, decreasing=TRUE)
  sum <- 0

  # Find where right tail is rel.tol^2 smaller than for largest strike
  pright <- pFMstable(strike[o[1]], paramObj, lower.tail=FALSE)
  if(pright <= 0) upto <- strike[o[1]] else upto <- qFMstable(
    max(min(.5, pright)*rel.tol^2, 1.e-300), paramObj, lower.tail=FALSE)

  # Work through strikes in descending order
  for (i in 1:length(strike)){
    s <- strike[o[i]]
    sum <- sum + integrate(f, lower=s, upper=upto, rel.tol=rel.tol)$value
    result[o[i]] <- sum
    upto <- s
  }
  lastcall <- result[o[length(strike)]]
  #if(lastcall-moments(1, paramObj)+s < .01 * lastcall) cat(paste(
  #  "callFMstable: For small strikes, consider using putFMstable",
  #  "and the relationship\ncall = put + mean - strike\n"))
  result
}
# Value of call and put options under finite moment logstable distribution
optionsFMstable <- function(strike, paramObj, rel.tol=1.e-10){
  if(!all(is.finite(strike) & strike>0)) stop(
    "callFMstable: All values of strike must be finite and strictly positive")
  call <- put <- double(length(strike))
  mean <- moments(1, paramObj)
  low <- strike < mean
  high <- !low
  if(any(low)){
    results <- putFMstable(strike[low], paramObj, rel.tol)
    put[low] <- results
    call[low] <- results + mean - strike[low]
  }
  if(any(high)){
    results <- callFMstable(strike[high], paramObj, rel.tol)
    call[high] <- results
    put[high] <- results + strike[high] - mean
  }
  return(list(put=put, call=call))
}

# Black-Scholes model for value of option
BSOptionValue <- function(spot, strike, expiry, volatility,
    intRate=0, carryCost=0, Call=TRUE){
  d1 <- (log(spot/strike)+(carryCost+.5*volatility*volatility)*expiry)/
    (volatility*sqrt(expiry))
  d2 <- d1 - volatility*sqrt(expiry)
  price <- if(Call){ spot*exp((carryCost-intRate)*expiry)*pnorm(d1) -
      strike*exp(-intRate*expiry)* pnorm(d2)
    } else {strike*exp(-intRate*expiry)*pnorm(-d2) -
      spot*exp((carryCost - intRate)*expiry)*pnorm(-d1)}
  return(price)
}
# tests
if(FALSE){
BSOptionValue(5, 5.1, .3, .2, Call=FALSE)
BSOptionValue(5, 5.1, .3, .2)
BSOptionValue(5, 5.2, .3, .2)
BSOptionValue(5, 5.3, .3, .2)
BSOptionValue(5, 5, .4, .3)
}

# Find implied volatilities to make prices consistent with Black-Scholes model
# May give vectors of strike, expiry and price
ImpliedVol <- function(spot, strike, expiry, price, intRate=0, carryCost=0,
    Call=TRUE, ImpliedVolLowerBound=.01, ImpliedVolUpperBound=1, tol=1.e-9){
  Nprice <- length(price)
  if(!(length(strike) %in% c(1,Nprice)))stop(
    "ImpliedVol: strike must be a scalar or of same length as price")
  if(!(length(expiry) %in% c(1,Nprice)))stop(
    "ImpliedVol: expiry must be a scalar or of same length as price")
  f <- function(vol){
    BSOptionValue(spot,Si,Ei,vol,intRate,carryCost,Call) - Pi
  }
  impVol <- rep(NA,Nprice)
  Si <- strike[1]
  Ei <- expiry[1]
  for (i in 1:Nprice){
    if(length(strike) > 1) Si <- strike[i]
    if(length(expiry) > 1) Ei <- expiry[i]
    Pi <- price[i]
    inrange <- f(ImpliedVolLowerBound)*f(ImpliedVolUpperBound) < 0
    if(is.na(inrange)) inrange <- FALSE
    if(inrange){impVol[i] <- uniroot(f,lower=ImpliedVolLowerBound,
      upper=ImpliedVolUpperBound,tol=tol)$root}
  }
  return(impVol)
}
# tests
if(FALSE){
  ImpliedVol(5,5.1,.3,0.1741752)
  ImpliedVol(5,5.1,.3,0.274152,Call=FALSE)
  ImpliedVol(5,rep(5.1,3),.3,c(0.1741752,.17,-.1))
  ImpliedVol(5,c(5.1,5.3),c(.3,.4),c(0.1741752,.17))
  ImpliedVol(5,5.2,c(.3,.4),c(0.1741752,.17))
  ImpliedVol(5,50,.7,0.001)
}

# Find parameters of lognormal distibution to match given mean and sd
lnorm.param <- function(mean, sd){
  cv <- sd/mean
  varlog <- log1p(cv*cv)
  return(list(meanlog=log(mean)-.5*varlog, sdlog=sqrt(varlog)))
}

######################################################
#Compute distribution function of stable distribution by numerical integration
# Do not bother with alpha=1, since can set oneminusalpha to be, say, 1.e-20.
# Leave out special cases alpha=2 and where distribution function is 0 or 1.

# Tamper with calculation of integrand to give good precision
# Note that the functions which compute integrands are always called with
#  a vector of length 21 in R.
# This tampering is based on the values which will be used when integrating
#  over a range -1 to 1.  These were obtained using the commands:
#    integrand <- function(theta){print(theta, digits=16);theta}
#    integrate(integrand, -1, 1)

####### Start of function to compute tail probabilities ########
tailsGstable <- function(x, logabsx, alpha, oneminusalpha, twominusalpha,
     beta, betaplus1, betaminus1, parametrization, lower.tail=TRUE){

  if(length(x) != 1) stop("tailsGstable: x must be of length 1")
  if(length(logabsx) != 1) stop("tailsGstable: logabsx must be of length 1")

  if(!missing(alpha)){
    # Case where alpha is specified
    if(length(alpha) != 1) stop("tailsGstable: alpha must be of length 1")
    if(!missing(oneminusalpha) || !missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    oneminusalpha <- 1 - alpha
    twominusalpha <- 2 - alpha
  } else if(!missing(oneminusalpha)){
    # Case where oneminusalpha is specified
    if(length(oneminusalpha) != 1) stop(
      "tailsGstable: oneminusalpha must be of length 1")
    if(!missing(twominusalpha)) stop(
      "Only specify one of alpha, oneminusalpha, twominusalpha")
    alpha <- 1 - oneminusalpha
    twominusalpha <- 1 + oneminusalpha
  } else {
    if(missing(twominusalpha))stop(
      "Specify one of alpha, oneminusalpha, twominusalpha")
    # Case where twominusalpha is specified
    if(length(twominusalpha) != 1) stop(
      "tailsGstable: twominusalpha must be of length 1")
    alpha <- 2 - twominusalpha
    oneminusalpha <- twominusalpha - 1
  }

  if(!missing(beta)){
    # Case where beta is specified
    if(length(beta) != 1) stop("tailsGstable: beta must be of length 1")
    if(!missing(betaplus1) || !missing(betaminus1)) stop(
      "Only specify one of beta, betaplus1, betaminus1")
    betaplus1 <- beta + 1
    betaminus1 <- beta - 1
  } else if(!missing(betaplus1)){
    # Case where betaplus1 is specified
    if(length(betaplus1) != 1) stop(
      "tailsGstable: betaplus1 must be of length 1")
    if(!missing(betaminus1)) stop(
      "Only specify one of beta, betaplus1, betaminus1")
    beta <- betaplus1 - 1
    betaminus1 <- betaplus1 - 2
  } else {
    if(missing(betaminus1))stop("Specify one of beta, betaplus1, betaminus1")
    # Case where betaminus1 is specified
    if(length(betaminus1) != 1) stop(
      "tailsGstable: betaminus1 must be of length 1")
    beta <- betaminus1 + 1
    betaplus1 <- betaminus1 + 2
  }

  INT.PTS <- c(0, -0.9739065285171717, 0.9739065285171717, -0.8650633666889845, 
    0.8650633666889845, -0.6794095682990244, 0.6794095682990244, 
    -0.4333953941292472, 0.4333953941292472, -0.1488743389816312, 
    0.1488743389816312, -0.9956571630258081, 0.9956571630258081, 
    -0.9301574913557082, 0.9301574913557082, -0.7808177265864169, 
    0.7808177265864169, -0.5627571346686047, 0.5627571346686047, 
    -0.2943928627014602, 0.2943928627014602)

  POWERSOFHALF <- c(1, cumprod(rep(.5, 60)))
  HPI <- pi/2
  LARGE <- 999

  integrand <- function(p){
    # p is the proportion of the integration range
    # Convention: global variables have names in capital letters
    # Test input data sometimes useful:  p <- INT.PTS*.5+.5
    # Find n such that integration range is 2^(-n) of full range
    n <- round(log(p[3]-p[2])/log(.5))
    if(n>53) stop("Integration interval too small")
    # Mid-point will have no rounding error, provided n < about 54
    # Check by printing:   1-POWERSOFHALF-1
    pmid <- p[1]
    qmid <- 1-pmid
    factor <- POWERSOFHALF[n+2]
    ps <- pmid+INT.PTS*factor
    qs <- qmid-INT.PTS*factor
  
    # Compute linear combinations of theta to high precision
    theta <- qs*THETA.LOW +ps*THETA.HIGH
    pctheta <- qs*PCTHETA.LOW +ps*PCTHETA.HIGH
    nctheta <- qs*NCTHETA.LOW +ps*NCTHETA.HIGH
    arg1 <- qs*ARG1.LOW +ps*ARG1.HIGH
    pcarg1 <- qs*PCARG1.LOW +ps*PCARG1.HIGH
    ncarg1 <- qs*NCARG1.LOW +ps*NCARG1.HIGH
    arg2 <- qs*ARG2.LOW +ps*ARG2.HIGH
    pcarg2 <- qs*PCARG2.LOW +ps*PCARG2.HIGH
    ncarg2 <- qs*NCARG2.LOW +ps*NCARG2.HIGH
  
    # Compute integrand to high precision
    # Integrand is exp(-w) where w=cos(theta-ALPHA*(theta-PHI0))*
    #    sin(ALPHA*(theta-PHI0))^(ALPHA/(1-ALPHA))*
    #    cos(theta)^(-1/(1-ALPHA))*X^(-ALPHA/(1-ALPHA)) 
    # Sometimes both sin(ALPHA*(theta-PHI0)) and X are negative. This is OK.
  
    # pctheta is theta+pi/2 
    # nctheta is theta-pi/2 
    # arg1 is alpha*(theta-phi0) 
    # pcarg1 is alpha*(theta-phi0) + pi/2
    # ncarg1 is alpha*(theta-phi0) - pi/2
    # arg2 is theta-alpha*(theta-phi0) 
    # pcarg2 is theta-alpha*(theta-phi0)+pi/2 
    # ncarg2 is theta-alpha*(theta-phi0)-pi/2 
  
    term1 <- l2 <- l3 <- rep(0, length(p))
    # Compute cos(arg2) as sine of one of
    #  theta+ALPHA*(theta-PHI0)+HPI and theta+ALPHA*(theta-PHI0)-HPI
    first <- abs(pcarg2) < abs(ncarg2)
    term1[first] <- sin(pcarg2[first])
    term1[!first] <- sin(-ncarg2[!first])
    # When arg1 is not near to unity, (say abs(arg1) > 1), 
    #  use sin(arg1)=1-2*sin(a)^2 where a is 0.5*(PI-arg1) or 0.5(arg1+PI)
    first <- abs(pcarg1) < abs(ncarg1)
    a <- ncarg1
    a[first] <- pcarg1[first]
    first <- abs(arg1) > 1
    l2[first] <- log1p(-2*sin(.5*a[first])^2)
    l2[!first] <- log(abs(sin(arg1[!first])))
    # When abs(theta) < 0.7, compute cos(theta) as 1-2*sin(0.5*theta)^5
    #  otherwise use sin(b) where b is theta+HPI or theta-HPI
    first <- abs(theta) < 0.7
    secondtest <- abs(pctheta) < abs(nctheta)
    second <- !first & secondtest
    third <- !first & !secondtest
    l3[first] <- log1p(-2*sin(.5*theta[first])^2)
    l3[second] <- log(sin(pctheta[second]))
    l3[third] <- log(sin(-nctheta[third]))
    w <- term1*exp((l3+ALPHA*(LOGABSX-l2))/ALPHAMINUS1)
    # Note that if w is infinite, we still want exp(-w) to give zero.
    # The next line is a safety precaution, knowing that exp(-LARGE) will give 0
    w[!is.finite(w)] <- LARGE
    # If subinterval is at an extreme of the integration region, 
    #  then ensure that the value for W nearest to the limit of the range
    #  of integration is within 1 of the value at the limit.
    # This is meant to stop the numerical procedure from completely missing
    #  the region where the integrand > 0, or is < 1.
    if(pmid==factor) w[2] <- if(LEFTLIMIT < RIGHTLIMIT){
      min(w[2], LEFTLIMIT+1) } else { max(w[2], LEFTLIMIT-1)}
    if(qmid==factor) w[3] <- if(LEFTLIMIT > RIGHTLIMIT){
      min(w[3], RIGHTLIMIT+1) } else { max(w[3], RIGHTLIMIT-1)}
    # Use complementary function if this will give more precise tail prob.
    if(COMPLEMENT){-expm1(-w)} else {exp(-w)}
  }

  # If x is outside the allowed exponents of real numbers, use +1 or -1 being
  #  careful to give correct sign.  Always input logabsx to full precision.
#  if(alpha > 1.5){alpha <- 2-twominusalpha; oneminusalpha <- twominusalpha-1} else
#   {if(alpha < .5){twominusalpha <- 2-alpha;oneminusalpha <- 1-alpha} else
#      {twominusalpha <- oneminusalpha+1; alpha <- 1-oneminusalpha}}
  k <- if(oneminusalpha>0){alpha}else{twominusalpha}
  # Change parametrization, if necessary
  if(parametrization < 2){
    a <- HPI*k
    # Compute tan(a) to high precision
    if(a<.8){tana <- tan(a)} else {tana <- 1/tan(HPI*abs(oneminusalpha))}
    btn <- beta*tana
    if(beta > 0.8){
      betaminus1 <- atan(betaminus1*tana/(1+btn*tana))/a
      beta <- betaminus1+1
      betaplus1 <- betaminus1+2
    } else { if(beta < -0.8){
      betaplus1 <- atan(betaplus1*tana/(1-btn*tana))/a
      beta <- betaplus1-1
      betaminus1 <- betaplus1-2
    } else {
      beta <- atan(btn)/a
      betaplus1 <- beta+1
      betaminus1 <- beta-1
      }
    }
    # If alpha>1 then change sign of beta
    if(oneminusalpha < 0 && parametrization < 2){
      beta <- -beta
      btn <- -btn
      temp <- betaminus1
      betaminus1 <- -betaplus1
      betaplus1 <- -temp
      }
    scale <- exp(-.5/alpha*log1p(btn^2)) # =(1+btn**2)**(-.5/alpha)
  }

  # Zolotarev's M = Nolan S0 (which is not used near alpha=0 or for extreme x)
  if(parametrization == 0){
    if(abs(btn) > 10*abs(x)){
      logabsx <- log1p(x/btn)-.5/alpha * log1p(1/(btn*btn))-
        oneminusalpha/alpha * log(abs(btn))
      x <- sign(btn)*exp(logabsx)
    } else {
      x <- (x+btn)*scale
      logabsx <- log(abs(x))
    }
  }

  # Zolotarev's A = Nolan S1
  if(parametrization == 1){
    x <- x*scale
    logabsx <- logabsx +log(scale)
  }

  # Zolotarev's C = Chambers, Mallows and Stuck
  if(parametrization == 2){
    scale <- 1.
  }

  # Code is not designed to work if x=0; alpha >= 2; alpha <=0; alpha=1;
  #  alpha<1, beta=1 & x<0; or alpha<1, beta=-1 & x>0
  if( x==0 || twominusalpha <= 0 || alpha <= 0 || oneminusalpha ==0 ||
    oneminusalpha > 0 && ((betaminus1==0 && x<0) || (betaplus1==0 && x>0))){
    print("Check parameters input to function tailsGstable")}
  phi0 <- -HPI*beta*k/alpha
  possibleWlimit <- abs(oneminusalpha)*abs(x/alpha)^(-alpha/oneminusalpha)
  if(oneminusalpha>0){  #Case when alpha<1
    if(x < 0){
  # THETA.LOW is the lower integration limit
  # PCTHETA.LOW is theta+pi/2 at the lower integration limit
  # NCTHETA.LOW is theta-pi/2 at the lower integration limit
  # ARG1.LOW is alpha*(theta-phi0) at the lower integration limit
  # PCARG1.LOW is alpha*(theta-phi0) + pi/2 at the lower integration limit
  # NCARG1.LOW is alpha*(theta-phi0) - pi/2 at the lower integration limit
  # ARG2.LOW is theta-alpha*(theta-phi0) at the lower integration limit
  # PCARG2.LOW is theta-alpha*(theta-phi0)+pi/2 at the lower integration limit
  # NCARG2.LOW is theta-alpha*(theta-phi0)-pi/2 at the lower integration limit

      THETA.LOW <- -HPI
      ARG1.LOW <- HPI*alpha*betaminus1
      THETA.HIGH <- phi0
      ARG1.HIGH <- 0
      LEFTLIMIT <- LARGE
      RIGHTLIMIT <- if(betaplus1==0){possibleWlimit}else{0}
    } else {
      THETA.LOW <- phi0
      ARG1.LOW <- 0
      THETA.HIGH <- HPI
      ARG1.HIGH <- HPI*alpha*betaplus1
      LEFTLIMIT <- if(betaminus1==0){possibleWlimit}else{0}
      RIGHTLIMIT <- LARGE
    }
    # Use complementary integrand whenever right tail
    # required or integrating over less than -HPI to HPI
    singleregion <- (phi0+HPI) == 0 || (phi0-HPI) == 0
    COMPLEMENT <- !(lower.tail && singleregion)
  } else {    #Case when alpha>1
    if(x < 0){
      THETA.LOW <- -HPI
      ARG1.LOW <- HPI*(-alpha+beta*twominusalpha)
      THETA.HIGH <- phi0
      ARG1.HIGH <- 0
      LEFTLIMIT <- if(betaplus1>0){0}else{possibleWlimit}
      RIGHTLIMIT <- LARGE
    } else {
      THETA.LOW <- phi0
      ARG1.LOW <- 0
       THETA.HIGH <-  HPI
      ARG1.HIGH <- HPI*(alpha+beta*twominusalpha)
      LEFTLIMIT <- LARGE
      RIGHTLIMIT <- if(betaminus1<0){0}else{possibleWlimit}
    }
    COMPLEMENT <- FALSE
  }

  PCTHETA.LOW <- THETA.LOW+HPI  
  NCTHETA.LOW <- THETA.LOW-HPI
  PCTHETA.HIGH <- THETA.HIGH+HPI  
  NCTHETA.HIGH <- THETA.HIGH-HPI
  PCARG1.LOW <- ARG1.LOW+HPI
  NCARG1.LOW <- ARG1.LOW-HPI
  PCARG1.HIGH <- ARG1.HIGH+HPI
  NCARG1.HIGH <- ARG1.HIGH-HPI
  ARG2.LOW <- THETA.LOW-ARG1.LOW
  PCARG2.LOW <- PCTHETA.LOW-ARG1.LOW
  NCARG2.LOW <- NCTHETA.LOW-ARG1.LOW
  ARG2.HIGH <- THETA.HIGH-ARG1.HIGH
  PCARG2.HIGH <- PCTHETA.HIGH-ARG1.HIGH
  NCARG2.HIGH <- NCTHETA.HIGH-ARG1.HIGH

  LOGABSX <- logabsx
  ALPHA <- alpha
  ALPHAMINUS1 <- -oneminusalpha

  result <- integrate(integrand, 0, 1, subdivisions=1000, stop.on.error=F, 
    abs.tol=1.e-50, rel.tol=1.e-14)

  contribution <- result$value*(THETA.HIGH-THETA.LOW)/pi

  if(x<0){
    lefttail <- contribution; righttail <- 1-lefttail
  }else{
    righttail <- contribution; lefttail <- 1-righttail}

  # If alpha<1 & lower.tail & singleregion then swap tails around
  swap <- if(oneminusalpha < 0) {FALSE} else {lower.tail && singleregion}
  if(swap){temp <- lefttail; lefttail <- righttail; righttail <- temp}
  list(left.tail.prob=lefttail, right.tail.prob=righttail, 
    est.error=result$abs.error, message=result$message)
}

# Compute quantiles of standard stable distributions
# Prototype code to be converted into C later.
# Logarithms of both tail probabilities assumed to be provided.
# Return both z and logz, with logz being meaningless if z <= 0
# Deals with only one quantile at a time.
# The stable distribution object passed must have standard location and scale
standardStableQuantile <- function(logp, logq, obj, browse=FALSE){
  if(length(logp) > 1 || length(logq) > 1) stop(
    "standardStableQuantile: logp and logq must both be scalars")
  if(abs(1-exp(logp)-exp(logq)) > 1.e-14) stop(
    "standardStableQuantile: logp and logq must be consistent")
  lower.tail <- logp < logq
  if(lower.tail){ if(logp == -Inf) stop(
    "standardStableQuantile: logp input as -Inf")
  } else { if(logq == -Inf) stop(
    "standardStableQuantile: logq input as -Inf")
  }

  if(browse)browser()

  # If alpha==2, return quantile immediately
  if(obj$twominusalpha ==0){
    z <- sqrt(2)* if(lower.tail) qnorm(logp, log.p=TRUE) else qnorm(logq,
      log.p=TRUE, lower.tail=FALSE)
    return(z)
  }

  # First stage: Find reasonable first approximation
  if(obj$alpha < .5){
    uselogz <- TRUE
    if(lower.tail){
      # First approximation to xi
      xi <- -logp
      # One Newton-Raphson iteration
      xi <- xi -(logp + .5*log(2*pi*obj$alpha*xi) +xi)/(.5/xi + 1)
      # Find logz for this xi
      logz <- log(obj$alpha)-obj$oneminusalpha/obj$alpha*log(xi/obj$oneminusalpha)
    } else {
      logCalpha <- log(gamma(obj$alpha)*sin(pi*obj$alpha)/pi)
      logz <- (logq-logCalpha) *(-1/obj$alpha)
    }  
  } else if (obj$alpha > 1.7 && min(logp, logq) > -20){
    # Case where quantile of normal distribution gives a useful starting point

# Speed would be improved by giving a better approximate root for
#    moderate alpha, moderate p

    uselogz <- FALSE
    z <- sqrt(2)* qnorm(logp, log.p=TRUE)
  } else if (!lower.tail){
    uselogz <- TRUE
    # Case of right tail for alpha >= .5
    logy <- (log(2*gamma(obj$alpha + 1)* sin(.5 * pi *obj$twominusalpha)/
      (pi*obj$alpha)) - logq)/obj$alpha
    if(obj$oneminusalpha ==0 || logy > 100) logz <- logy else{
      logz <- logy + log1p(expm1(obj$oneminusalpha * logy)/
       (tan(pi*obj$oneminusalpha/2) * exp(logy)))
    }
  } else {
    # Case of left tail for alpha >= 0.5
    uselogz <- FALSE
    # First approximation to xi
    xi <- -logp
    # One Newton-Raphson iteration
    xi <- xi -(logp + .5*log(2*pi*obj$alpha*xi) +xi)/(.5/xi + 1)
    # Find z for this xi
    if(obj$oneminusalpha == 0){ # i.e. alpha = 1
      z <- -(1+log(.5*pi*xi))/(.5*pi)
    } else if(obj$alpha < 1.9){ # i.e. not near alpha==2
      ang <- .5*pi*obj$oneminusalpha
      cosa <- if(obj$alpha < 1.5) cos(ang) else sin(.5*pi*obj$twominusalpha)
      z <- expm1(  log(obj$alpha/cosa)+
        (obj$oneminusalpha/obj$alpha)*log(obj$oneminusalpha/
        (xi* sin(ang)))  )/tan(ang)
    } else  {	#i.e. when alpha is near 2
      ang <- .5*pi*obj$twominusalpha
      z <- -expm1(  log(obj$alpha/sin(ang))+
        (obj$oneminusalpha/obj$alpha)*log(-obj$oneminusalpha/
        (xi* cos(ang)))  )*tan(ang)
    }
  }

  # Temporary intermediate stage:  Check that approximation is OK
  if(browse){
    if(uselogz){
      if(logz > 700) return(Inf)
      z <- exp(logz)
    }
    if(!is.finite(z)) return(Inf)
    frac.error <- if(lower.tail){
      pEstable(z, obj, lower.tail=TRUE)/exp(logp) - 1
    } else {
      pEstable(z, obj, lower.tail=FALSE)/exp(logq) - 1
    }
    print(paste("Fractional error:", frac.error), quote=FALSE)
  }

  # Second stage: Refine approximation by Newton-Raphson
  if(uselogz){
    if(lower.tail){
      adjustmentlog <- function(logz){
        tls <- tailsEstable(exp(logz), obj)
        (logp - tls$logF)*exp(-logz + tls$logF - tls$logdensity)
      }
    } else {
      adjustmentlog <- function(logz){
        tls <- tailsEstable(exp(logz), obj)
        (tls$logrighttail - logq)*exp(-logz + tls$logrighttail - tls$logdensity)
      }
    }
    for (i in 1:10){
      adj <- adjustmentlog(logz)
      logz <- logz + adj
      if(logz > 700) return(Inf)
      if(abs(adj) < 1.e-8* max(1, abs(logz))) break
    }
    if(browse) print(paste("Log z scale:",i,"iterations used"), quote=FALSE)
    z <- exp(logz)
  } else {
    if(lower.tail){
      adjustment <- function(z){
        tls <- tailsEstable(z, obj)
        (logp - tls$logF)*exp(tls$logF - tls$logdensity)
      }
    } else {
      adjustment <- function(z){
        tls <- tailsEstable(z, obj)
        (tls$logrighttail - logq)*exp(tls$logrighttail - tls$logdensity)
      }
    }
    for (i in 1:10){
      adj <- adjustment(z)
      z <- z + adj
      if(abs(adj) < 1.e-8 * max(1, abs(z))) break
      if(!is.finite(z) || abs(z) > 1.e300) return(NA)
    }
    if(browse) print(paste("Using z scale:",i,"iterations used"), quote=FALSE)
    logz <- if(z>0) log(z) else NA
  }
  return(z)
}

# Quantile-finding using R function standardStableQuantile
# Find quantiles of a stable distribution
qEstable <- function(p, stableParamObj, log=FALSE, lower.tail=TRUE){
  if(class(stableParamObj) != "stableParameters")stop(paste("qEstable:",
    "Parameter stableParamObj must be of class stableParameters"))
  if(!all(is.finite(p))) stop("qEstable: Components of p must all be finite")
  if(log){
    logp <- p
    if(any(logp >= 0)) stop("qEstable: Logs of probabilities must be < 0")
    complement <- log1p(-exp(p))
  } else {
    if(any(p >= 1)) stop("qEstable: Require p < 1")
    logp <- log(p)
    complement <- log1p(-p)
  }
  scale <- exp(stableParamObj$logscale)
  results <- double(length(p))
  standard <- stableParamObj
  standard$location <- 0
  standard$logscale <- 0
  standard$created.by <- "function qEstable"
  for (i in 1:length(p)){
    if(lower.tail){
      r <- standardStableQuantile(logp[i], complement[i], standard)
    } else {
      r <- standardStableQuantile(complement[i], logp[i], standard)
    }
    results[i] <- r * scale + stableParamObj$location
  }
  return(results)
}
# Find quantiles of a log stable distribution
qFMstable <- function(p, stableParamObj, lower.tail=TRUE){
  if(!all(is.finite(p))) stop("qFMstable: Components of p must all be finite")
  if(class(stableParamObj) != "stableParameters")stop(paste("qFMstable:",
    "Parameter stableParamObj must be of class stableParameters"))
  return(exp(-qEstable(p, stableParamObj, lower.tail=!lower.tail, log=FALSE)))
}
