#
# bspec: Bayesian spectral inference.
#
# S3 Version.
#
# http://cran.r-project.org/package=bspec
#
#
# see also:
#
#   Roever, C., Meyer, R. and Christensen, N. (2011):
#   Modelling coloured residual noise in gravitational-wave signal processing.
#   Classical and Quantum Gravity, 28(1):015010
#   URL http://dx.doi.org/10.1088/0264-9381/28/1/015010
#
#   Arxiv preprint 0804.3853 [stat.ME]
#   URL http://arxiv.org/abs/0804.3853
#
#
# and:
#
#   Roever, C. (2011):
#   A Student-t based filter for robust signal detection.
#   Physical Review D, 84(12):122004
#   URL http://dx.doi.org/10.1103/PhysRevD.84.122004
#
#   Arxiv preprint 1109.0442 [physics.data-an]
#   URL http://arxiv.org/abs/1109.0442
#

print.bspec <- function(x, ...)
{
  if (x$two.sided) cat(" 'bspec' posterior spectrum (two-sided).\n")
  else cat(" 'bspec' posterior spectrum (one-sided).\n")
  if (is.element("temperature", names(x))){
    cat(paste(" (tempered: T = ",x$temperature,")\n",sep=""))
  }
  cat(paste(" frequency range     : ",
            as.character(signif(x$freq[1],4)),
            "--",
            as.character(signif(x$freq[length(x$freq)],4)),
            "\n",sep=""))
  cat(paste(" number of parameters: ",length(x$scale),"\n",sep=""))
  expect <- expectation(x)
  cat(" finite expectations : ")
  if (all(is.finite(expect))) cat("all\n")
  else if (any(is.finite(expect))) cat("some\n")
  else cat("none\n")
  vari <- variance(x)
  cat(" finite variances    : ")
  if (all(is.finite(vari))) cat("all\n")
  else if (any(is.finite(vari))) cat("some\n")
  else cat("none\n")
  cat(" call: "); print(x$call)
  invisible(x)
}



print.bspecACF <- function(x, ...)
{
  if (x$type=="correlation")
    cat(paste("Autocorrelations of bspec object '",x$bspec,"', by lag\n", sep=""))
  else
    cat(paste("Autocovariances of bspec object '",x$bspec,"', by lag\n", sep=""))
  printvec <- x$acf
  names(printvec) <- x$lag
  print(printvec)
  invisible(x)
}



bspec.default <- function(x, priorscale=1, priordf=0,
                          intercept=TRUE, two.sided=FALSE, ...)
# x                   : time series object (uni- or multivariate)
# priorscale, priordf : a vector or a function of frequency
# intercept           : flag indicating whether to include zero frequency
# two.sided           : flag indicating whether to refer to one- or two-sided spectrum.
#                       only effect within _this_ function is interpretation of prior scale.
#                       Note that the 'two.sided' flag is "inherited" by other functions
#                       via their default arguments. See e.g. 'expectation.bspec()'.
{
  # some initial checks:
  if (is.vector(x) | (is.ts(x) & !is.mts(x)))
    N <- length(x)
  else if (is.mts(x) | is.matrix(x) | is.data.frame(x)) 
    N <- nrow(x)
  else warning("incompatible argument 'x'")
  FTlength <- (N %/% 2) + 1   # size of (nonredundant) FT output
  Neven <- ((N %% 2) == 0)    # indicator for even N
  stopifnot(is.function(priorscale) ||
            ((length(priorscale)==1) |
            (intercept & (length(priorscale)==FTlength)) |
            ((!intercept) & (length(priorscale)==(FTlength-1)))))
  stopifnot(is.function(priordf) ||
            ((length(priordf)==1) |
            (intercept & (length(priordf)==FTlength)) |
            ((!intercept) & (length(priordf)==(FTlength-1)))))
  stopifnot(is.function(priorscale) || (all(is.finite(priorscale)) & all(priorscale>0)),
            is.function(priordf) || (all(is.finite(priordf)) & all(priordf>=0)))
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("argument 'x' is not a time-series object, default conversion 'as.ts(x)' applied.")
  }
  deltat <- 1 / tsp(x)[3]
  deltaf <- 1 / (N*deltat)
  t0 <- tsp(x)[1]           # time stamp corresponding to 1st observation
  kappa <- c(0, rep(1,FTlength-2), ifelse(Neven,0,1))
  # 1-D case:
  if (!is.mts(x)) {
    # Fourier transform:
    y <- fft(as.vector(x))
    # (`fft()' yields unnormalised FT)
    nonredundant <- 1:FTlength
    # vector of (N/2 + 1) cosine coefficients:
    a <- (1+kappa) * sqrt(deltat/N) * Re(y[nonredundant])
    # vector of (N/2 + 1) sine coefficients:
    b <- -(1+kappa) * sqrt(deltat/N) * Im(y[nonredundant])
    datassq <- a^2 + b^2
    datadf <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2))
  }
  # multidimensional case:
  else {
    # Fourier transform:
    #y <- apply(x, 2, fft)
    y <- mvfft(x)
    # (`fft()' yields unnormalised FT)
    nonredundant <- 1:FTlength
    # vector of (N/2 + 1) cosine coefficients:
    a <- sqrt(deltat/N) * Re(y[nonredundant,])
    for (i in 1:ncol(a)) a[,i] <- (1+kappa)*a[,i]
    # vector of (N/2 + 1) sine coefficients:
    b <- -sqrt(deltat/N) * Im(y[nonredundant,])
    for (i in 1:ncol(b)) b[,i] <- (1+kappa)*b[,i]
    datassq <- apply(a^2 + b^2, 1, sum)
    datadf <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2)) * ncol(x)
  }
  # vector of corresponding frequencies:
  freq <- (0:(FTlength-1)) * deltaf
  
  # set up a-priori scale:
  if (is.function(priorscale))
    priorscalevec <- priorscale(freq)
  else if (length(priorscale)==1)
    priorscalevec <- rep(priorscale, FTlength)
  else if (!intercept)
    priorscalevec <- c(0, priorscale)
  else priorscalevec <- priorscale
  if (two.sided){
    #  /!\  different interpretation of prior scale  /!\
    priorscalevec <- priorscalevec * (1+kappa)
  }
  
  # set up a-priori degrees-of-freedom:
  if (is.function(priordf))
    priordfvec <- priordf(freq)
  if (length(priordf)==1)
    priordfvec <- rep(priordf, FTlength)
  else if (!intercept)
    priordfvec <- c(0, priordf)
  else priordfvec <- priordf
  
  #arg <- seq(from=tsp(x)[1], to=tsp(x)[2], le=500)
  #trigo <- arg*0
  #for (i in 1:length(a))
  # trigo<-trigo+(a[i]*cos(2*pi*freq[i]*(arg-t0)) + b[i]*sin(2*pi*freq[i]*(arg-t0)))
  #trigo <- trigo / sqrt(N*deltat)
  #plot(x,type="b")
  #lines(arg, trigo, col="red")

  if (! intercept) {
    freq <- freq[-1]
    priorscalevec <- priorscalevec[-1]
    priordfvec <- priordfvec[-1]
    datassq <- datassq[-1]
    datadf <- datadf[-1]
    kappa <- kappa[-1]
  }  
  
  # determine posterior distribution's parameters:
  # (posterior distn. of 1-sided spectrum S_1(f_j) ==  sigmasquared_j!)
  #  -->  S_2(f_j)  =  sigmasquared_j / (1+kappa(j))  =  S_1(f_j) / (1+kappa(j))
  scale <- (priordfvec*priorscalevec + datassq) / (priordfvec + datadf)
  df <- datadf + priordfvec
  
  result <- list(freq = freq,
                 scale = scale,
                 df = df,
                 priorscale = priorscalevec,
                 priordf = priordfvec,
                 datassq = datassq,
                 datadf = datadf,
                 N = N, deltat = deltat, deltaf = deltaf,
                 start = t0,
                 call = match.call(expand.dots=FALSE),
                 two.sided = two.sided)
  class(result) <- "bspec"
  return(result)
}



expectation.bspec <- function(x, two.sided=x$two.sided, ...)
{
  expect <- rep(0, length(x$scale))
  finite <- x$df > 2
  expect[finite] <- x$scale[finite] * (x$df[finite]/(x$df[finite]-2))
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
    if (x$freq[1] != 0) index <- index-1
    expect[index] <- expect[index] / 2
  }
  expect[!finite] <- Inf
  return(expect)
}



variance.bspec <- function(x, two.sided=x$two.sided, ...)
{
  vari <- rep(0, length(x$scale))
  finite <- (x$df > 4)
  vari[finite] <- (2*x$df[finite]^2*x$scale[finite]^2) / ((x$df[finite]-2)^2 * (x$df[finite]-4))
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
    if (x$freq[1] != 0) index <- index-1
    vari[index] <- vari[index] / 4
  }
  vari[!finite] <- Inf
  return(vari)
}



sample.bspec <- function(x, size=1, two.sided=x$two.sided, ...)
{
  rinvchisq <- function(n=1, df=1, scale=1)
  # sample n times from an Inverse-Chi-Squared distribution
  # with scale "scale" and degrees-of-freedom "df".
  {
    return(scale * (df/rchisq(n=n,df=df)))
  }
  sigmasq <- matrix(rinvchisq(n = size*length(x$freq),
                              df = x$df, scale = x$scale),
                    ncol=size)
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, x$N %/% 2 + 1)
    if (x$freq[1] != 0) index <- index-1
    sigmasq[index,] <- sigmasq[index,] / 2
  }
  if (size==1) sigmasq <- as.vector(sigmasq)
  return(sigmasq)
}



quantile.bspec <- function(x, probs = c(0.025,0.5,0.975),
                           two.sided=x$two.sided, ...)
{
  qinvchisq <- function(p=0.5, df=1, scale=1)
  {
    return((df*scale)/qchisq(1-p,df=df))
  }
  sigmasq <- matrix(qinvchisq(p = rep(probs,each=length(x$freq)),
                              df = x$df, scale = x$scale),
                    ncol=length(probs))
  if (two.sided){
    Neven <- ((x$N %% 2) == 0)
    index <- 2:ifelse(Neven, x$N %/% 2, x$N %/% 2 + 1)
    if (x$freq[1] != 0) index <- index-1
    sigmasq[index,] <- sigmasq[index,] / 2
  }
  if (length(probs)==1) sigmasq <- as.vector(sigmasq)
  else colnames(sigmasq) <- as.character(signif(probs))
  return(sigmasq)
}



plot.bspec <- function(x, two.sided=x$two.sided, ...)
{
  oldpar <- par(no.readonly=TRUE)
  quant <- quantile(x, probs=c(0.025,0.5,0.975),
                    two.sided=two.sided)
  plot(range(x$freq), range(quant),type="n", axes=FALSE,
       log="y", xlab="", ylab="", ...)
  matlines(rbind(x$freq,x$freq), pch=rep(NA,2),
           t(quant[,c(1,3)]), col="darkgrey", lty="solid")
  points(x$freq, quant[,2], pch=5)
  expect <- expectation(x, two.sided=two.sided)
  if (any(is.finite(expect)))
    points(x$freq, expect, pch=20)
  axis(1)
  mtext(expression("frequency "*italic(f)),
        side=1, line=par("mgp")[1], cex=par("cex.lab"))
  axis(2)
  if (two.sided)
    mtext(expression("(two-sided) posterior spectrum "*italic(S(f))),
          side=2, line=par("mgp")[1], cex=par("cex.lab"))
  else
    mtext(expression("(one-sided) posterior spectrum "*italic(S(f))),
          side=2, line=par("mgp")[1], cex=par("cex.lab"))
  sqrtrange <- par("usr")[3:4]/2
  ticks <- axTicks(side=4, usr=sqrtrange, log=par("ylog"))
  axlabel <- NULL
  for (i in 1:length(ticks))
    axlabel <- c(axlabel, eval(bquote(expression((.(ticks[i]))^2))))     
  axis(4, at=ticks^2, labels=axlabel)
  box()
  result <- list("freq"=x$freq,
                 "spectrum"=cbind("0.5 %"=quant[,1],
                                  "median"=quant[,2],
                                  "99.5 %"=quant[,3],
                                  "mean"=expect))
  invisible(result)
}



ppsample.bspec <- function(x, start=x$start, ...)
# sample from posterior predictive distribution
{
  kappa <- c(0, rep(1,(x$N %/% 2) -1), ifelse(((x$N %% 2) == 0),0,1))
  # sample spectrum:
  spectrumsample <- sample(x, size=1, two.sided=TRUE)
  # sample data (in Fourier domain):
  Acoefsample <- rnorm(n=length(spectrumsample), mean=0, sd=sqrt(spectrumsample/(1+kappa)))
  Bcoefsample <- rnorm(n=length(spectrumsample), mean=0, sd=sqrt(spectrumsample/(1+kappa)))
  real <- sqrt(x$N/x$deltat) * Acoefsample
  imag <- -sqrt(x$N/x$deltat) * Bcoefsample
  real <- c(real, rev(real[kappa>0]))
  imag <- c(imag, -rev(imag[kappa>0]))
  noiseFT <- complex(x$N, real=real, imaginary=imag)
  # transform to time domain:
  noiseTS <- Re(fft(noiseFT, inverse=TRUE)) / x$N
  noiseTS <- ts(noiseTS, start=0, deltat=x$deltat)
  return(noiseTS)
}



acf.bspec <- function(x, spec=NULL,
                      type=c("covariance", "correlation"),
                      two.sided=x$two.sided, ...)
{
  if (is.null(spec)){
    spec <- expectation(x, two.sided=FALSE)
    vars <- variance(x, two.sided=FALSE)
    if (all(is.finite(vars)))
      estimate.errors <- TRUE
    else {
      stderr <- rep(Inf, x$N%/%2 + 1)
      estimate.errors <- FALSE
    }
  }
  else {
    stderr <- rep(0, x$N%/%2 + 1)
    estimate.errors <- FALSE
    if (two.sided){
      if ((x$N%%2)==0)
        multi <- c(1, rep(2,(x$N %/% 2)-1), 1)
      else
        multi <- c(1, rep(2,(x$N %/% 2)-1))
      spec <- multi * spec
    }        
  }
  type <- match.arg(type)
  lags <- (0:(x$N%/%2)) * x$deltat
  kappa <- c(0, rep(1,(x$N %/% 2) -1), ifelse(((x$N %% 2) == 0),0,1))
  if (x$freq[1] != 0) kappa <- kappa[-1]
  if (all(is.finite(spec))){
    autocov1 <- function(deltat)
    {
      expect <- sum(spec*(1+kappa)/2*cos(2*pi*x$freq * deltat)) / (x$N*x$deltat)
      return(c(expect, NA))
    }
    autocov2 <- function(deltat)
    {
      coef <- (1+kappa) * 0.5 * cos(2*pi*x$freq * deltat) / (x$N*x$deltat)
      expect <- sum(spec*coef) 
      stdev <- sqrt(sum(vars*coef^2))
      return(c(expect,stdev))
    }
    if (estimate.errors){
      dummy <- apply(matrix(lags,ncol=1),1,autocov2)
      acf <- dummy[1,]
      stderr <- dummy[2,]
    }
    else {
      dummy <- apply(matrix(lags,ncol=1),1,autocov1)
      acf <- dummy[1,]
    }
  }
  else acf <- rep(Inf, length(lags))
  result <- list(lag = lags,
                 acf = acf,
                 stderr = stderr,
                 type = type,
                 N = x$N,
                 bspec = deparse(substitute(x)))
  if (type=="correlation")
    result$acf <- result$acf/result$acf[1]
  class(result) <- "bspecACF"
  return(result)
}



expectation.bspecACF <- function(x, ...)
{
  if (x$type == "covariance")
    result <- x$acf
  else {
    result <- rep(NA, length(x$acf))
    warning("'expectation()' only defined for autocovariances, not for autocorrelations")
  }
  return(result)
}



variance.bspecACF <- function(x, ...)
{
  if (x$type == "covariance")
    result <- x$stderr^2
  else {
    result <- rep(NA, length(x$acf))
    warning("'variance()' only defined for autocovariances, not for autocorrelations")
  }
  return(result)
}



plot.bspecACF <- function(x, ci=0.95,
                          type = "h", xlab = NULL, ylab = NULL,
                          main = NULL, ci.col = "blue", ...)
{
  if (is.null(xlab)) xlab <- "lag"
  if (is.null(ylab)) ylab <- ifelse(x$type=="correlation", "autocorrelation", "autocovariance")
  plot(x$lag, x$acf, type=type,
       xlab=xlab, ylab=ylab, ...)
  abline(h=0)
  if ((x$type=="covariance") && all(is.finite(x$stderr)) && (any(x$stderr>0))) {
    z <- qnorm(1-(1-ci)/2)
    cilines <- cbind(x$acf-z*x$stderr, x$acf+z*x$stderr)
    matlines(x$lag, cilines, col=ci.col, lty="dashed")
    
    z <- 1/sqrt(1-ci)
    cilines <- cbind(x$acf-z*x$stderr, x$acf+z*x$stderr)
    matlines(x$lag, cilines, col="red", lty="dashed")
  }
  invisible(x)
}



temper.bspec <- function(x, temperature=2.0, likelihood.only=TRUE, ...)
{
  tempered <- list(freq = x$freq,
                   scale = NA,
                   df = NA,
                   priorscale = x$priorscale,
                   priordf = x$priordf,
                   datassq = x$datassq,
                   datadf = x$datadf,
                   N = x$N,
                   deltat = x$deltat,
                   deltaf = x$deltaf,
                   start = x$start,
                   call = x$call,
                   two.sided = x$two.sided)
  if (temperature != 1) {
    tempered <- c(tempered, "temperature"=temperature)
    if (likelihood.only){
      tempered$scale <- (x$priordf*x$priorscale + (x$datadf/temperature)*x$datassq) / (x$priordf + x$datadf/temperature)
      tempered$df    <- x$priordf + (x$datadf)/temperature
    }
    else {
      stopifnot(temperature < (min(x$df)+2)/2)
      # ...for T < (nu+2)/2...
      tempered$scale <- (x$df*x$scale) / (2 + x$df - 2*temperature)
      tempered$df    <- (x$df+2)/temperature - 2
    }
  }
  else{
    tempered$scale <- (x$priordf*x$priorscale + x$datassq) / (x$priordf + x$datadf)
    tempered$df    <- x$datadf + x$priordf
  }
  class(tempered) <- "bspec"
  return(tempered)
}



temperature.bspec <- function(x, ...)
{
  if (is.element("temperature", names(x)))
    result <- x$temperature
  else
    result <- 1.0
  return(result)
}



dprior.bspec <- function(x, theta,
                         two.sided=x$two.sided, log=FALSE, ...)
# (log-) prior density
# `theta' is the 2-sided PSD (S_2(f_j))
{
  stopifnot(length(theta) == length(x$freq))
  if (any(theta<=0))
    prior <- -Inf
  else{
    if (two.sided) {
      Neven <- ((x$N %% 2) == 0)
      index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
      if (x$freq[1] != 0) index <- index-1
      theta[index] <- theta[index] * 2
    } 
    halfdf <- x$priordf/2
    prior <- sum(halfdf*log(halfdf) - lgamma(halfdf) + halfdf*log(x$priorscale) - (halfdf+1)*log(theta) - (x$priordf*x$priorscale)/(2*theta))
  }
  if (!log)
    prior <- exp(prior)
  return(prior)
}



likelihood.bspec <- function(x, theta,
                             two.sided=x$two.sided, log=FALSE, ...)
# (log-) likelihood
# `theta' is the 2-sided PSD (S_2(f_j))
{
  stopifnot(length(theta) == length(x$freq))
  if (any(theta<=0))
    likeli <- -Inf
  else {
    if (two.sided) {
      Neven <- ((x$N %% 2) == 0)
      index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
      if (x$freq[1] != 0) index <- index-1
      theta[index] <- theta[index] * 2
    } 
    likeli <- -0.5*sum(x$datadf)*log(2*pi) + sum(-(x$datadf/2)*log(theta) - x$datassq/(2*theta))
  }
  if (! log)
    likeli <- exp(likeli)
  return(likeli)
}



marglikelihood.bspec <- function(x, log=FALSE, ...)
# (log-) likelihood
{
  #likeli <- -0.5*sum(x$datadf)*log(2*pi) + sum(-(x$datadf/2)*log(theta) - x$datassq/(2*theta))
  likeli <- -sum(((x$priordf+x$datadf)/2) * log(1 + x$datassq/(x$priordf*x$scale)))
  likeli <- likeli + sum(lgamma((x$priordf+x$datadf)/2) - lgamma(x$priordf/2) - (x$datadf/2)*(log(x$priordf)+log(pi)) - (x$datadf/2)*log(x$scale))
  if (! log)
    likeli <- exp(likeli)
  return(likeli)
}



dposterior.bspec <- function(x, theta,
                             two.sided=x$two.sided, log=FALSE, ...)
# (log-) posterior density
# `theta' is the 2-sided PSD (S_2(f_j))
{
  stopifnot(length(theta) == length(x$freq))
  if (any(theta<=0))
    posterior <- -Inf
  else {
    if (two.sided) {
      Neven <- ((x$N %% 2) == 0)
      index <- 2:ifelse(Neven, x$N %/% 2, (x$N %/% 2) + 1)
      if (x$freq[1] != 0) index <- index-1
      theta[index] <- theta[index] * 2
    } 
    halfdf <- x$df/2
    posterior <- sum(halfdf*log(halfdf) - lgamma(halfdf) + halfdf*log(x$scale) - (halfdf+1)*log(theta) - (x$df*x$scale)/(2*theta))
  }
  if (!log)
    posterior <- exp(posterior)
  return(posterior)
}



matchedfilter <- function(data, signal,
                          noisePSD,
                          timerange=NA,
                          reconstruct=TRUE,
                          two.sided=FALSE)
# This is algorithm II (Appendix A3e) from the Student-t filtering paper.
{
  # check 'data' properties:
  stopifnot(is.ts(data) | is.vector(data))
  if (!is.ts(data)) {
    data <- as.ts(data)
    warning("argument 'data' is not a time-series object; default conversion 'as.ts(data)' applied.")
  }
  N <- length(data)
  FTlength <- (N %/% 2) + 1   # size of (nonredundant) FT output
  Neven <- ((N %% 2) == 0)    # indicator for even N
  deltat <- 1 / tsp(data)[3]
  deltaf <- 1 / (N*deltat)
  freq <- (0:(FTlength-1)) * deltaf
  kappa <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2))
  nonredundant <- 1:FTlength
  dataStart <- tsp(data)[1]
  dataFT <- fft(as.vector(data))
  if (any(is.na(timerange)))
    timerange <- dataStart + c(0,N*deltat) + c(1,-1) * 0.1 * (N*deltat)
  stopifnot(diff(timerange)>0, any((time(data)>=timerange[1]) & (time(data)<=timerange[2])))
  
  # check 'signal' properties:
  stopifnot(((is.ts(signal)) && (tsp(signal)[3]==tsp(data)[3]))
            | (is.vector(signal) && (length(signal)<=N))
            | (is.matrix(signal) && (nrow(signal)<=N)))
  if (is.ts(signal)) {
    signalStart <- tsp(signal)[1]
    if (is.mts(signal)) {
      dn <- dimnames(signal)
      k <- ncol(signal)
      attributes(signal) <- NULL
      signal <- matrix(signal, ncol=k, dimnames=dn)
    }
    else {
      k <- 1
      signal <- matrix(signal, ncol=k)
    }
  }
  else {
    signalStart <- 0.0
    if (is.vector(signal)) {
      k <- 1
      signal <- matrix(signal, ncol=k)
    }
    else k <- ncol(signal)
  }

  if ((k > 1) && !isTRUE(all.equal(as.vector(cor(signal) - diag(k)), rep(0,k^2))))
    warning("Please check: signal components (columns of 'signal' argument) need to be orthogonal.")
  
  stopifnot(signalStart > -N*deltat, signalStart <= 0)
  # do zero-padding if necessary:
  if (nrow(signal) < N)
    signal <- rbind(signal, matrix(0, ncol=ncol(signal), nrow=N-nrow(signal)))

  # check 'noisePSD' properties:
  stopifnot(is.function(noisePSD) || (is.vector(noisePSD) && (length(noisePSD)==FTlength)))
  if (is.function(noisePSD))
    S1 <- noisePSD(freq)
  else 
    S1 <- noisePSD
  if (two.sided)
    S1 <- kappa * S1
  stopifnot(all(S1 > 0.0)) # (?)
  S1ext <- c(S1, rev(S1[kappa==2]))

  # Fourier-transform signal basis waveforms:
  signalFT <- mvfft(signal)
  # (result is an (N x k) matrix).

  # compute signal bases' norms:
  normVec <- (deltat/N) * apply(signalFT, 2, function(x)sum(abs(x)^2/S1ext))
  
  # convolutions of signal waveforms with data:
  convMat <- matrix(NA, nrow=N, ncol=k)
  for (i in 1:k)
    convMat[,i] <- dataFT * Conj(signalFT[,i]) / S1ext

  # apply inverse FT:
  FTMat <- Re(mvfft(convMat, inverse=TRUE))

  # compute logarithmic profile likelihood:
  maxLLR <- (deltat/N)^2 * (FTMat^2 %*% (1/normVec))

  # vector of time points corresponding to filter output:
  timevec <- dataStart-signalStart + (0:(N-1))*deltat
  first <- timevec >= dataStart + N*deltat
  timevec[first] <- timevec[first] - N*deltat
  
  # rearrange output to match input data:
  oldindex <- 1:N
  oldindex <- c(oldindex[first], oldindex[!first])
  timevec <- timevec[oldindex]
  maxLLR <- maxLLR[oldindex]
  inRange <- (timevec>=timerange[1] & timevec<=timerange[2])
  
  imax <- which.max(maxLLR * inRange)
  tHat <- timevec[imax]

  if (reconstruct) {  # determine best-fitting signal:
    betaHat <- (deltat/N) * as.vector(FTMat[oldindex[imax],]) / normVec
    signalRec <- as.vector(signalFT %*% betaHat)
    signalRec <- signalRec * exp(-2*pi*complex(real=0,imaginary=1)*c(freq,rev(-freq[kappa==2])) * (tHat-dataStart+signalStart))
    signalRec <- Re(fft(signalRec, inverse=TRUE))/N
    signalRec <- ts(signalRec, start=dataStart, deltat=deltat)
  }
  else {
    betaHat <- rep(NA,k)
    signalRec <- NULL
  }
  
  return(list("maxLLR"=maxLLR[imax],
              "maxLLRseries"=ts(maxLLR, start=timevec[1], deltat=deltat),
              "timerange"=timerange,
              "betaHat"=betaHat,
              "tHat"=tHat,
              "reconstruction"=signalRec,
              "call" = match.call(expand.dots=FALSE)))
}



studenttfilter <- function(data, signal,
                           noisePSD, df=10,
                           timerange=NA,
                           deltamax = 1e-6,
                           itermax = 100,
                           reconstruct=TRUE,
                           two.sided=FALSE)
# This is algorithm IV (Appendix A3e) from the Student-t filtering paper.
{
  # check 'data' properties:
  stopifnot(is.ts(data) | is.vector(data))
  if (!is.ts(data)) {
    data <- as.ts(data)
    warning("argument 'data' is not a time-series object; default conversion 'as.ts(data)' applied.")
  }
  N <- length(data)
  FTlength <- (N %/% 2) + 1   # size of (nonredundant) FT output
  Neven <- ((N %% 2) == 0)    # indicator for even N
  deltat <- 1 / tsp(data)[3]
  deltaf <- 1 / (N*deltat)
  freq <- (0:(FTlength-1)) * deltaf
  kappa <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2))
  nonredundant <- 1:FTlength
  dataStart <- tsp(data)[1]
  dataFT <- fft(as.vector(data))
  if (any(is.na(timerange)))
    timerange <- dataStart + c(0,N*deltat) + c(1,-1) * 0.1 * (N*deltat)
  stopifnot(diff(timerange)>0,
            any((time(data)>=timerange[1]) & (time(data)<=timerange[2])))
  
  # check 'signal' properties:
  stopifnot(((is.ts(signal)) && (tsp(signal)[3]==tsp(data)[3]))
            | (is.vector(signal) && (length(signal)<=N))
            | (is.matrix(signal) && (nrow(signal)<=N)))
  if (is.ts(signal)) {
    signalStart <- tsp(signal)[1]
    if (is.mts(signal)) {
      dn <- dimnames(signal)
      k <- ncol(signal)
      attributes(signal) <- NULL
      signal <- matrix(signal, ncol=k, dimnames=dn)
    }
    else {
      k <- 1
      signal <- matrix(signal, ncol=k)
    }
  }
  else {
    signalStart <- 0.0
    if (is.vector(signal)) {
      k <- 1
      signal <- matrix(signal, ncol=k)
    }
    else k <- ncol(signal)
  }
  stopifnot(signalStart > -N*deltat, signalStart <= 0)

  if ((k > 1) && !isTRUE(all.equal(as.vector(cor(signal) - diag(k)), rep(0,k^2))))
    warning(paste("Please check: signal components (columns of 'signal' argument) need to be orthogonal."))
  
  # do zero-padding if necessary:
  if (nrow(signal) < N)
    signal <- rbind(signal, matrix(0, ncol=ncol(signal), nrow=N-nrow(signal)))

  # check 'noisePSD' properties:
  stopifnot(is.function(noisePSD) || (is.vector(noisePSD) && (length(noisePSD)==FTlength)))
  if (is.function(noisePSD))
    S1 <- noisePSD(freq)
  else 
    S1 <- noisePSD
  if (two.sided)
    S1 <- kappa * S1
  stopifnot(all(S1 > 0.0))
  S1ext <- c(S1, rev(S1[kappa==2])) # ('extended' PSD vector)

  # check 'df' properties:
  stopifnot(is.function(df) || (is.vector(df) && (length(df)==FTlength | length(df)==1)))
  if (is.function(df))
    nu <- df(freq)
  else {
    if (length(df)==1)
      nu <- rep(df, FTlength)
    else
      nu <- df
  }
  stopifnot(all(nu > 0.0))

  stopifnot(itermax > 1, deltamax > 0)
  
  # Fourier-transform signal basis waveforms:
  signalFT <- mvfft(signal)
  # (result is an (N x k) matrix).
  
  # vector of time points / -shifts corresponding to filter output:
  timevec <- dataStart-signalStart + (0:(N-1))*deltat
  first <- timevec >= dataStart + N*deltat
  timevec[first] <- timevec[first] - N*deltat  
  inRange <- (timevec>=timerange[1] & timevec<=timerange[2])
  
  # compute (unnormalized) log-likelihood under "noise-only" model (H_0):
  LL0 <- -sum(((nu+2)/2) * log(1 + (1/nu) * abs(dataFT[1:FTlength])^2 / ((N/(4*deltat)) * S1)))

  EMprogress <- NULL
  deltaLLR <- 1
  LLRprev <- 0
  S1adapted <- S1ext
  EMiter <- 1
  while ((deltaLLR > deltamax) & (EMiter <= itermax)) { # EM-iterations:
    # compute signal bases' norms:
    normVec <- (deltat/N) * apply(signalFT, 2, function(x)sum(abs(x)^2/S1adapted))

    # convolutions of signal waveforms with data:
    convMat <- matrix(NA, nrow=N, ncol=k)
    for (i in 1:k)
      convMat[,i] <- dataFT * Conj(signalFT[,i]) / S1adapted

    # apply inverse FT:
    FTMat <- Re(mvfft(convMat, inverse=TRUE))

    # compute logarithmic profile likelihood:
    maxLLR <- (deltat/N)^2 * (FTMat^2 %*% (1/normVec))

    # determine ML (subject to time constraint):
    imax <- which.max(maxLLR * inRange)
    tHat <- timevec[imax]

    # determine best-fitting signal:
    betaHat <- (deltat/N) * as.vector(FTMat[imax,]) / normVec
    signalRec <- as.vector(signalFT %*% betaHat)
  
    # determine (conditional) noise residuals:
    dataFTshifted <- dataFT * exp(-2*pi*complex(real=0,imaginary=1)*c(freq,rev(-freq[kappa==2])) * (-1) * (tHat-dataStart+signalStart))
    residual2 <- abs(dataFTshifted[1:FTlength] - signalRec[1:FTlength])^2

    # compute log-likelihood under "signal" model (H_1):
    LL1 <- -sum(((nu+2)/2) * log(1 + (1/nu) * residual2 / ((N/(4*deltat)) * S1)))
    LLR <- LL1 - LL0
    deltaLLR <- LLR - LLRprev
    LLRprev <- LLR

    # keep track of EM-algorithm's progress:
    EMprogress <- rbind(EMprogress,
                        c("iteration"=EMiter, "LLR"=LLR))
    
    # adapt the PSD:
    S1adapted <- (nu/(nu+2)) * S1 + (2/(nu+2)) * (2*deltat/N) * residual2
    S1adapted <- c(S1adapted, rev(S1adapted[kappa==2]))

    EMiter <- EMiter + 1
  }

  if (reconstruct) {
    # time-shift the reconstructed signal, convert to time-domain:
    signalRec <- signalRec * exp(-2*pi*complex(real=0,imaginary=1)*c(freq,rev(-freq[kappa==2])) * (tHat-dataStart+signalStart))
    signalRec <- Re(fft(signalRec, inverse=TRUE))/N
    signalRec <- ts(signalRec, start=dataStart, deltat=deltat)
  }
  else
    signalRec <- NULL
  
  return(list("maxLLR"=LLR,
              "timerange"=timerange,
              "betaHat"=betaHat,
              "tHat"=tHat,
              "reconstruction"=signalRec,
              "EMprogress"=EMprogress,
              "call" = match.call(expand.dots=FALSE)))
}



empiricalSpectrum <- function(x, two.sided=FALSE)
# returns the (by default one-sided) "empirical" spectrum of a time series
# (spectral power based on Fourier transform)
{
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("argument 'x' is not a time-series object, default conversion 'as.ts(x)' applied.")
  }
  N <- length(x)
  FTlength <- (N %/% 2) + 1   # size of (nonredundant) FT output
  Neven <- ((N %% 2) == 0)    # indicator for even N
  deltat <- 1 / tsp(x)[3]
  deltaf <- 1 / (N*deltat)
  freq <- (0:(FTlength-1)) * deltaf
  kappa <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2))
  # Fourier transform:
  xFFT <- fft(as.vector(x))
  # (`fft()' yields unnormalised FT)
  nonredundant <- 1:FTlength
  spec <- kappa * (deltat/N) * abs(xFFT[nonredundant])^2
  if (two.sided) spec <- spec / kappa
  result <- list("frequency"=freq,
                 "power"=spec,
                 "kappa"=kappa,
                 "two.sided"=two.sided)
  return(result)
}



snr <- function(x, psd, two.sided=FALSE)
# signal-to-noise ratio (SNR) 
{
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("argument 'x' is not a time-series object, default conversion 'as.ts(x)' applied.")
  }
  N <- length(x)
  FTlength <- (N %/% 2) + 1   # size of (nonredundant) FT output
  stopifnot(is.function(psd) || (is.vector(psd) && (length(psd)==FTlength)))
  Neven <- ((N %% 2) == 0)    # indicator for even N
  deltaT <- 1 / tsp(x)[3]
  deltaF <- 1 / (N*deltaT)
  freq <- (0:(FTlength-1)) * deltaF
  kappa <- c(1, rep(2,FTlength-2), ifelse(Neven,1,2))
  nonredundant <- 1:FTlength
  # determine spectral density:
  if (is.function(psd)) spec <- psd(freq)
  else spec <- psd
  # re-scale to 1-sided (if given as 2-sided):
  if (two.sided) spec <- kappa * spec
  # determine standard deviations in Gaussian (Whittle) model:
  sigma <- sqrt( (N/deltaT) * (spec/(kappa^2)) )
  # do signal FT:
  xFT <- fft(as.vector(x))
  rho <- sqrt(sum((Re(xFT[nonredundant])/sigma)^2) + sum((Im(xFT[nonredundant])/sigma)^2))
  return(rho)
}



welchPSD <- function(x, seglength,
                     two.sided=FALSE,
                     windowfun=tukeywindow,
                     method=c("mean","median"),
                     windowingPsdCorrection=TRUE,...)
# Spectrum estimation by averaging over overlapping, windowed data segments
# following
#   Welch, P. D. (1967): "The use of Fast Fourier Transform
#   for the estimation of power spectra..."
#   http://dx.doi.org/10.1109/TAU.1967.1161901
# (always using 50% overlap for now)
#
# Arguments:
#   x         :  a time series object
#   seglength :  the length of overlapping segments (in time units)
#   two.sided :  logical flag indicating whether 2-sided (insted of 1-sided) is desired
#   windowfun :  a windowing function
#   windowingPsdCorrection: flag indicating whether to correct (re-scale) PSD estimate for windowing effect
#   ...       :  additional arguments to the windowing function
{
  if (!is.ts(x)) {
    x <- as.ts(x)
    warning("argument 'x' is not a time-series object, default conversion 'as.ts(x)' applied.")
  }
  N <- length(x)       # total number of samples
  tspx <- tsp(x)       # properties of time series x
  deltaT <- deltat(x)  # sampling cadence (=1/rate)
  if (seglength > diff(tspx[1:2])) {
    warning("length of 'x' is less than requested segment length!")
    return()
  }
  L <- floor(seglength/deltaT)  # number of samples in each segment
  if (L %% 2 != 0) L <- L-1     # make sure L is even
  if (L < 4) {
    warning("requested segment length leaves less than 4 samples per segment!")
    return()
  }
  win <- windowfun(L, ...)
  if (windowingPsdCorrection) {
    winRms <- sqrt(mean(win^2)) # windowing root-mean-square
    win <- win / winRms
  }
  K <- N %/% (L/2) - 1 # number of overlapping segments
  specMat <- matrix(NA, nrow=L/2+1, ncol=K)
  for (i in 1:K) {     # loop over segments:
    segment <- ts(x[((i-1)*(L/2)+1):((i-1)*(L/2)+L)],
                  start=0, deltat=deltaT)
    espec <- empiricalSpectrum(segment*win, two.sided=two.sided)
    specMat[,i] <- espec$power
  }
  # each column of 'specMat' now holds an "empirical" spectrum
  # of one of the K segments.
  # Compute average:
  method <- match.arg(method)
  if (method=="mean")
    spec <- rowMeans(specMat)
  else if (method=="median") {
    medianCorrection <- rep(2.0 / qchisq(0.5,df=2), L/2+1)
    medianCorrection[espec$kappa==1] <- 1.0 / qchisq(0.5,df=1)
    spec <- apply(specMat,1,median) * medianCorrection
  }
  
  result <- list("frequency"=espec$frequency,
                 "power"=spec,
                 "kappa"=espec$kappa,
                 "two.sided"=two.sided,
                 "segments"=K)
}



tukeywindow <- function(N, r=0.1)
# a.k.a. "cosine tapered" or "split cosine bell window."
# (r is the fraction of the window
#  in which it behaves sinusoidally;
#  (1-r) is the `flat' fraction.)
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0, r>=0, r<=1)
  r2 <- r/2
  result <- rep(1, N)
  i <- (0:(N-1))
  if (r>0) {
    leftTail  <- (i <= r2*N)
    rightTail <- (i >= (1-r2)*N)
    result[leftTail]  <- 0.5*(1-cos(pi*(i[leftTail]/(r2*N))))
    result[rightTail] <-  0.5*(1-cos(pi*((N-i[rightTail])/(r2*N))))
  }
  return(result)
}



squarewindow <- function(N)
# a.k.a. "rectangular" or "boxcar window".
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0)
  return(rep(1.0, N))
}



hannwindow <- function(N)
# cosine shape; zero derivative at both ends.
# same as cosinewindow(N, a=2).
# same as tukeywindow(N, r=1).
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0)
  i <- (0:(N-1))
  result <- 0.5 * (1-cos((i/N) * 2*pi))
  return(result)
}



welchwindow <- function(N)
# Parabola shape.
# a.k.a. "Riesz", "Bochner" or "Parzen window".
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0)
  i <- (0:(N-1))
  result <- 1 - ((i-N/2)/(N/2))^2
  return(result)
}



trianglewindow <- function(N)
# triangular shape
# a.k.a. "Bartlett" or "Fejer window".
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0)
  i <- (0:(N-1))
  result <- 1 - abs((i-N/2)/(N/2))
  return(result)
}



hammingwindow <- function(N, alpha=0.543478261)
# non-zero offset at ends, but zero derivative.
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0, alpha>=0, alpha<=1)
  i <- (0:(N-1))
  result <- alpha - (1-alpha) * cos((i/N) * 2*pi)
  return(result)
}



cosinewindow <- function(N, alpha=1)
# sine shape; non-zero derivative at both ends (for alpha=1.0).
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0, alpha>0)
  i <- (0:(N-1))
  result <- sin((i/N) * pi)^alpha
  return(result)
}



kaiserwindow <- function(N, alpha=3)
# a.k.a. "Kaiser-Bessel window".
{
  stopifnot(length(N)==1, N==round(N), is.finite(N), N>0, is.finite(alpha), alpha>0)
  i <- (0:(N-1))
  result <- besselI(pi*alpha*sqrt(1-(2*i/N-1)^2), nu=0) / besselI(pi*alpha, nu=0)
  return(result)
}



is.bspec <- function(object)
{
  is.element("bspec", class(object))
}



is.bspecACF <- function(object)
{
  is.element("bspecACF", class(object))
}



one.sided.bspec <- function(x, ...)
{
  x$two.sided <- FALSE
  return(x)
}



two.sided.bspec <- function(x, ...)
{
  x$two.sided <- TRUE
  return(x)
}


##################################

bspec <- function(x, ...)
{
  UseMethod("bspec")
}


expectation <- function(x, ...)
{
  UseMethod("expectation")
}


variance <- function(x, ...)
{
  UseMethod("variance")
}


one.sided <- function(x, ...)
{
  UseMethod("one.sided")
}


two.sided <- function(x, ...)
{
  UseMethod("two.sided")
}


sample <- function(x, ...)
{
  UseMethod("sample")
}


sample.default <- base::sample
formals(sample.default) <- c(formals(sample.default), alist(... = ))


ppsample <- function(x, ...)
{
  UseMethod("ppsample")
}


acf <- function(x, ...)
{
  UseMethod("acf")
}

acf.default <- stats::acf


temper <- function(x, ...)
{
  UseMethod("temper")
}


temperature <- function(x, ...)
{
  UseMethod("temperature")
}


dprior <- function(x, ...)
{
  UseMethod("dprior")
}


dposterior <- function(x, ...)
{
  UseMethod("dposterior")
}


likelihood <- function(x, ...)
{
  UseMethod("likelihood")
}


marglikelihood <- function(x, ...)
{
  UseMethod("marglikelihood")
}
