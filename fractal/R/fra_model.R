#######################################################################################
# FRACTAL funtionality for long memory models
#
#   lmACF
#   lmConfidence
#   lmConvert
#   lmModel
#   lmSDF
#   lmSimulate
#
#######################################################################################

##
# lmACF
##

"lmACF" <- function(x, lag.max=32, type="correlation")
{
  if (!is(x,"lmModel"))
    stop("x must be an object of class \"lmModel\"")
  if (x$model == "dfbm")
    stop("FBM model not supported in ACF estimator")

  checkScalarType(type,"character")
  checkScalarType(lag.max,"integer")

  type    <- match.arg(lowerCase(type),c("correlation","covariance","partial"))
  partial <- is.element(type,"partial")

  if (x$model == "ppl"){

    # define local functions
    "Bm" <- function(m, alpha){
	  q <- 4 * m + 1
	  (pi^(alpha + q) * ((q + 1) * q * (alpha + q + 2) -
	    (alpha + q) * pi^2))/((alpha + q) * (alpha + q + 2) * gamma(q + 2))
    }

    "DeltaTau" <- function(tau, alpha, dterms){
      x0 <- pi * tau + pi/2
      result <- 0

      # compute values for n=1 case to get going
      an <- alpha * x0^(alpha - 1)
      An <- 2

      for(n in seq(from=1, by=2, length=dterms)){

        result <- result + an * An
        An <- 2 * (n + 2) * (pi/2)^(n + 1) - (n + 2) * (n + 1) * An
          an <- (an * (alpha - n) * (alpha - n - 1))/((n + 2) * (n + 1) * x0^2)
      }

      result * (-1)^(tau + 1)
    }

    temp    <- 0:lag.max
    temp[1] <- x$variance

    if (lag.max > 0){

      IandDeltaTau    <- 1:lag.max
      IandDeltaTau[1] <- sum(sapply(x$bterms:0, Bm, x$alpha))

      if(lag.max > 1){
        lags <- 2:lag.max
        IandDeltaTau[lags] <- sapply(lags - 1, DeltaTau, x$alpha, x$dterms)
      }

      temp[2:(lag.max + 1)] <- (2 * x$Cs * cumsum(IandDeltaTau))/(2 * pi * (1:lag.max))^(x$alpha + 1)
    }

    exponent <- x$alpha
  }
  else if (x$model == "fdp"){

  	temp    <- seq(0,lag.max)
    temp[1] <- x$variance

    if(lag.max){
      positive.lags <- seq(lag.max)
      temp[2:(lag.max + 1)] <- (positive.lags + x$delta - 1) / (positive.lags - x$delta)
      temp <- cumprod(temp)
    }

    exponent <- x$delta

  }
  else if (x$model == "fgn"){

    temp <- 0:lag.max
    temp <- (x$variance * (abs(temp + 1)^(2 * x$HG) -
      2 * abs(temp)^(2 * x$HG) + abs(temp - 1)^(2 * x$HG)))/2
    exponent <- x$HG
  }

  if (type == "partial"){
  	data <- ACVStoPACS(temp)
  	tag  <- "PACS"
  }
  else if (type == "correlation"){
  	data <- temp / temp[1]
  	tag  <- "ACF"
  }
  else{
  	data <- temp
  	tag  <- "ACVF"
  }

  wmtsa::create.signalSeries(data,
    position=list(from=ifelse1(partial,1,0), by=1, length=length(data), units="lag"),
    title.data=paste(upperCase(x$model), "(", exponent, ") ", tag, sep=""),
    documentation=paste(tag, " for ", upperCase(x$model), "(", exponent, ") process", sep=""),
    units="ACF")
}

##
# lmConfidence
##

"lmConfidence" <- function(x, model, conf.level=0.95, parm.known=FALSE, n.rep=100000)
{
  if (!is(model,"lmModel"))
    stop("model must be an object of class \"lmModel\"")
  if (model$model == "dfbm")
    stop("FBM model not supported")

  # convert lmModel exponent to alpha
  alpha <- lmConvert(model, to="alpha")

  if (alpha <= -1 || alpha > 0)
    stop("alpha must be on (-1,0],",
      "delta on [0, 0.5), ",
      "or HG on [0.5, 1)")

  x  <- wmtsa::create.signalSeries(x)
  N  <- length(x)
  a  <- N^(1 + alpha)
  sd <- sqrt((variance(as.numeric(x), unbiased=FALSE) * a)/(a - 1))

  q <- if (parm.known || N >= 40000)
    qnorm(1 - (1 - conf.level) / 2)
  else{
    Ns <- c(50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700,
      800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000,
      8000, 9000, 10000, 20000, 30000, 40000)
    if (conf.level == 0.9)
      approx(Ns, c(1.9124, 1.8922, 1.8726, 1.8555, 1.8438,
        1.8287, 1.7686, 1.7433, 1.7238, 1.7122, 1.703,
        1.6987, 1.6911, 1.6915, 1.6863, 1.667, 1.6639,
        1.6603, 1.6542, 1.6558, 1.6543, 1.655, 1.6528,
        1.652, 1.644854, 1.644854, 1.644854), N)$y
    else if (conf.level == 0.95)
      approx(Ns, c(2.476, 2.4335, 2.3909, 2.3623, 2.3362,
        2.3189, 2.2035, 2.1429, 2.1139, 2.0893, 2.073,
        2.063, 2.0508, 2.0506, 2.0395, 2.0113, 2.0013,
        1.9879, 1.9877, 1.9827, 1.9776, 1.9756, 1.9797,
        1.9764, 1.9695, 1.967, 1.959964), N)$y
    else if (conf.level == 0.99)
      approx(Ns, c(3.8463, 3.7429, 3.6501, 3.5771, 3.5152,
        3.4515, 3.1702, 3.035, 2.96, 2.9017, 2.8719,
        2.8462, 2.8192, 2.7948, 2.7804, 2.6976, 2.6684,
        2.6516, 2.643, 2.6275, 2.6219, 2.615, 2.6143,
        2.6112, 2.5933, 2.5885, 2.575829), N)$y
    else{
      ps <- c((1 - conf.level)/2, 1 - (1 - conf.level)/2)
      mean(abs(quantile(rnorm(n.rep) * N^((rnorm(n.rep) *
        sqrt(6))/(pi * sqrt(N))), ps)))
    }
  }

  half.width <- (q * sd) / N^(0.5 + 0.5 * alpha)
  xbar <- mean(x)
  c(xbar - half.width, xbar + half.width)
}

##
# lmConvert
##

"lmConvert" <- function(x, to="delta")
{
  alpha.to.beta <- function(alpha)
  {
    if(alpha <= -1 || alpha >= 0)
      stop(paste(
        "\n-1 < x < 0 required to convert alpha to beta"
        ))
 - alpha - 1
  }
  alpha.to.delta <- function(alpha)
 - alpha/2
  alpha.to.HB <- function(alpha)
  {
    if(alpha <= -3 || alpha >= -1)
      stop(paste(
        "\n-3 < x < -1 required to convert alpha to HB"
        ))
    -0.5 - alpha/2.
  }
  alpha.to.HG <- function(alpha)
  {
    if(alpha <= -1 || alpha >= 1)
      stop(paste(
        "\n-1 < x < 1 required to convert alpha to HG")
        )
    0.5 - alpha/2.
  }
  beta.to.alpha <- function(beta)
  {
    if(beta <= -1 || beta >= 0)
      stop(paste("\n-1 < x < 0 required for beta"))
 - beta - 1
  }
  beta.to.delta <- function(beta)
  {
    if(beta <= -1 || beta >= 0)
      stop(paste("\n-1 < x < 0 required for beta"))
    0.5 * (beta + 1)
  }
  beta.to.HG <- function(beta)
  {
    if(beta <= -1 || beta >= 0)
      stop(paste("\n-1 < x < 0 required for beta"))
    0.5 * beta + 1
  }
  delta.to.alpha <- function(delta)
  -2 * delta
  delta.to.beta <- function(delta)
  {
    if(delta <= 0 || delta >= 0.5)
      stop(paste(
        "\n0.0 < x < 0.5 required to convert delta to beta"
        ))
    2 * delta - 1
  }
  delta.to.HB <- function(delta)
  {
    if(delta <= 0.5 || delta >= 1.5)
      stop(paste(
        "\n0.5 < x < 1.5 required to convert delta to HB"
        ))
    -0.5 + delta
  }
  delta.to.HG <- function(delta)
  {
    if(delta <= -0.5 || delta >= 0.5)
      stop(paste(
        "\n-0.5 < x < 0.5 required to convert delta to HG"
        ))
    0.5 + delta
  }
  HB.to.alpha <- function(HB)
  {
    if(HB <= 0 || HB >= 1)
      stop(paste("\n0 < x < 1 required for HB"))
    -1 - 2 * HB
  }
  HB.to.delta <- function(HB)
  {
    if(HB <= 0 || HB >= 1)
      stop(paste("\n0 < x < 1 required for HB"))
    0.5 + HB
  }
  HG.to.alpha <- function(HG)
  {
    if(HG <= 0 || HG >= 1)
      stop(paste("\n0 < x < 1 required for HG"))
    1 - 2 * HG
  }
  HG.to.beta <- function(HG)
  {
    if(HG <= 0.5 || HG >= 1)
      stop(paste(
        "\n0.5 < x < 1 required to convert HG to beta")
        )
    2 * HG - 2
  }
  HG.to.delta <- function(HG)
  {
    if(HG <= 0 || HG >= 1)
      stop(paste("\n0 < x < 1 required for HG"))
    -0.5 + HG
  }

  if (!is(x,"lmModel"))
    stop("x must be an object of class \"lmModel\"")

  #TODO: what model is beta?
  from <- switch(x$model,
    ppl = "alpha",
    foo = "beta",
    fdp = "delta",
    fgn = "HG",
    dfbm = "HB")


  valid.parm <- c("alpha", "beta", "delta", "HB", "HG")

  to <- match.arg(to, valid.parm)

  from.i <- charmatch(from, valid.parm, nomatch=0)
  if(from.i == 0)
    stop(paste("\n\"from\" argument is", from,
      "but should match one of the following:\n", paste(
      valid.parm, collapse=", ")))
  to.i <- charmatch(to, valid.parm, nomatch=0)
  if(to.i == 0)
    stop(paste("\n\"to\" argument is", to,
      "but should match one of the following:\n", paste(
      valid.parm, collapse=", ")))

  x <- x[[from]]

  switch(from.i,
    switch(to.i,
      x,
      alpha.to.beta(x),
      alpha.to.delta(x),
      alpha.to.HB(x),
      alpha.to.HG(x)),
    switch(to.i,
      beta.to.alpha(x),
      x,
      beta.to.delta(x),
      stop(paste("\ncannot convert beta to HB")),
      beta.to.HG(x)),
    switch(to.i,
      delta.to.alpha(x),
      delta.to.beta(x),
      x,
      delta.to.HB(x),
      delta.to.HG(x)),
    switch(to.i,
      HB.to.alpha(x),
      stop(paste("\ncannot convert HB to beta")),
      HB.to.delta(x),
      x,
      stop(paste("\ncannot convert HB to HG"))),
    switch(to.i,
      HG.to.alpha(x),
      HG.to.beta(x),
      HG.to.delta(x),
      stop(paste("\ncannot convert HG to HB")),
      x))
}

##
# lmModel
##

"lmModel" <- function(model, variance.=1.0, delta=0.45, alpha= -0.9, HG=0.95, HB=0.95,
  innovations.var=NULL, Cs=NULL, bterms=10, dterms=10, M=100)
{
  # For a stationary process, the variance parameter is interpreted as
  # the process variance; for a nonstationary process whose dth order
  # differences are stationary, it is interpreted as the variance of
  # the stationary dth order difference process.  There are two ways
  # of indirectly setting the variance parameter for FD and PPL
  # processes.  The first way is to specify the innovations variance
  # by setting innovations.var (if the user does so, this specification
  # takes precedence over any setting (default or user-specified) of
  # the variance parameter.  The second way (valid only for PPL processes)
  # is to specify the Cs parameter (if the user does so, this takes
  # precedence over any settings of the variance and innovations.var
  # parameters).

  # define local functions
  fdp.var.to.ivar <- function(var,delta) {
   del.s <- if(delta < 0.5) delta else delta - floor(delta+0.5)
   var * (gamma(1-del.s))^2 / gamma(1-2*del.s)}

  fdp.ivar.to.var <- function(ivar,delta) {
   del.s <- if(delta < 0.5) delta else delta - floor(delta+0.5)
   ivar * gamma(1-2*del.s)/(gamma(1-del.s))^2}

# currently not used:
#
#  ppl.ivar.to.var <- function(ivar,alpha){ # function valid for all alpha
#   ppl.Cs.to.ivar(ivar * (4 * exp(1))^alpha,alpha)}

  ppl.Cs.to.var <- function(Cs,alpha){ # function requires alpha > - 5

    # define local functions
    ppl.term.m1.to.m3 <- function(n, alpha){
     tn <- 2*n
     ((-1)^n)*(pi^tn)/(factorial(tn-1)*(alpha+1+tn))
    }

    ppl.term.m3.to.m5 <- function(n){
     tn <- 2*n
     ((-1)^(n+1))*(pi^(tn+2))*(1-2^tn)/(factorial(tn-1)*(alpha+3+tn))
    }

     if(alpha > -1) Cs/((2^alpha) * (alpha + 1))
     else
     {if(alpha==-1) Cs * 6.593110554818031
      else
      {if(alpha > -3){Cs*(1.0+0.5*(sum(ppl.term.m1.to.m3(seq(1,15,2),alpha))+sum(ppl.term.m1.to.m3(seq(2,16,2),alpha))))/((2^(alpha-2))*(alpha+1))}
       else
       {if(alpha==-3) Cs * 185.3064454167528
        else
        {if(alpha > -5) {Cs*(1 - (pi^2/((alpha+2)*(alpha+3))) + (sum(ppl.term.m3.to.m5(seq(1,15,2)))+sum(ppl.term.m3.to.m5(seq(2,16,2))))/(2*(alpha+2)*(alpha+3)))/((2^(alpha-4))*(alpha+1))}
    	else
    	{if(alpha==-5) Cs * 5508.610242714385
         	else stop("sorry - cannot currently handle alpha < -5")}}}}}
  }

  checkScalarType(model,"character")
  checkScalarType(variance.,"numeric")
  if (variance. <= 0.0)
    stop("variance must be positive")
  checkScalarType(delta,"numeric")
  checkScalarType(alpha,"numeric")
  checkScalarType(HG,"numeric")
  checkScalarType(HB,"numeric")
  checkScalarType(bterms,"integer")
  checkScalarType(dterms,"integer")
  checkScalarType(M,"integer")

  model <- match.arg(model,c("ppl","fdp","fgn","dfbm"))

  # check model parameter range
  if (model == "fgn" && (HG >= 1 || HG <= 0))
    stop("HG must be on the interval (0,1)")
  if (model == "dfbm" && (HB >= 1 || HB <= 0))
    stop("HB must be on the interval (0,1)")

  # if the model is ppl, check to see if user has supplied a value
  # for Cs and if Cs is positive; if so, use Cs to set variance.
  # and innovations.var; if not, check to see if use has supplied
  # a value for innovations.var and if innovations.var is positive;
  # if so, use innovations.var to set variance. and Cs; otherwise,
  # use variance. to set innovations.var and Cs
  if (model == "ppl")
	{if(!is.null(Cs))
	 {# user has specified Cs, so check that it is OK
    	  checkScalarType(Cs,"numeric")
          if (Cs <= 0.0) stop("Cs must be positive")
          variance. <- ppl.Cs.to.var(Cs,alpha)
          innovations.var <- Cs/(2*exp(1))^alpha}
         else
         {if(!is.null(innovations.var))
          {# user has innovations.var, so check that it is OK
           checkScalarType(,"numeric")
           if (innovations.var <= 0.0) stop("innovations.var must be positive")
           Cs <- ((2*exp(1))^alpha)*innovations.var
           variance. <- ppl.Cs.to.var(Cs,alpha)}
	  else
          {# neither Cs nor innovations.var has been set, so use variance.
           # to set both of these
           Cs <- variance./ppl.Cs.to.var(1.0,alpha)
           innovations.var <- Cs/(2*exp(1))^alpha}}}

  # if the model is fdp, check to see if user has supplied a value
  # for innovations.var and, if so, use it to set variance if all is OK
  if ((model == "fdp") && !is.null(innovations.var)) {

    # innovations.var has been specified, so check that it is OK
    checkScalarType(innovations.var,"numeric")
    if (innovations.var <= 0.0)
      stop("innovations.var must be positive")

    # since the innovations variance has been supplied, we need to calculate
    # variance rather than use its default or user-supplied value
    variance. <- fdp.ivar.to.var(innovations.var,delta)
  }
  else{
   innovations.var <- fdp.var.to.ivar(variance., delta)
  }

  z <- c(list(model=model),
    switch(model,
      ppl=list(alpha=alpha, Cs=Cs, bterms=bterms, dterms=dterms, innovations.var=innovations.var, variance=variance.),
      fdp=list(delta=delta, innovations.var=innovations.var, variance=variance.),
      fgn=list(HG=HG,M=M, variance=variance.),
      dfbm=list(HB=HB,M=M, variance=variance.)))

  oldClass(z) <- "lmModel"

  z
}

##
# lmSDF
##

"lmSDF" <- function(x, sampling.interval=1, n.freq=NULL,
  n.sample=NULL, with.Nyquist=NULL)
{
  # deparse model and model parameters
  if (!is(x,"lmModel"))
    stop("x must be an object of class \"lmModel\"")
  checkScalarType(sampling.interval, "numeric")

  ######
  #TODO: much of this seems much to complicated for a simple idea
  #ask Don how he wants this to be reformulated

  if (is.null(n.freq)){
    if (!is.null(n.sample))
      n.freq <- n.sample %/% 2
    else
      n.freq <- 32
  }

  if (is.null(with.Nyquist)){
    if (!is.null(n.sample))
      with.Nyquist <- (n.sample %% 2L) == 0L
    else
      with.Nyquist <- TRUE
  }

  delta.f <- 1 / (2 * n.freq + ifelse(with.Nyquist, 0, 1))
  freqs   <- seq(n.freq) * delta.f

  #
  ####### end region to fix

  if (x$model == "ppl"){

    sdf <- x$Cs * sampling.interval * freqs^x$alpha
    exponent <- x$alpha
  }
  else if (x$model == "fdp"){

    sdf <- (x$innovations.var * sampling.interval) / (2 * sin(pi * freqs))^(2 * x$delta)
    exponent <- x$delta
  }
  else if (is.element(x$model, c("fgn","dfbm"))){

    exponent <- ifelse1(x$model == "fgn", x$HG, x$HB)

  	"first.sum" <- function(f, H, M){
      sum(abs(f + ((-M):M))^(-2 * H - 1))
    }

    "second.sum" <- function(f, H, M){
      tH   <- 2 * H
      tHp1 <- tH + 1
      tHp2 <- tH + 2
      tHp3 <- tH + 3
      tHp4 <- tH + 4
      return(1/(tH * (- f + M + 1)^tH) + 1/(tH * (f + M + 1)^tH) - (tHp1 *
		tHp2 * tHp3)/(720 * (- f + M + 1)^tHp4) - (tHp1 * tHp2 *
		tHp3)/(720 * (f + M + 1)^tHp4) + tHp1/(12 * (- f +
		M + 1)^tHp2) + tHp1/(12 * (f + M + 1)^tHp2) + 1/(2 *
		(- f + M + 1)^tHp1) + 1/(2 * (f + M + 1)^tHp1))
    }

    CH <- (gamma(2 * exponent + 1) * sin(pi * exponent))/(2 * pi)^(2 * exponent + 1)

    if (x$model == "fgn"){
      sdf <- 4 * sampling.interval * x$variance * CH * sin(pi * freqs)^2 *
        (sapply(freqs, first.sum, exponent, x$M) + sapply(freqs, second.sum, exponent, x$M))
    }
    else{
      sdf <- sampling.interval * x$variance * CH * (sapply(freqs, first.sum, exponent, x$M) +
        sapply(freqs, second.sum, exponent, x$M))
    }
  }

  by  <- from  <- delta.f / sampling.interval
  tag <- upperCase(x$model)

  if (!as.logical(length(sdf)))
    stop("SDF length is zero. Check model parameters.")

  wmtsa::create.signalSeries(sdf,
    position=list(from=from, by=by, units="frequency"),
    title.data=paste(tag, "(", exponent, ") SDF", sep=""),
    documentation=paste("SDF for ", tag, "(", exponent, ") process", sep=""),
    units="SDF")
}

##
# lmSimulate
##

"lmSimulate" <- function(x, sampling.interval=1,
  mean=0, n.sample=128, generate.Sj=FALSE, Sj=NULL, rn=NULL)
{
  if (is(x,"DaviesHarte"))
    Sj <- x

  if (is.null(Sj)){
    if (!is(x,"lmModel"))
      stop("x must be an object of class \"lmModel\"")

    checkScalarType(n.sample,"integer")
    if (n.sample <= 0)
      stop("number of samples must be a positive number")
    checkScalarType(generate.Sj,"logical")
  }
  else if (!is(x,"DaviesHarte"))
    checkVectorType(Sj,"numeric")

  checkScalarType(mean,"numeric")
  checkScalarType(sampling.interval,"numeric")
  if (sampling.interval <= 0)
      stop("sampling interval must be a positive number")
  if (!is.null(Sj) && generate.Sj)
    return(Sj)

  # define local functions
  "DHMgenSj" <- function(acvf)
  {
  	acvf <- acvf@data
    MM <- length(acvf) - 1
    if(!(MM == nextDyadic(MM)))
      stop("Length of acvf must be 1 greater than a power of 2")
    Sj <- Re(fft(c(acvf, rev(acvf[2:MM]))))[1:(MM + 1)]
    if(!(all(Sj >= 0)))
      stop("some of the S_j's are negative")
    Sj
  }

  "DHMgenSim" <- function(Sj, rn=NULL)
  {
  	if (is.null(rn))
  	  rn <- rnorm(2 * length(Sj) - 2)

    M         <- 2 * length(Sj) - 2
    N         <- M/2
    Yj        <- 0:N
    Yj[1]     <- sqrt(Sj[1]) * rn[1]
    Yj[N + 1] <- sqrt(Sj[N + 1]) * rn[M]
    js        <- 2:N
    Yj[js] <- sqrt(0.5 * Sj[js]) * complex(real=rn[2 * (1:(N - 1))],
      imaginary=rn[2 * (1:(N - 1)) + 1])
    Re(fft(c(Yj, Conj(rev(Yj[js])))))[1:N]/sqrt(M)
  }

  if (is.null(Sj)){

    if (x$model == "ppl"){

      if (x$alpha <= -1)
        stop("alpha must be greater than -1")
      exponent <- x$alpha
      Sj <- DHMgenSj(lmACF(x, nextDyadic(n.sample), type="covariance"))
    }
    else if (x$model == "fdp"){

      "FDsimGenSj" <- function(func, x, n.sample=1024){

        # migrate delta to stationary region
        # NOTE: this doesn't work if delta <= -1, in which case we need
        #       to difference a simulated process a certain number of times
        #       (as a result, the undifferenced process will have to be longer)
         x$delta <- if(x$delta < 0.5) x$delta else x$delta - floor(x$delta+0.5)

        func(lmACF(x, nextDyadic(n.sample), type="covariance"))
      }

  	  Sj <- FDsimGenSj(DHMgenSj, x, n.sample)
  	  exponent <- x$delta
    }
    else if (is.element(x$model, c("fgn","dfbm"))){

      # Force FGN model on ACF since FBM ACF is not supported.
      # In this case, A cumulative summation over the Sj's is
      # performed at the end of this function to migrate the
      # FGN result to FBM
      model <- ifelse1(x$model == "dfbm", lmModel("fgn", HG=x$HB, M=x$M), x)
      Sj <- DHMgenSj(lmACF(model, nextDyadic(n.sample), type="covariance"))

      exponent <- ifelse1(x$model == "fgn", x$HG, x$HB)
    }

    oldClass(Sj) <- "DaviesHarte"
    attr(Sj, "n.sample") <- n.sample
    attr(Sj, "model")    <- x
    attr(Sj, "exponent") <- exponent
  }

  if (generate.Sj)
  	return(Sj)
  else{

  	# unpack Sj attributes
  	Sj.attr  <- attributes(Sj)
  	n.sample <- Sj.attr$n.sample
  	model    <- Sj.attr$model
  	exponent <- Sj.attr$exponent

  	from <- 0
  	tag  <- upperCase(model$model)

  	if (model$model == "fdp"){

      data <- DHMgenSim(Sj, rn=rn)[1:n.sample]

	    if (exponent >= 0.5){

        data <- cumsum(data)
        remaining.order <- floor(exponent + 0.5) - 1

        if (remaining.order == 1)
          data <- cumsum(data)
        else if(remaining.order > 1)
          for(i in 1:remaining.order)
            data <- cumsum(data)
      }
  	}
  	else if (model$model == "dfbm"){
  	  data <- cumsum(DHMgenSim(Sj, rn=rn)[1:n.sample]) + mean
  	}
  	else
  	  data <- DHMgenSim(Sj, rn=rn)[1:n.sample] + mean

    wmtsa::create.signalSeries(data,
      position=list(from=from, by=sampling.interval, units="time"),
      title.data=paste(tag, "(", exponent, ") ", "simulation", sep=""),
      documentation=paste("simulation of ", tag, "(", exponent, ") process", sep=""))
  }
}
