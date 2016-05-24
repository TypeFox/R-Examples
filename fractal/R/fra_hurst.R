################################################
# FRACTAL Hurst coefficient estimators and
# related statistical functions
#
#  Functions:
#
#    dispersion
#    DFA
#    HDEst
#    hurstACVF
#    hurstBlock :: aggAbs, aggVar, diffVar, Higu
#    hurstSpec  :: pgramReg, pgramRegMod, RobInt
#    RoverS
#
################################################

##
# dispersion
##

"dispersion" <- function(x, front=FALSE)
{
  checkScalarType(front,"logical")
  n.sample <- length(x)
  if (n.sample < 32)
    stop("Time series must contain at least 32 points")

  scale <- as.integer(logScale(scale.min=1, scale.max=n.sample, scale.ratio=2))
  sd <- unlist(lapply(scale, function(m, x, n.sample, front){
     nscale   <- n.sample %/% m
     dyad     <- nscale * m
     n.offset <- ifelse1(front, n.sample - dyad, 0)
     stdev(aggregateData(x[seq(1 + n.offset, length=dyad)],
         by=m, FUN=mean))
  },
  x=x,
  n.sample=n.sample,
  front=front))

  list(sd=sd, scale=scale)
}

##
# DFA
##

"DFA" <- function(x, detrend="poly1", sum.order=0, overlap=0,
  scale.max=trunc(length(x)/2), scale.min=NULL, scale.ratio=2,
  verbose=FALSE)
{
  # define local functions
  "polyfit.model" <- function(polyfit.order){
    if (polyfit.order < 0)
      stop("Polynomial fit order must be positive")
    if (polyfit.order == 0)
      return("x ~ 1")
    if (polyfit.order == 1)
      return("x ~ 1 + t")
    else{
      return(paste(c("x ~ 1", "t", paste("I(t^", seq(2, polyfit.order), ")", collapse = " + ", sep = "")),
                        collapse=" + ", sep = ""))
    }
  }

  "regression.poly" <- function(x, model=polyfit.model(1)){

     order <- length(attr(terms(formula(x~1)),"order"))

     t <- x@positions
     x <- x@data

     if (order == 0)
       return(sum((x - mean(x))^2))

     if (order == 1){

       x <- x - mean(x)
       t <- t - mean(t)

       slope <- sum(t*x)/sum(t^2)
       return (sum((x - slope*t)^2))
     }

     fit <- lm(model, data=data.frame(list(t=t, x=x)))
     return(sum(fit$residuals^2))
  }

  "regression.bridge" <- function(x, ...){

     x <- x@data
     N <- length(x)
     bridge <- seq(x[1], x[N], length=N)
     x   <- x - bridge
     return(sum(x^2))
  }

  "regression.none" <- function(x, ...)
     return(sum(x@data^2))

  # check input argument types and lengths
  checkScalarType(detrend,"character")
  checkScalarType(scale.ratio,"numeric")
  checkScalarType(verbose,"logical")
  checkScalarType(overlap,"numeric")

  # obtain series name
  data.name <- deparse(substitute(x))

  # convert sequence to signalSeries class
  x <- wmtsa::create.signalSeries(x)

  # check overlap argument
  if ((overlap < 0) | (overlap >= 1))
   stop("Overlap factor must be in the range [0,1)")

  # map the detrending method
  detrend <- lowerCase(detrend)

  if (substring(detrend,1,4) == "poly"){

	polyfit.order <- as.integer(substring(detrend,5))

        if (is.na(polyfit.order)){
	 stop("Improperly formed detrending string")
        }
        else if (polyfit.order < 0){
	 stop("Polynomial fit order must be positive")
	}

	# form detrending model
    model <- formula(polyfit.model(polyfit.order))

    if (verbose) {
      cat("Detrending model: ");
      print(model)
    }

	# set regression function
    regressor   <- regression.poly
    modstr      <- as.character(model)[2:3]
    regress.str <- paste(modstr,collapse=" ~ ")
  }
  else if (charmatch(detrend, "bridge", nomatch=FALSE)){
    regressor   <- regression.bridge
    regress.str <- "Bridge detrended"
  }
  else if (charmatch(detrend, "none", nomatch=FALSE)){
    regressor  <- regression.none
    regress.str <- "None"
  }
  else stop("Detrending method is not supported")

  # coerce orders to integers
  sum.order <- trunc(sum.order)

  # perform cumulative summation(s) or difference(s)
  if (sum.order > 0){
     for (i in seq(sum.order))
       x <- cumsum(x)
  }
  else if (sum.order < 0){
     for (i in seq(sum.order))
       x <- diff(x)
  }

  # obtain sequence length (after possibly
  # taking differences)
  N <- length(x)

  # create scale for data analysis
  if (is.null(scale.min))
    scale.min <- ifelse1(substring(detrend,1,4) == "poly", 2*(polyfit.order + 1), min(N/4, 4))

  checkScalarType(scale.min, "numeric")
  checkScalarType(scale.max, "numeric")

  # DFA requires integer scales
  scale.min <- trunc(scale.min)
  scale.max <- trunc(scale.max)

  if (scale.min > scale.max)
    stop("Scale minimum cannot exceed scale maximum")
  if (any(c(scale.min, scale.max) < 0))
    stop("Scale minimum and maximum must be positive")
  if (any(c(scale.min, scale.max) > N))
    stop(paste("Scale minimum and maximum must not exceed length",
	"of (differenced/cummulatively summed) time series"))

  # initialize output
  scale <- logScale(scale.min, scale.max, scale.ratio=scale.ratio, coerce=trunc)
  scale <- scale[scale > 1]
  rmse  <- rep(0.0, length(scale))

  # calculate DFA
  for (i in seq(along=scale)){

    sc       <- scale[i]
    noverlap <- trunc(overlap * sc)

    if (verbose)
      cat("Processing scale", sc,"\n")

    cumsum <- 0.0
    count  <- 0

    index <- seq(sc)

    # sum the mse over all blocks
    while(index[sc] <= N){
      rmse[i] <- rmse[i] + regressor(x[index], model=model)
      count   <- count + 1
      index   <- index + sc - noverlap
    }

    # normalize by the number points
    # in all blocks. this normalization
    # factor may differ from the number
    # of points in the time series (N)
    # since N may not be evenly divisible
    # by the current block size. finally,
    # take the sqrt() to form the root
    # mean square error
    rmse[i] <- sqrt(rmse[i] / (count * sc))
  }

  # fit the log-log data
  # weight the coefficients starting from 1 down to 1/N
  xx <- log(scale)
  yy <- log(rmse)
  ww <- 1/seq(along=xx)
  w <- NULL # quell R CMD check scoping issue
  
  logfit   <- lm(y ~ 1 + x, data=data.frame(x=xx, y=yy, w=ww), weights=w)
  exponent <- logfit$coefficients["x"]

  # return the result as a fractal block exponent
  fractalBlock(
    domain       ="time",
    estimator    ="Detrended Fluctuation Analysis",
    exponent     =exponent,
    exponent.name="H",
    scale        =scale,
    stat         =rmse,
    stat.name    ="RMSE",
    detrend      =regress.str,
    overlap      =overlap,
    data.name    =data.name,
    sum.order    =sum.order,
    series       =asVector(x),
    logfit       =logfit)
}

##
# HDEst
##

"HDEst" <- function(NFT, sdf, A=0.3, delta=6/7)
{
  checkScalarType(NFT,"integer")
  checkScalarType(A,"numeric")
  checkScalarType(delta,"numeric")

  # define frequencies
  omega <- ((2 * pi)/NFT) * (1:NFT)

  # notation as in reference: construct X matrix
  L  <- round(A * NFT^delta)
  c1 <- rep(1, L)
  c2 <- log(abs(2 * sin(omega[1:L]/2)))
  c3 <- omega[1:L]^2/2
  X  <- cbind(c1, c2, c3)

  # construct least squares product matrix and find 3rd row:
  TX    <- t(X)
  ProdX <- solve(TX %*% X) %*% TX
  b     <- ProdX[3,]

  # evaluate K.hat, C.hat, and optimum m, as in Reference;
  # use first L non-zero frequencies:
  K.hat <- as.numeric(b %*% log(sdf[1:L]))
  if (K.hat == 0)
    stop("Method fails: K.hat=0")
  C.hat <- ((27/(128 * pi^2))^0.2) * (K.hat^2)^-0.2

  round(C.hat * NFT^0.8)
}

##
# hurstACVF
##

"hurstACVF" <- function(x, Ascale=1000, lag.min=1, lag.max=64)
{
  # check input argument types and lengths
  checkScalarType(Ascale,"numeric")
  checkScalarType(lag.min,"integer")
  checkScalarType(lag.max,"integer")

  # estimates long memory parameters by linear regression of scaled asinh
  # of ACVF vs log(lag) over intermediate lag values.
  x <- wmtsa::create.signalSeries(x)@data
  n.sample <- length(x)

  # evaluate autocovariance function (ACVF)
  acfx <- acf(x, n.sample, "covariance")$acf

  # include lag=0 point
  yy <- c(1, acfx)

  # linear regression over lag.min to lag.max
  fitcoef <- lsfit(log(lag.min:lag.max), asinh(Ascale * yy[(lag.min + 1):
    (lag.max + 1)])/asinh(Ascale))$coef

  # relate slope to long memory parameters (beta=slope of ACVF,
  # alpha=parameter in corresponing PPL model, HG=corresponding
  # Hurst parameter for stationary FGN model.
  beta  <- fitcoef[2] * asinh(Ascale)
  alpha <- -1 - beta
  HG    <- 1 + beta/2

  drop(list(beta=beta, alpha=alpha, HG=HG))
}

##
# hurstBlock
##

"hurstBlock" <- function(x, method="aggAbs", scale.min=8, scale.max=NULL,
  scale.ratio=2, weight=function(x) rep(1,length(x)), fit=lm)
{
  # convert input to signalSeries class
  x <- wmtsa::create.signalSeries(x)

  # initialize variables
  n.sample <- length(x)
  if (is.null(scale.max))
    scale.max <- n.sample

  # obtain series name
  data.name <- deparse(substitute(x))

  # check inputs
  checkScalarType(method, "character")
  method <- match.arg(lowerCase(method), c("aggabs","aggvar","diffvar","higuchi"))
  checkScalarType(scale.ratio, "integer")
  if (scale.ratio < 2)
    stop("scale.ratio must be greater than one")
  checkScalarType(scale.min, "integer")
  if (scale.min < 2)
    stop("scale.min must be greater than one")
  checkScalarType(scale.max, "integer")
  if (scale.max > n.sample)
    stop("scale.max cannot exceed series length")
  if (!is.function(weight))
    stop("weight must be a function")
  if (!is.function(fit))
    stop("fit must be a linear regression function")
  fitstr    <- deparse(substitute(fit))
  supported <- c("lm","lmsreg","ltsreg")
  if (!is.element(fitstr,supported	))
    stop("Supported linear regression functions are:", supported)
  lmfit <- is.element(fitstr,"lm")

  # form the scale vector
  scale <- logScale(scale.min, scale.max, scale.ratio=scale.ratio, coerce=trunc)

  # define block statistic
  statfunc <- ifelse1(method == "aggabs",
    function(x) mean(abs(x)),
    function(x) variance(x, unbiased=FALSE))

  statistic <- unlist(lapply(scale,
    function(blockWidth, x, statfunc, n.sample, higuchi){

      if (higuchi){
        return(mean(abs(aggregateData(x, by=1, FUN=mean, moving=blockWidth))) *
          ((n.sample - 1) / blockWidth))
      }

      nBlocks  <- n.sample %/% blockWidth
      nUpper   <- nBlocks * blockWidth
      xcut     <- as.numeric(x[1:nUpper])

      # reshape recentered time series as a matrix
      tsMat <- matrix(xcut - mean(xcut), blockWidth, nBlocks)

      # evaluate the statistic of each column's/block's mean value
      statfunc(colMeans(tsMat))
    },
    x=x@data, statfunc=statfunc, n.sample=n.sample,
    higuchi=is.element(method,"higuchi")))

  # set up for linear regression of variability vs group size:
  if (is.element(method, c("aggabs","aggvar","higuchi"))){
    xvar <- scale
    yvar <- statistic
  }
  else if (is.element(method,"diffvar")){

    # take first difference, and exclude any points
    # with non-positive differences
    dStat <- -diff(statistic)
    xvar  <- scale[dStat > 0][-length(scale)]
    yvar  <- dStat[dStat > 0]
  }
  w <- NULL # quell R CMD check scoping issue
  
  # fit the log-log data weight the coefficients
  logfit <- ifelse1(lmfit,
    fit(y ~ 1 + x, data=data.frame(x=log(xvar), y=log(yvar), w=weight(yvar)), weights=w),
    fit(y ~ 1 + x, data=data.frame(x=log(xvar), y=log(yvar))))
  exponent <- logfit$coefficients["x"]

  # estimate H from the slope of the regression fit:
  H <- switch(method,
    aggabs =1 + exponent,
    aggvar =1 + exponent/2,
    diffvar=1 + exponent/2,
    higuchi=exponent + 2)

  # return the result as a fractal block exponent
  fractalBlock(
    domain       ="time",
    estimator    ="Hurst coefficient via Aggregated Series",
    exponent     =H,
    exponent.name="H",
    scale        =xvar,
    stat         =yvar,
    stat.name    =method,
    detrend      =NULL,
    overlap      =0,
    data.name    =data.name,
    sum.order    =0,
    series       =asVector(x),
    logfit       =logfit)
}

##
# hurstSpec
##

"hurstSpec" <- function(x, method="standard", freq.max=0.25, dc=FALSE, n.block=NULL,
  weight=function(x) rep(1,length(x)), fit=lm, sdf.method="direct", ...)
{
  # check input arguments
  checkScalarType(method,"character")
  method <- match.arg(lowerCase(method),c("standard","smoothed","robinson"))
  checkScalarType(freq.max,"numeric")
  if (freq.max <= 0.0 || freq.max >= 0.5)
    stop("freq.max must be on the normalized frequency range (0,0.5)")
  if (!is.function(fit))
    stop("fit must be a linear regression function")
  fitstr    <- deparse(substitute(fit))
  supported <- c("lm","lmsreg","ltsreg")
  if (!is.element(fitstr,supported	))
    stop("Supported linear regression functions are:", supported)
  lmfit <- is.element(fitstr,"lm")

  # obtain series name
  data.name <- deparse(substitute(x))

  # if SDF smoothing is requested, force removal of DC component
  if (is.element(method,"smoothed"))
    dc <- FALSE

  # calculate single-sided SDF estimate with a normalized
  # spectral resolution of 1/N
  sdf <- sapa::SDF(x, method=sdf.method, single.sided=TRUE, npad=length(x), ...)
  attr(sdf,"series.name") <- data.name

  # associate freq.max normalized frequency
  # with corresponding index into SDF
  # frequency vector
  df <- attr(sdf, "deltaf")
  freq.indx.max <- ceiling(freq.max/df)
  freq.indx.min <- ifelse1(dc, 1, 2)

  if (is.element(method,c("standard","smoothed"))){

    # set up variables for linear regression:
    xvar <- log10(freq.indx.min:freq.indx.max)
    yvar <- log10(sdf[freq.indx.min:freq.indx.max])

    if (is.element(method,"smoothed")){

  	  # check n.block, the number of logarithmic paritions
  	  # of the SDF. if null, use the floor(log2(n.sample))
  	  # to maintain similarity to DWT decompositions
  	  if (is.null(n.block))
  	    n.block <- as.integer(floor(logb(length(x),base=2)))
  	  checkScalarType(n.block,"integer")

      # partition SDF into nonoverlapping dyadic length
      # blocks and obtain block averages
      lxvar <- length(xvar)
      delx  <- (xvar[lxvar] - xvar[1])/n.block
      xbox  <- c(xvar[1] + ((1:n.block) - 0.5) * delx)
      xbdy  <- as.numeric(c(xvar[1], (xbox + 0.5 * delx)))

      ybox <- vector(mode="double", length=n.block)
      for(k in 1:n.block){
        ybox[k] <- mean(yvar[(xvar >= xbdy[k]) & (xvar <= xbdy[k + 1])])
      }

      xvar <- xbox
      yvar <- ybox
    }
  
    w <- NULL # quell R CMD check scoping issue
  
    logfit <- ifelse1(lmfit,
      fit(y ~ 1 + x, data=data.frame(x=xvar, y=yvar, w=weight(yvar)), weights=w),
      fit(y ~ 1 + x, data=data.frame(x=xvar, y=yvar)))

    # perform trimmed least squares fit to get 1 - 2 H
    exponent <- logfit$coefficients["x"]
    H <- (1 - exponent) / 2
  }
  else if (is.element(method,"robinson")){

    n.freq <- length(sdf)

    # form cumulative sum of periodogram
    #spec.low <- sdf[freq.indx.min:freq.indx.max]
    spec.low <- sdf[1:freq.indx.max]
    csumspec <- cumsum(spec.low)

    # Set initial value of q, and a large enough delH to ensure
    # at least one iteration
    parmq <- 0.6
    delH  <- 1

    # Evaluate Robinson estimate of H
    while(abs(delH) > 0.001){

      qm    <- floor(parmq * freq.indx.max)
      NewH1 <- 1 - log(csumspec[qm]/csumspec[freq.indx.max])/(2 * log(parmq))

      if (NewH1 >= 0.75){

        # Since H >= 0.75, use suboptimal expression;
        # take avg of derivative over three nearest frequencies
        NewH1 <- 1 - (((pi * qm)/n.freq) * mean(spec.low[(
          qm - 2):qm]))/csumspec[qm]

        if (NewH1 >= 0.75){

          # terminate with current value;
          # it is the best we can do for now
          AvgH <- NewH1

          # set flag
          parmq <- 1
        }
        else{

          # do another iteration using a smaller q
          AvgH  <- NewH1 - 1
          parmq <- 0.9 * parmq
        }
      }
      else{

        # Since H < 0.75, we can iterate
        # to optimal q(H). First calculate H(Old q)
        NewH2 <- 0.75 - exp((0.01829 - parmq)/0.1292)

        # average previous H and H(Old q)
        AvgH <- mean(c(NewH1, NewH2))

        # obtain new q=q(AvgH)
        parmq <- 0.01829 - 0.1292 * log(0.75 - AvgH)
      }

      delH <- NewH1 - AvgH
      xvar <- yvar <- logfit <- NULL
    }

    H <- NewH1

    # estimate stdev of the inferred H
    #StdofH <- sqrt(Hvar(H, parmq)/freq.indx.max)

    # And drop H, optimum parameter q,
    # and the estimated std of H
    #drop(list(Hparm=NewH1, Qopt=parmq, StdH=StdofH))
  } # robinson

  # return Hurst coefficient estimate
  names(H) <- "H"

  # return the result as a fractal block exponent
  fractalBlock(
    domain        = "frequency",
    estimator     = "Hurst coefficient via regression of nonparametric SDF estimate",
    exponent      = H,
    exponent.name = "H",
    scale         = xvar,
    stat          = yvar,
    stat.name     = switch(method,standard="SDF",smoothed="mean of SDF blocks",robinson="Robinson Integration"),
    detrend       = NULL,
    overlap       = NULL,
    data.name     = data.name,
    sum.order     = 0,
    series        = asVector(x),
    logfit        = logfit,
    sdf           = sdf)
}

##
# RoverS
##

"RoverS" <- function(x, n.block.min=2, scale.ratio=2, scale.min=8)
{
  checkScalarType(n.block.min,"integer")
  checkScalarType(scale.ratio,"numeric")
  checkScalarType(scale.min,"numeric")

  # define local functions
  "Rstat" <- function(cSumSeries, timeValue, shiftValue){

    # Provide storage for the differenced cumulative sum:
    difSum <- vector(mode="double", length=shiftValue + 1)

    # Evaluate the difSum vector and the scalar difference dSumk:
    difSum <- cSumSeries[timeValue:(timeValue + shiftValue)] - cSumSeries[timeValue]
    dSumk <- cSumSeries[timeValue + shiftValue] - cSumSeries[timeValue]

    # Evaluate max and min range:
    RMax <- max(difSum - (c(timeValue:(timeValue + shiftValue)) * dSumk)/shiftValue)
    Rmin <- min(difSum - (c(timeValue:(timeValue + shiftValue)) * dSumk)/shiftValue)

    # return the R-statistic
    return(RMax - Rmin)
  }

  # convert data, if necessary, to signalSeries class:
  x <- wmtsa::create.signalSeries(x)

  # find length of series:
  n.sample <- length(x@data)

  # and find number of data blocks that can be used:
  nBlocksMax <- floor(n.sample/n.block.min)

  # find max number of times number of blocks can be increased,
  # and provide storage for result vectors:
  kmax <- floor(1 + log(n.sample/(n.block.min * scale.min))/log(scale.ratio))
  RoS <- vector(mode="double", length=kmax * nBlocksMax)
  scale <- vector(mode="double", length=kmax * nBlocksMax)

  # do cumulative sum of time series once and for all:
  cSumser <- cumsum(x@data)

  # set up initial values for while loop:
  numBlocks <- n.block.min
  k <- 1
  indx <- 0

  # and begin loop to generate R/S statistic vs k:
  while(k <= kmax) {

    nBlocks <- floor(n.sample/numBlocks)
    nUpper  <- nBlocks * numBlocks
    imax    <- floor(nBlocks - 1 - 1/numBlocks)

    for (i in (1:imax)){

      indx <- indx + 1
      t <- numBlocks * i + 1
      R <- Rstat(cSumser, t, numBlocks)
      S <- variance(x[(t + 1):(t + numBlocks)], unbiased=FALSE)
      S <- sqrt(S)
      scale[indx] <- numBlocks
      RoS[indx] <- R/S
    }

    # increase numBlocks and k; then (maybe) go around again:

    numBlocks <- numBlocks * scale.ratio
    k <- k + 1
  }

  # do linear regression over all points that have positive R/S:
  xvar <- log10(scale[RoS > 0])
  yvar <- log10(RoS[RoS > 0])

  # remove NA values
  xbad <- which(is.na(xvar))
  ybad <- which(is.na(yvar))

  bad <- union(xbad, ybad)

  if (length(bad)){
   xvar <- xvar[-bad]
   yvar <- yvar[-bad]
  }
  fitcoef <- lsfit(xvar, yvar)$coef[1:2]

  # and evaluate H:
  RoverSH <- fitcoef[2]
  RoverSH
}


###
## pgramReg
###
#
#"pgramReg" <- function(sdf, freq.indx.min=2, freq.indx.max=32, plot=FALSE)
#{
#  # Convert data, if necessary, to signalSeries class:
#  x    <- create.signalSeries(sdf)
#  spec <- x@data
#
#  # set up variables for linear regression:
#  xvar <- log10(freq.indx.min:freq.indx.max)
#  yvar <- log10(spec[freq.indx.min:freq.indx.max])
#
#  fitcoef <- lsfit(xvar, yvar)$coef[1:2]
#
#  # return Hurst coefficient estimate
#  (1 - fitcoef[2]) / 2
#}


#"pgramRegMod" <- function(sdf, freq.indx.min=2, freq.indx.max=32, Nbox=8, plot=FALSE)
#{
#  # convert data, if necessary, to signalSeries class:
#  x    <- create.signalSeries(sdf)
#  spec <- x@data
#
#  # set up variables for linear regression:
#  xvar <- log10(freq.indx.min:freq.indx.max)
#  yvar <- log10(sdf[freq.indx.min:freq.indx.max])
#  ybox <- vector(mode="double", length=Nbox)
#
#  # assign to boxes and get box averages:
#  lxvar <- length(xvar)
#  delx  <- (xvar[lxvar] - xvar[1])/Nbox
#  xbox  <- c(xvar[1] + ((1:Nbox) - 0.5) * delx)
#  xbdy  <- c(xvar[1], (xbox + 0.5 * delx))
#  for(k in 1:Nbox){
#    ybox[k] <- mean(yvar[(xvar >= as.numeric(xbdy[k])) & (xvar <=
#      as.numeric(xbdy[k + 1]))])
#  }
#
#  # perform trimmed least squares fit to get 1 - 2 H; solve for H
#  fitcoef <- ltsreg(xbox, ybox)$coef[1:2]
#
#  if(plot == T) {
#    plot(xbox, ybox, xlab="log10(Frequency Index)", ylab =
#      "log10(SDF)")
#    panel.abline(fitcoef)
#  }
#
#  # return Hurst coefficient estimate
#  (1 - fitcoef[2]) / 2
#}

##
# RobInt
##

#"RobInt" <- function(sdf, freq.max=0.2)
#{
#  #TODO: what is Hvar(q=1) supposed to be?
#  # define local functions
#  "Hvar" <- function(H, q)
#  {
#    # Evaluates variance estimate of distribution of parameter H
#    # Ref: I. Lobato and P. M. Robinson,
#    temp1 <- (1 + 1/q - 2 * (q^(1 - 2 * H)))/((log(q))^2)
#    (temp1 * ((1 - H)^2))/(3 - 4 * H)
#  }
#
#  if (!is(sdf, "SDF"))
#    stop("sdf must be an object of class \"SDF\"")
#  if (!attr(sdf,"single.sided"))
#    stop("sdf must be single-sided")
#
#  checkScalarType(freq.max,"numeric")
#  if (freq.max <= 0.0 || freq.max >= 0.5)
#    stop("freq.max must be on the normalized frequency range (0,0.5)")
#
#  # associate freq.max normalized frequency
#  # with corresponding index into SDF
#  # frequency vector
#  df <- attr(sdf, "deltaf")
#  freq.indx.max <- ceiling(freq.max/df)
#  spec  <- as.vector(sdf)
#  n.freq <- length(spec)
#
#  #TODO: does RobInt need to be a periodogram?
#  # form cumulative sum of periodogram
#  spec.low <- spec[seq(freq.indx.max)]
#  csumspec <- cumsum(spec.low)
#
#  # Set initial value of q, and a large enough delH to ensure
#  # at least one iteration
#  parmq <- 0.6
#  delH  <- 1
#
#  # Evaluate Robinson estimate of H
#  while(abs(delH) > 0.001){
#
#    qm    <- floor(parmq * freq.indx.max)
#    NewH1 <- 1 - log(csumspec[qm]/csumspec[freq.indx.max])/(2 * log(parmq))
#
#    if (NewH1 >= 0.75){
#
#      # Since H >= 0.75, use suboptimal expression;
#      # take avg of derivative over three nearest frequencies
#      NewH1 <- 1 - (((pi * qm)/n.freq) * mean(spec.low[(
#        qm - 2):qm]))/csumspec[qm]
#
#      if (NewH1 >= 0.75){
#        # terminate with current value;
#        # it is the best we can do for now
#        AvgH <- NewH1
#
#        # set flag
#        parmq <- 1
#      }
#      else{
#
#        # do another iteration using a smaller q
#        AvgH  <- NewH1 - 1
#        parmq <- 0.9 * parmq
#      }
#    }
#    else{
#
#        # Since H < 0.75, we can iterate
#        # to optimal q(H). First calculate H(Old q)
#        NewH2 <- 0.75 - exp((0.01829 - parmq)/0.1292)
#
#        # average previous H and H(Old q)
#        AvgH <- mean(c(NewH1, NewH2))
#
#        # obtain new q=q(AvgH)
#        parmq <- 0.01829 - 0.1292 * log(0.75 - AvgH)
#    }
#
#    delH <- NewH1 - AvgH
#  }
#
#
#  # estimate stdev of the inferred H
#  StdofH <- sqrt(Hvar(NewH1, parmq)/freq.indx.max)
#
#  # And drop H, optimum parameter q,
#  # and the estimated std of H
#  drop(list(Hparm=NewH1, Qopt=parmq, StdH=StdofH))
#}


###
## disWhittle
###
#
#"disWhittle" <- function(sdf, NFT, freq.indx.min=2, freq.indx.max=floor((NFT + 1)/2))
#{
#  # convert datato signalSeries class
#  x <- create.signalSeries(sdf)
#  Sx <- x@data
#
#  # define function to evaluate the weighted sum of the spectrum
#  Isumd <- function(delta, Sx, n.freq, Mindx, Maxdx)
#  {
#  	xatt <- attributes(Sx)
#  	f    <- xatt$frequency
#  	N    <- xatt$n.sample
#  	Np   <- (N - 1) %/% 2
#
#  #  f     <- c(1:(Maxdx - Mindx + 1)) / n.freq
#  #  Nstar <- length(f)
#
#    Isum1 <- mean(as.numeric(Sx) * (4 * (sin(pi * f))^2)^delta)
#    log(Isum1) + 2 * delta * (((N/(2 * Np)) - 1) *
#      log(2) - log(2 * N)/(2 * Np))
#  }
#
#  # call optimize() function to find the best value of delta
#  dfinal <- optimize(Isumd, lower=0, upper=0.5, tol=0.0001,
#    maximum=F, Sx=Sx[freq.indx.min:freq.indx.max], n.freq=NFT,
#    Mindx=freq.indx.min, Maxdx=freq.indx.max)
#
#  # return Hurst coefficient estimate
#  dfinal$minimum + 0.5
#}
#
#######################################################################################
#
#"Whittle" <- function(sdf, NFT, freq.indx.min=2, freq.indx.max=floor((NFT + 1)/2))
#{
#  # Convert data, if necessary, to signalSeries class:
#  x <- create.signalSeries(sdf)
#  spec <- x@data
#
#  #  Here is the function to evaluate the weighted sum of the spectrum:
#  Isumd <- function(dparm, spec, n.freq, Mindx, Maxdx)
#  {
#    f <- c(1:(Maxdx - Mindx + 1))/n.freq
#    Isumd <- sum(spec * (4 * (sin(pi * f))^2)^dparm)/n.freq
#    Isumd
#  }
#
#  # Just call the S-Plus "optimize" function to find the best
#  #  value of delta:
#  dfinal <- optimize(Isumd, lower=0, upper=0.5, tol=0.001,
#    maximum=F, spec=spec[freq.indx.min:freq.indx.max], n.freq=NFT,
#    Mindx=freq.indx.min, Maxdx=freq.indx.max)
#
#  # return Hurst coefficient estimate
#  dfinal$minimum + 0.5
#}
