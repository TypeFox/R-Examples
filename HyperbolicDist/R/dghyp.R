### Function to calculate the density of the
### generalized hyperbolic distribution
dghyp <- function(x, Theta)
{
  if(length(Theta) == 4) Theta <- c(1, Theta)
  if(length (Theta) != 5){
    stop("parameter vector must contain 5 values") 
  }
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]
  
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (abs(beta) >= alpha) {
    stop("absolute value of beta must be less than alpha")
  }
  gamma <- sqrt(alpha^2 - beta^2)

  ## Argument of Bessel K function in numerator
  y <- alpha*sqrt(delta^2 + (x - mu)^2)
  bx <- beta*(x - mu)
  ## Deal with underflow in ratio of Bessel K functions
  ## besselK underflows for x > 740
  ## Use exponentially scaled besselK
  if(delta*gamma > 700){# underflow in constant part
    expTerm <- exp(delta*gamma - y + bx)
    besselRatio <- besselK(x = y, nu = lambda -1/2, expon.scaled = TRUE)/
      besselK(x = delta*gamma, nu = lambda, expon.scaled = TRUE)  
    expAndBessel <- expTerm*besselRatio
  }else{
    expAndBessel <- ifelse(y > 700|bx > 700, # underflow in variable part
                           exp(delta*gamma - y + bx)*
                           besselK(x = y, nu = lambda -1/2,
                                   expon.scaled = TRUE)/
                           besselK(x = delta*gamma, nu = lambda,
                                   expon.scaled = TRUE),
                          exp(bx)*besselK(x = y, nu = lambda -1/2)/
                          besselK(x = delta*gamma, nu = lambda))
  }
  dens <- (y/alpha)^(lambda - 1/2)*((gamma/delta)^lambda)*
    alpha^(1/2-lambda)*expAndBessel/sqrt(2*pi)
  dens
} ## End of dghyp()


### Cumulative distribution function of the generalized hyperbolic
### New version intended to give guaranteed accuracy
### Uses exponential approximation in tails
### SafeIntegrate() is used over six parts in the middle of the distribution
### Calls ghypBreaks to determine the breaks
###
### DJS 07/01/07
pghyp <- function(q, Theta, small = 10^(-6), tiny = 10^(-10),
                    deriv = 0.3, subdivisions = 100,
                    accuracy = FALSE, ...){

  if(length(Theta) == 4) Theta <- c(1, Theta)
  if(length (Theta) != 5){
    stop("parameter vector must contain 5 values") 
  }
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]
  
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (abs(beta) >= alpha) {
    stop("absolute value of beta must be less than alpha")
  }

  bks <- ghypBreaks(Theta, small, tiny, deriv, ...)
  xTiny <- bks$xTiny
  xSmall <- bks$xSmall
  lowBreak <- bks$lowBreak
  highBreak <- bks$highBreak
  xLarge <- bks$xLarge
  xHuge <- bks$xHuge
  modeDist <- bks$modeDist

  qSort <- sort(q)
  qTiny <- which(qSort < xTiny)
  qSmall <- which(qSort < xSmall)
  qLow <- which(qSort < lowBreak)
  qLessEqMode <- which(qSort <= modeDist)
  qGreatMode <- which(qSort > modeDist)
  qHigh <- which(qSort > highBreak)
  qLarge <- which(qSort > xLarge)
  qHuge <- which(qSort > xHuge)

  ## Break indices into 8 groups: beware of empty groups
  if (length(qLow) > 0) qLessEqMode <- qLessEqMode[qLessEqMode > max(qLow)]
  if (length(qHigh) > 0) qGreatMode <- qGreatMode[qGreatMode < min(qHigh)]
  if (length(qSmall) > 0) qLow <- qLow[qLow > max(qSmall)]
  if (length(qLarge) > 0) qHigh <- qHigh[qHigh < min(qLarge)]
  if (length(qTiny) > 0) qSmall <- qSmall[qSmall > max(qTiny)]
  if (length(qHuge) > 0) qLarge <- qLarge[qLarge < min(qHuge)]
  intFun <- rep(NA, length(q))
  if (length(qTiny) > 0) intFun[qTiny] <- 0
  if (length(qHuge) > 0) intFun[qHuge] <- 1
  intErr <- rep(NA, length(q))
  if (length(qTiny) > 0) intErr[qTiny] <- tiny
  if (length(qHuge) > 0) intErr[qHuge] <- tiny


  ## Use safeIntegrate function between xTiny and xHuge in 6 sections
  dghypInt <- function(q){ 
    dghyp(q, Theta)
  }

  ## Calculate integrals and errors to cut points
  resSmall <- safeIntegrate(dghypInt, xTiny, xSmall, subdivisions, ...)
  resLarge <- safeIntegrate(dghypInt, xLarge, xHuge, subdivisions, ...)
  intSmall <- resSmall$value
  intLarge <- resLarge$value
  errSmall <- tiny + resSmall$abs.error
  errLarge <- tiny + resLarge$abs.error
  resLow <- safeIntegrate(dghypInt, xSmall, lowBreak, subdivisions, ...)
  resHigh <- safeIntegrate(dghypInt, highBreak, xLarge, subdivisions, ...)
  intLow <- intSmall + resLow$value
  intHigh <- intLarge + resHigh$value
  errLow <- errSmall + resLow$abs.error
  errHigh <- errLarge + resHigh$abs.error

  for (i in qSmall){
    intRes <- safeIntegrate(dghypInt, xTiny, qSort[i], subdivisions, ...)
    intFun[i] <- intRes$value
    intErr[i] <- intRes$abs.error + tiny
  }
  for (i in qLarge){
    intRes <- safeIntegrate(dghypInt, qSort[i], xHuge, subdivisions, ...)
    intFun[i] <-  1- intRes$value 
    intErr[i] <- intRes$abs.error + tiny
  }
  for (i in qLow){
    intRes <- safeIntegrate(dghypInt, xSmall, qSort[i], subdivisions, ...)
    intFun[i] <- intRes$value + intSmall
    intErr[i] <- intRes$abs.error + errSmall
  }
  for (i in qHigh){
    intRes <- safeIntegrate(dghypInt, qSort[i], xLarge, subdivisions, ...)
    intFun[i] <-  1- intRes$value - intLarge
    intErr[i] <- intRes$abs.error + errLarge
  }
  for (i in qLessEqMode){
    intRes <- safeIntegrate(dghypInt, lowBreak, qSort[i], subdivisions, ...)
    intFun[i] <- intRes$value + intLow
    intErr[i] <- intRes$abs.error + errLow
  }
  for (i in qGreatMode){
    intRes <- safeIntegrate(dghypInt, qSort[i], highBreak, subdivisions, ...)
    intFun[i] <-  1- intRes$value - intHigh
    intErr[i] <- intRes$abs.error + errLarge
  }

  if (!accuracy){
    return(intFun[rank(q)])
  }else{
    return(list(value=intFun[rank(q)], error=intErr[rank(q)]))
  }
} ## End of pghyp()

### Cumulative distribution function of the generalized hyperbolic
### New version intended to give guaranteed accuracy
### qghyp using breaks as for pghyp and splines as in original qghyp
###
### DJS 06/09/06
qghyp <- function(p, Theta, small = 10^(-6), tiny = 10^(-10),
                    deriv = 0.3, nInterpol = 100, subdivisions = 100, ...){ 
  if(length(Theta) == 4) Theta <- c(1, Theta)
  if(length (Theta) != 5){
    stop("parameter vector must contain 5 values") 
  }
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]
  
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (abs(beta) >= alpha) {
    stop("absolute value of beta must be less than alpha")
  }

  bks <- ghypBreaks(Theta, small, tiny, deriv, ...)
  xTiny <- bks$xTiny
  xSmall <- bks$xSmall
  lowBreak <- bks$lowBreak
  highBreak <- bks$highBreak
  xLarge <- bks$xLarge
  xHuge <- bks$xHuge
  modeDist <- bks$modeDist

  yTiny <- pghyp(xTiny, Theta)
  ySmall <- pghyp(xSmall, Theta)
  yLowBreak <- pghyp(lowBreak, Theta)
  yHighBreak <- pghyp(highBreak, Theta)
  yLarge <- pghyp(xLarge, Theta)
  yHuge <- pghyp(xHuge, Theta)
  yModeDist <- pghyp(modeDist, Theta)

  pSort <- sort(p)
  pSmall <- which(pSort < pghyp(xSmall, Theta))
  pTiny <- which(pSort < pghyp(xTiny, Theta))
  pLarge <- which(pSort > pghyp(xLarge, Theta))
  pHuge <- which(pSort > pghyp(xHuge, Theta))
  pLow <- which(pSort < pghyp(lowBreak, Theta))
  pHigh <- which(pSort > pghyp(highBreak, Theta))
  pLessEqMode <- which(pSort <= pghyp(modeDist, Theta))
  pGreatMode <- which(pSort > pghyp(modeDist, Theta))

  ## Break indices into 8 groups: beware of empty groups
  if (length(pLow) > 0) pLessEqMode <- pLessEqMode[pLessEqMode > max(pLow)]
  if (length(pHigh) > 0) pGreatMode <- pGreatMode[pGreatMode < min(pHigh)]
  if (length(pSmall) > 0) pLow <- pLow[pLow > max(pSmall)]
  if (length(pLarge) > 0) pHigh <- pHigh[pHigh < min(pLarge)]
  if (length(pTiny) > 0) pSmall <- pSmall[pSmall > max(pTiny)]
  if (length(pHuge) > 0) pLarge <- pLarge[pLarge < min(pHuge)]
  qSort <- rep(NA, length(pSort))
  if (length(pTiny) > 0) qSort[pTiny] <- -Inf
  if (length(pHuge) > 0) qSort[pHuge] <- Inf


  if (length(pTiny) > 0){
    for (i in pTiny){
      zeroFun<-function(x){ 
        pghyp(x,Theta) - pSort[i] 
      }
      interval <- c(xTiny - (xSmall - xTiny),xTiny)
      while(zeroFun(interval[1])*zeroFun(interval[2])>0) {
        interval[1] <- interval[1] - (xSmall - xTiny)
      }
      qSort[i] <- uniroot(zeroFun,interval)$root
    }
  }
  if (length(pSmall) > 0){
    xValues <- seq(xTiny, xSmall, length = nInterpol)
    pghypValues <- pghyp(xValues, Theta, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    pghypSpline <- splinefun(xValues, pghypValues)
    for(i in pSmall){ 
      zeroFun<-function(x){ 
        pghypSpline(x) - pSort[i] 
      }
      if (zeroFun(xTiny) >= 0){
        qSort[i] <- xTiny
      }else{
        if (zeroFun(xSmall) <= 0){
          qSort[i] <- xSmall
        }else{
          qSort[i] <- uniroot(zeroFun, interval = c(xTiny,xSmall), ...)$root
        }
      }
    }
  }
  if (length(pLow) > 0){
    xValues <- seq(xSmall, lowBreak, length = nInterpol)
    pghypValues <- pghyp(xValues, Theta, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    pghypSpline <- splinefun(xValues, pghypValues)
    for(i in pLow){ 
      zeroFun<-function(x){ 
        pghypSpline(x) - pSort[i] 
      }
      if (zeroFun(xSmall) >= 0){
        qSort[i] <- xSmall
      }else{
        if (zeroFun(lowBreak) <= 0){
          qSort[i] <- lowBreak
        }else{
          qSort[i] <- uniroot(zeroFun, interval = c(xSmall,lowBreak), ...)$root
        }
      }
    }
  }
  if (length(pLessEqMode) > 0){
    xValues <- seq(lowBreak, modeDist, length = nInterpol)
    pghypValues <- pghyp(xValues, Theta, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    pghypSpline <- splinefun(xValues, pghypValues)
    for(i in pLessEqMode){ 
      zeroFun<-function(x){ 
        pghypSpline(x) - pSort[i] 
      }
      if (zeroFun(lowBreak) >= 0){
        qSort[i] <- lowBreak
      }else{
        if (zeroFun(modeDist) <= 0){
          qSort[i] <- modeDist
        }else{
          qSort[i] <-
            uniroot(zeroFun, interval = c(lowBreak,modeDist), ...)$root
        }
      }
    }
  }
  if (length(pGreatMode) > 0){
    xValues <- seq(modeDist, highBreak, length = nInterpol)
    pghypValues <- pghyp(xValues, Theta, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    pghypSpline <- splinefun(xValues, pghypValues)
    for(i in pGreatMode){ 
      zeroFun<-function(x){ 
        pghypSpline(x) - pSort[i] 
      }
      if (zeroFun(modeDist) >= 0){
        qSort[i] <- modeDist
      }else{
        if (zeroFun(highBreak) <= 0){
          qSort[i] <- highBreak
        }else{
          qSort[i] <-
            uniroot(zeroFun, interval = c(modeDist,highBreak), ...)$root
        }
      }
    }
  }
  if (length(pHigh) > 0){
    xValues <- seq(highBreak, xLarge, length = nInterpol)
    pghypValues <- pghyp(xValues, Theta, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    pghypSpline <- splinefun(xValues, pghypValues)
    for(i in pHigh){ 
      zeroFun<-function(x){ 
        pghypSpline(x) - pSort[i] 
      }
      if (zeroFun(highBreak) >= 0){
        qSort[i] <- highBreak
      }else{
        if (zeroFun(xLarge) <= 0){
          qSort[i] <- xLarge
        }else{
          qSort[i] <-
            uniroot(zeroFun, interval = c(highBreak,xLarge), ...)$root
        }
      }
    }
  }
  if (length(pLarge) > 0){
    xValues <- seq(xLarge, xHuge, length = nInterpol)
    pghypValues <- pghyp(xValues, Theta, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    pghypSpline <- splinefun(xValues, pghypValues)
    for(i in pLarge){ 
      zeroFun<-function(x){ 
        pghypSpline(x) - pSort[i] 
      }
      if (zeroFun(xLarge) >= 0){
        qSort[i] <- xLarge
      }else{
        if (zeroFun(xHuge) <= 0){
          qSort[i] <- xHuge
        }else{
          qSort[i] <-
            uniroot(zeroFun, interval = c(xLarge,xHuge), ...)$root
        }
      }
    }
  }        
  if (length(pHuge) > 0){
    for (i in pHuge){
      zeroFun<-function(x){ 
        pghyp(x,Theta) - pSort[i] 
      }
      interval <- c(xHuge,xHuge + (xHuge - xLarge))
      while(zeroFun(interval[1])*zeroFun(interval[2])>0) {
        interval[1] <- interval[1] + (xHuge - xLarge)
      }
      qSort[i] <- uniroot(zeroFun,interval)$root
    }
  }

  return(qSort[rank(p)]) 
} # End of qghyp()

### Derivative of the density
ddghyp <- function(x, Theta){
  if(length(Theta) == 4) Theta <- c(1, Theta)
  if(length (Theta) != 5){
    stop("parameter vector must contain 5 values") 
  }
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]

  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (abs(beta) >= alpha) {
    stop("absolute value of beta must be less than alpha")
  }
  ## Terms for simplification of programming
  t1 <- sqrt(delta^2 + (x - mu)^2)
  t2 <- sqrt(alpha^2 - beta^2)
  t3 <- besselK(x = alpha*t1, nu = lambda - 0.5)
  t4 <- besselK(x = alpha*t1, nu = lambda + 0.5)
  
  ddghyp <- (t3*(beta*delta^2 + (2*lambda - 1)*(x - mu) + beta*(x - mu)^2) -
             t4*alpha*t1*(x - mu))*
               exp(beta*(x - mu))*t1^(lambda - (5/2))*t2^lambda/
                          (sqrt(2*pi)*alpha^(lambda -  1/2)*delta^lambda*
                           besselK(x = delta*t2, nu = lambda))
                      
  ddghyp
} ## End of ddghyp()

### Function to set up breaks for pghyp and qghyp
ghypBreaks <- function(Theta, small = 10^(-6), tiny = 10^(-10),
                         deriv = 0.3, ...){
  if(length(Theta) == 4) Theta <- c(1, Theta)
  if(length (Theta) != 5){
    stop("parameter vector must contain 5 values") 
  }
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]

  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (abs(beta) >= alpha) {
    stop("absolute value of beta must be less than alpha")
  }
  
  xTiny <- ghypCalcRange(Theta, tiny, density = TRUE)[1] 
  xSmall <- ghypCalcRange(Theta, small, density = TRUE)[1]
  xLarge <- ghypCalcRange(Theta, small, density = TRUE)[2]
  xHuge <- ghypCalcRange(Theta, tiny, density = TRUE)[2]

  modeDist <- ghypMode(Theta)
  ## Determine break points, based on size of derivative
  xDeriv <- seq(xSmall, modeDist, length.out = 101)
  derivVals <- ddghyp(xDeriv,Theta)
  maxDeriv <- max(derivVals)
  minDeriv <- min(derivVals)
  breakSize <- deriv*maxDeriv
  breakFun <- function(x){ 
    ddghyp(x, Theta) - breakSize
  }
  if ((maxDeriv < breakSize)||(derivVals[1] > breakSize)){
    lowBreak <- xSmall
  }else{
    whichMaxDeriv <- which.max(derivVals)
    lowBreak <- uniroot(breakFun,c(xSmall,xDeriv[whichMaxDeriv]))$root
  }

  xDeriv <- seq(modeDist, xLarge, length.out = 101)
  derivVals <- -ddghyp(xDeriv,Theta)
  maxDeriv <- max(derivVals)
  minDeriv <- min(derivVals)
  breakSize <- deriv*maxDeriv
  breakFun <- function(x){ 
    - ddghyp(x, Theta) - breakSize
  }
  if ((maxDeriv < breakSize)||(derivVals[101] > breakSize)){
    highBreak <- xLarge
  }else{
    whichMaxDeriv <- which.max(derivVals)
    highBreak <- uniroot(breakFun,c(xDeriv[whichMaxDeriv],xLarge))$root
  }
  breaks <- c(xTiny,xSmall,lowBreak,highBreak,xLarge,xHuge,modeDist)
  breaks <- list(xTiny = breaks[1], xSmall =  breaks[2],
                 lowBreak =  breaks[3], highBreak =  breaks[4],
                 xLarge =  breaks[5], xHuge =  breaks[6],
                 modeDist =  breaks[7])
  return(breaks)
} ## End of ghypBreaks()


### Function to generate random observations from a
### (generalized) hyperbolic distribution using the
### mixing property of the generalized inverse
### Gaussian distribution.
rghyp <- function(n, Theta){
  if(length(Theta) == 4) Theta <- c(1,Theta)
  if(length (Theta) != 5){
    stop("parameter vector must contain 5 values") 
  }
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]
   
  if (alpha <= 0) {
    stop("alpha must be positive")
  }
  if (delta <= 0) {
    stop("delta must be positive")
  }
  if (abs(beta) >= alpha) {
    stop("absolute value of beta must be less than alpha")
  }

  chi <- delta^2
  psi <- alpha^2 - beta^2

  if(lambda == 1){
    X <- rgig1(n, c(lambda,chi,psi))
  } else{
    X <- rgig(n, c(lambda,chi,psi))
  }
  
  sigma <- sqrt(X)
  Z <- rnorm(n)
  Y <- mu + beta*sigma^2 + sigma*Z

  Y
} ## End of rghyp()
