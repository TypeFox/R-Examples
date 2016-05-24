### Functions for the hyperbolic distribution
### Density of the hyperbolic distribution
dhyperb <- function(x, Theta, KNu = NULL, logPars = FALSE){ 
  if (length(Theta) != 4){
    stop("parameter vector must contain 4 values")
  }
  Theta <- as.numeric(Theta)
  if(logPars == TRUE){ 
    hyperbPi <- Theta[1] 
    lZeta <- Theta[2] 
    lDelta <- Theta[3] 
    mu <- Theta[4] 

    if(is.null(KNu)){ 
      KNu <- besselK(exp(lZeta), nu = 1)     
    } 
     
    hyperbDens <- (2*exp(lDelta)* 
                    sqrt(1 + hyperbPi^2)*KNu)^(-1)* 
                    exp(-exp(lZeta)* 
            (sqrt(1 + hyperbPi^2)* 
                    sqrt(1 + ((x - mu)/
            exp(lDelta))^2) - 
                    hyperbPi*(x - mu)/
            exp(lDelta))) 
  } 
  else{ 
    hyperbPi <- Theta[1] 
    zeta <- Theta[2] 
    delta <- Theta[3] 
    mu <- Theta[4] 
     
    if(is.null(KNu)){ 
      KNu <- besselK(zeta, nu=1)     
    }   

    hyperbDens <- (2*delta* 
            sqrt(1 + hyperbPi^2)*KNu)^(-1)* 
                    exp(-zeta*(sqrt(1 + hyperbPi^2)* 
                    sqrt(1 + ((x - mu)/delta)^2) - 
                    hyperbPi*(x - mu)/delta)) 
  } 
  as.numeric(hyperbDens) 
} ## End of dhyperb()


### Cumulative distribution function of the hyperbolic distribution
### New version intended to give guaranteed accuracy
### Uses exponential approximation in tails
### Integrate() is used over four parts in the middle of the distribution
### This version calls hyperbBreaks to determine the breaks
###
### DJS 05/09/06
phyperb <- function(q, Theta, small = 10^(-6), tiny = 10^(-10),
                    deriv = 0.3, subdivisions = 100,
                    accuracy = FALSE, ...){

  if (length(Theta) != 4){
   stop("parameter vector must contain 4 values")
  }
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3]
  mu <- Theta[4]
  if(zeta <= 0) stop("zeta must be positive")
  if(delta <= 0) stop("delta must be positive")
  ## standardise distribution to delta=1, mu =0
  q <- (q - mu)/delta
  delta <- 1
  mu <- 0
  Theta[3] <- delta
  Theta[4] <- mu
  KNu <- besselK(zeta, nu = 1)
  phi <- as.numeric(hyperbChangePars(1, 3, Theta)[1])
  gamma <- as.numeric(hyperbChangePars(1, 3, Theta)[2])
  const <- 1/(2*(1 + hyperbPi^2)^(1/2)*KNu)

  bks <- hyperbBreaks(Theta, small, tiny, deriv, ...)
  xTiny <- bks$xTiny
  xSmall <- bks$xSmall
  lowBreak <- bks$lowBreak
  highBreak <- bks$highBreak
  xLarge <- bks$xLarge
  xHuge <- bks$xHuge
  modeDist <- bks$modeDist

  qSort <- sort(q)
  qSmall <- which(qSort < xSmall)
  qTiny <- which(qSort < xTiny)
  qLarge <- which(qSort > xLarge)
  qHuge <- which(qSort > xHuge)
  qLow <- which(qSort < lowBreak)
  qHigh <- which(qSort > highBreak)
  qLessEqMode <- which(qSort <= modeDist)
  qGreatMode <- which(qSort > modeDist)

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

  ## Approximate pdf by exponential between xTiny and xSmall
  if (length(qSmall) > 0){
    intFun[qSmall] <- (1/phi)*const*exp(phi*qSort[qSmall])
  }
  ## Approximate pdf by exponential between xLarge and xHuge
  if (length(qLarge) > 0){
    intFun[qLarge] <- 1 - (1/gamma)*const*exp(-gamma*qSort[qLarge])
  }
  ## Use integrate function between xSmall and xLarge
  dhypInt <- function(q){ 
    dhyperb(q, Theta, KNu)
  }

  ## Calculate integrals and errors to cut points
  intSmall <- (1/phi)*const*exp(phi*xSmall)
  intLarge <- (1/gamma)*const*exp(-gamma*xLarge)
  resLow <- safeIntegrate(dhypInt, xSmall, lowBreak, subdivisions, ...)
  resHigh <- safeIntegrate(dhypInt, highBreak, xLarge, subdivisions, ...)
  intLow <- intSmall + resLow$value
  intHigh <- intLarge + resHigh$value
  errSmall <- tiny + (1/phi)*const*exp(phi*xSmall)*
                      (1 - (phi + gamma)/(4*abs(xSmall)))
  errLarge <- tiny + (1/gamma)*const*exp(-gamma*xLarge)*
                         (1 - (phi + gamma)/(4*xLarge))
  errLow <- errSmall + resLow$abs.error
  errHigh <- errLarge + resHigh$abs.error

  for (i in qLow){
    intRes <- safeIntegrate(dhypInt, xSmall, qSort[i], subdivisions, ...)
    intFun[i] <- intRes$value + intSmall
    intErr[i] <- intRes$abs.error + errSmall
  }
  for (i in qHigh){
    intRes <- safeIntegrate(dhypInt, qSort[i], xLarge, subdivisions, ...)
    intFun[i] <-  1- intRes$value - intLarge
    intErr[i] <- intRes$abs.error + errLarge
  }
  for (i in qLessEqMode){
    intRes <- safeIntegrate(dhypInt, lowBreak, qSort[i], subdivisions, ...)
    intFun[i] <- intRes$value + intLow
    intErr[i] <- intRes$abs.error + errLow
  }
  for (i in qGreatMode){
    intRes <- safeIntegrate(dhypInt, qSort[i], highBreak, subdivisions, ...)
    intFun[i] <-  1- intRes$value - intHigh
    intErr[i] <- intRes$abs.error + errLarge
  }

  if (!accuracy){
    return(intFun[rank(q)])
  }else{
    return(list(value=intFun[rank(q)], error=intErr[rank(q)]))
  }
} ## End of phyperb()

### qhyperb using breaks as for phyperb and splines as in original qhyperb
###
### DJS 06/09/06
qhyperb <- function(p, Theta, small = 10^(-6), tiny = 10^(-10),
                    deriv = 0.3, nInterpol = 100, subdivisions = 100, ...){ 
  if (length(Theta) != 4){
   stop("parameter vector must contain 4 values")
  }
  Theta <- as.numeric(Theta)
  if(Theta[2] <= 0) stop("zeta must be positive")
  if(Theta[3] <= 0) stop("delta must be positive")
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  ## Find quantiles of standardised hyperbolic: adjust later
  delta <- 1
  mu <- 0
  ThetaStand <- c(hyperbPi,zeta,delta,mu)
  KNu <- besselK(zeta, nu = 1)
  phi <- as.numeric(hyperbChangePars(1, 3, ThetaStand)[1])
  gamma <- as.numeric(hyperbChangePars(1, 3, ThetaStand)[2])
  const <- 1/(2*delta*(1 + hyperbPi^2)^(1/2)*KNu)

  bks <- hyperbBreaks(ThetaStand, small, tiny, deriv, ...)
  xTiny <- bks$xTiny
  xSmall <- bks$xSmall
  lowBreak <- bks$lowBreak
  highBreak <- bks$highBreak
  xLarge <- bks$xLarge
  xHuge <- bks$xHuge
  modeDist <- bks$modeDist

  yTiny <- phyperb(xTiny, ThetaStand)
  ySmall <- phyperb(xSmall, ThetaStand)
  yLowBreak <- phyperb(lowBreak, ThetaStand)
  yHighBreak <- phyperb(highBreak, ThetaStand)
  yLarge <- phyperb(xLarge, ThetaStand)
  yHuge <- phyperb(xHuge, ThetaStand)
  yModeDist <- phyperb(modeDist, ThetaStand)

  pSort <- sort(p)
  pSmall <- which(pSort < phyperb(xSmall, ThetaStand))
  pTiny <- which(pSort < phyperb(xTiny, ThetaStand))
  pLarge <- which(pSort > phyperb(xLarge, ThetaStand))
  pHuge <- which(pSort > phyperb(xHuge, ThetaStand))
  pLow <- which(pSort < phyperb(lowBreak, ThetaStand))
  pHigh <- which(pSort > phyperb(highBreak, ThetaStand))
  pLessEqMode <- which(pSort <= phyperb(modeDist, ThetaStand))
  pGreatMode <- which(pSort > phyperb(modeDist, ThetaStand))

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

  ## Approximate pdf by exponential between xTiny and xSmall
  if (length(pSmall) > 0){
    qSort[pSmall] <- (1/phi)*log(phi*pSort[pSmall]/const)
  }
  ## Approximate pdf by exponential between xLarge and xHuge
  if (length(pLarge) > 0){
    qSort[pLarge] <- -(1/gamma)*log(gamma*(1-pSort[pLarge])/const)
  }

  if (length(pLow) > 0){
    xValues <- seq(xSmall, lowBreak, length = nInterpol)
    phyperbValues <- phyperb(xValues, ThetaStand, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    phyperbSpline <- splinefun(xValues, phyperbValues)
    for(i in pLow){ 
      zeroFun<-function(x){ 
        phyperbSpline(x) - pSort[i] 
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
    phyperbValues <- phyperb(xValues, ThetaStand, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    phyperbSpline <- splinefun(xValues, phyperbValues)
    for(i in pLessEqMode){ 
      zeroFun<-function(x){ 
        phyperbSpline(x) - pSort[i] 
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
    phyperbValues <- phyperb(xValues, ThetaStand, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    phyperbSpline <- splinefun(xValues, phyperbValues)
    for(i in pGreatMode){ 
      zeroFun<-function(x){ 
        phyperbSpline(x) - pSort[i] 
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
    phyperbValues <- phyperb(xValues, ThetaStand, small, tiny, deriv,
                            subdivisions = subdivisions, accuracy = FALSE)
    phyperbSpline <- splinefun(xValues, phyperbValues)
    for(i in pHigh){ 
      zeroFun<-function(x){ 
        phyperbSpline(x) - pSort[i] 
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

   return(qSort[rank(p)]*Theta[3] + Theta[4]) 
} # End of qhyperb()

### Function to generate random observations from a
### hyperbolic distribution using the
### mixing property of the generalized inverse
### Gaussian distribution and Dagpunar's algorithm
### for the generalized inverse Gaussian
rhyperb <- function(n, Theta){
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3]
  mu <- Theta[4]

  if(zeta <= 0) stop("zeta must be positive")
  if(delta <= 0) stop("delta must be positive")

  alpha <- zeta*sqrt(1 + hyperbPi^2)/delta
  beta <- zeta*hyperbPi/delta
  chi <- delta^2
  psi <- alpha^2 - beta^2

  X <- rgig1(n, c(1,chi,psi))

  sigma <- sqrt(X)
  Z <- rnorm(n)
  Y <- mu + beta*sigma^2 + sigma*Z

  Y
} ## End of rhyperb()

### Derivative of the density
ddhyperb <- function(x, Theta, KNu = NULL, ...){
  if (length(Theta) != 4){
   stop("parameter vector must contain 4 values")
  }
  if(Theta[2] <= 0) stop("zeta must be positive")
  if(Theta[3] <= 0) stop("delta must be positive")
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3] 
  mu <- Theta[4]
  if(is.null(KNu)){ 
    KNu <- besselK(zeta, nu=1)     
  } 
  ddhyp <- dhyperb(x, Theta,  KNu)*(zeta*hyperbPi/delta -
             zeta*sqrt(1 + hyperbPi^2)*(x - mu)/
             (sqrt(1 + ((x - mu)/delta)^2)*delta^2))
  ddhyp
} ## End of ddhyperb()

### Function to set up breaks for phyperb and qhyperb
hyperbBreaks <- function(Theta, small = 10^(-6), tiny = 10^(-10),
                         deriv = 0.3, ...){
  if (length(Theta) != 4){
   stop("parameter vector must contain 4 values")
  }
  if(Theta[2] <= 0) stop("zeta must be positive")
  if(Theta[3] <= 0) stop("delta must be positive")
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  ## Find quantiles of standardised hyperbolic: adjust later
  delta <- 1
  mu <- 0
  ThetaStand <- c(hyperbPi,zeta,delta,mu)

  KNu <- besselK(zeta, nu = 1)
  phi <- as.numeric(hyperbChangePars(1, 3, ThetaStand)[1])
  gamma <- as.numeric(hyperbChangePars(1, 3, ThetaStand)[2])

  const <- 1/(2*(1 + hyperbPi^2)^(1/2)*KNu)

  xSmall <- 1/phi*log(small*phi/const) 
  xTiny <- 1/phi*log(tiny*phi/const) 
  xLarge <- -1/gamma*log(small*gamma/const)
  xHuge <- -1/gamma*log(tiny*gamma/const)

  modeDist <- hyperbMode(ThetaStand)
  ## Determine break points, based on size of derivative
  xDeriv <- seq(xSmall, modeDist, length.out = 101)
  derivVals <- ddhyperb(xDeriv, ThetaStand)
  maxDeriv <- max(derivVals)
  minDeriv <- min(derivVals)
  breakSize <- deriv*maxDeriv
  breakFun <- function(x){ 
    ddhyperb(x, ThetaStand, KNu) - breakSize
  }
  if ((maxDeriv < breakSize)||(derivVals[1] > breakSize)){
    lowBreak <- xSmall
  }else{
    whichMaxDeriv <- which.max(derivVals)
    lowBreak <- uniroot(breakFun, c(xSmall,xDeriv[whichMaxDeriv]))$root
  }

  xDeriv <- seq(modeDist, xLarge, length.out = 101)
  derivVals <- -ddhyperb(xDeriv, ThetaStand)
  maxDeriv <- max(derivVals)
  minDeriv <- min(derivVals)
  breakSize <- deriv*maxDeriv
  breakFun <- function(x){ 
    - ddhyperb(x, ThetaStand, KNu) - breakSize
  }
  if ((maxDeriv < breakSize)||(derivVals[101] > breakSize)){
    highBreak <- xLarge
  }else{
    whichMaxDeriv <- which.max(derivVals)
    highBreak <- uniroot(breakFun, c(xDeriv[whichMaxDeriv], xLarge))$root
  }
  breaks <-
    Theta[3]*c(xTiny,xSmall,lowBreak,highBreak,xLarge,xHuge,modeDist) +
      Theta[4]
  breaks <- list(xTiny = breaks[1], xSmall =  breaks[2], lowBreak =  breaks[3],
       highBreak =  breaks[4], xLarge =  breaks[5], xHuge =  breaks[6],
       modeDist =  breaks[7])
  return(breaks)
} ## End of hyperbBreaks()



