dvg <- function (x, vgC = 0, sigma = 1, theta = 0, nu = 1,
  param = c(vgC,sigma,theta,nu), log = FALSE,
  tolerance = .Machine$double.eps ^ 0.5, ...) {

  ## check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  if (log == TRUE) {
    stop ("This function is not yet implemented")
  } else {
    if (abs(nu - 2) < tolerance) {
      vgDens <- ifelse(abs(x - vgC) < tolerance, Inf,
        ((2*exp(theta*(x - vgC)/sigma^2))/(nu^(1/nu)*sqrt(2*pi)*
          sigma*gamma(1/nu)))*((abs(x - vgC)/sqrt(2*(sigma^2)/nu + theta^2))^
          (1/nu - 1/2))*besselK(x = ((1/sigma^2)*abs(x - vgC)*
          sqrt((2*sigma^2/nu) + theta^2)), nu = (1/nu - 1/2)))
    } else {
      if(nu < 2) {
        if (nu^(1/nu)- 0 < tolerance) {
          vgDens <- NA
        } else {
          vgDens <- ifelse(abs(x - vgC) < tolerance, gamma(1/nu - 1/2)/(sigma*
            sqrt(2*pi)*nu^(1/nu)*gamma(1/nu))*((2*sigma^2/
            (2*sigma^2/nu + theta^2))^(1/nu - 1/2)),
            ((2*exp(theta*(x - vgC)/sigma^2))/(nu^(1/nu)*sqrt(2*pi)*
            sigma*gamma(1/nu)))*((abs(x - vgC)/sqrt(2*(sigma^2)/nu + theta^2))^
            (1/nu - 1/2))*besselK(x = ((1/sigma^2)*abs(x - vgC)*
            sqrt((2*sigma^2/nu) + theta^2)), nu = (1/nu - 1/2)))
        }
      }
      if (nu > 2) {
        vgDens <- ifelse(abs(x - vgC) < tolerance, Inf,
          ((2*exp(theta*(x - vgC)/sigma^2))/(nu^(1/nu)*sqrt(2*pi)*
          sigma*gamma(1/nu)))*((abs(x - vgC)/sqrt(2*(sigma^2)/nu + theta^2))^
          (1/nu - 1/2))*besselK(x = ((1/sigma^2)*abs(x - vgC)*
          sqrt((2*sigma^2/nu) + theta^2)), nu = (1/nu - 1/2)))
      }
    }
    vgDens <- ifelse(is.nan(vgDens), 0, vgDens)
  }

  return(vgDens)
}

pvg <- function (q, vgC = 0, sigma = 1, theta = 0, nu = 1,
                 param = c(vgC,sigma,theta,nu), lower.tail = TRUE,
                 log.p = FALSE, small = 10^(-6), tiny = 10^(-10),
                 deriv = 0.3, subdivisions = 100, accuracy = FALSE, ...) {

  ## check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  if (lower.tail == TRUE) {
    if (log.p == FALSE){
      bks <- vgBreaks(param = param, small = small, tiny = tiny,
                      deriv = deriv, ...)
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
      if (length(qLow) > 0)
        qLessEqMode <- qLessEqMode[qLessEqMode > max(qLow)]
      if (length(qHigh) > 0)
        qGreatMode <- qGreatMode[qGreatMode < min(qHigh)]
      if (length(qSmall) > 0)
        qLow <- qLow[qLow > max(qSmall)]
      if (length(qLarge) > 0)
        qHigh <- qHigh[qHigh < min(qLarge)]
      if (length(qTiny) > 0)
        qSmall <- qSmall[qSmall > max(qTiny)]
      if (length(qHuge) > 0)
        qLarge <- qLarge[qLarge < min(qHuge)]
      intFun <- rep(NA, length(q))
      if (length(qTiny) > 0)
        intFun[qTiny] <- 0
      if (length(qHuge) > 0)
        intFun[qHuge] <- 1
      intErr <- rep(NA, length(q))
      if (length(qTiny) > 0)
        intErr[qTiny] <- tiny
      if (length(qHuge) > 0)
        intErr[qHuge] <- tiny
      dvgInt <- function (q) {
        dvg(q, param = param, log = FALSE)
      }
      resSmall <- safeIntegrate(dvgInt, xTiny, xSmall,
                                subdivisions, ...)
      resLarge <- safeIntegrate(dvgInt, xLarge, xHuge,
                                subdivisions, ...)
      intSmall <- resSmall$value
      intLarge <- resLarge$value
      errSmall <- tiny + resSmall$abs.error
      errLarge <- tiny + resLarge$abs.error
      resLow <- safeIntegrate(dvgInt, xSmall, lowBreak, subdivisions, ...)
      resHigh <- safeIntegrate(dvgInt, highBreak, xLarge, subdivisions, ...)
      intLow <- intSmall + resLow$value
      intHigh <- intLarge + resHigh$value
      errLow <- errSmall + resLow$abs.error
      errHigh <- errLarge + resHigh$abs.error
      for (i in qSmall) {
        intRes <- safeIntegrate(dvgInt, xTiny, qSort[i], subdivisions, ...)
        intFun[i] <- intRes$value
        intErr[i] <- intRes$abs.error + tiny
      }
      for (i in qLarge) {
        intRes <- safeIntegrate(dvgInt, qSort[i], xHuge, subdivisions, ...)
        intFun[i] <- 1 - intRes$value
        intErr[i] <- intRes$abs.error + tiny
      }
      for (i in qLow) {
        intRes <- safeIntegrate(dvgInt, xSmall, qSort[i], subdivisions, ...)
        intFun[i] <- intRes$value + intSmall
        intErr[i] <- intRes$abs.error + errSmall
      }
      for (i in qHigh) {
        intRes <- safeIntegrate(dvgInt, qSort[i], xLarge, subdivisions, ...)
        intFun[i] <- 1 - intRes$value - intLarge
        intErr[i] <- intRes$abs.error + errLarge
      }
      for (i in qLessEqMode) {
        intRes <- safeIntegrate(dvgInt, lowBreak, qSort[i], subdivisions, ...)
        intFun[i] <- intRes$value + intLow
        intErr[i] <- intRes$abs.error + errLow
      }
      for (i in qGreatMode) {
        intRes <- safeIntegrate(dvgInt, qSort[i], highBreak, subdivisions, ...)
        intFun[i] <- 1 - intRes$value - intHigh
        intErr[i] <- intRes$abs.error + errLarge
      }
      if (!accuracy) {
        return(intFun[rank(q)])
      } else {
        return(list(value = intFun[rank(q)], error = intErr[rank(q)]))
      }
    }

    if (log.p == TRUE) {
      stop("This function is not yet implemented")
    }
  }
  if (lower.tail == FALSE) {
    stop("This function is not yet implemented")
    if (log.p == FALSE) {
    }
    if (log.p == TRUE) {
    }
  }
}

qvg <- function (p, vgC = 0, sigma = 1, theta = 0, nu = 1,
                 param = c(vgC,sigma,theta,nu), lower.tail = TRUE,
                 log.p = FALSE, small = 10^(-6), tiny = 10^(-10),
                 deriv = 0.3, nInterpol = 100, subdivisions = 100, ...) {

  ## check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  if (lower.tail == TRUE) {
    if (log.p == FALSE){
      bks <- vgBreaks(param = param, small = small, tiny = tiny,
                      deriv = deriv, ...)
      xTiny <- bks$xTiny
      xSmall <- bks$xSmall
      lowBreak <- bks$lowBreak
      highBreak <- bks$highBreak
      xLarge <- bks$xLarge
      xHuge <- bks$xHuge
      modeDist <- bks$modeDist
      yTiny <- pvg(q = xTiny, param = param)
      ySmall <- pvg(q = xSmall, param = param)
      yLowBreak <- pvg(q = lowBreak, param = param)
      yHighBreak <- pvg(q= highBreak, param = param)
      yLarge <- pvg(q = xLarge, param = param)
      yHuge <- pvg(q = xHuge, param = param)
      yModeDist <- pvg(q = modeDist, param = param)
      pSort <- sort(p)
      pSmall <- which(pSort < pvg(q = xSmall, param = param))
      pTiny <- which(pSort < pvg(q = xTiny, param = param))
      pLarge <- which(pSort > pvg(q = xLarge, param = param))
      pHuge <- which(pSort > pvg(q = xHuge, param = param))
      pLow <- which(pSort < pvg(q = lowBreak, param = param))
      pHigh <- which(pSort > pvg(q = highBreak, param = param))
      pLessEqMode <- which(pSort <= pvg(q = modeDist, param = param))
      pGreatMode <- which(pSort > pvg(q = modeDist, param = param))
      if (length(pLow) > 0)
        pLessEqMode <- pLessEqMode[pLessEqMode > max(pLow)]
      if (length(pHigh) > 0)
        pGreatMode <- pGreatMode[pGreatMode < min(pHigh)]
      if (length(pSmall) > 0)
        pLow <- pLow[pLow > max(pSmall)]
      if (length(pLarge) > 0)
        pHigh <- pHigh[pHigh < min(pLarge)]
      if (length(pTiny) > 0)
        pSmall <- pSmall[pSmall > max(pTiny)]
      if (length(pHuge) > 0)
        pLarge <- pLarge[pLarge < min(pHuge)]
      qSort <- rep(NA, length(pSort))
      if (length(pTiny) > 0)
        qSort[pTiny] <- -Inf
      if (length(pHuge) > 0)
        qSort[pHuge] <- Inf
      if (length(pTiny) > 0) {
        for (i in pTiny) {
          zeroFun <- function(x) {
            pvg(q = x, param = param) - pSort[i]
          }
          interval <- c(xTiny - (xSmall - xTiny), xTiny)
          while (zeroFun(interval[1]) * zeroFun(interval[2]) > 0) {
            interval[1] <- interval[1] - (xSmall - xTiny)
          }
          qSort[i] <- uniroot(zeroFun, interval)$root
        }
      }

      if (length(pSmall) > 0) {
        xValues <- seq(xTiny, xSmall, length = nInterpol)
        pvgValues <- pvg(q = xValues, param = param, small = small,
                         tiny = tiny, deriv = deriv,
                         subdivisions = subdivisions, accuracy = FALSE)
        pvgSpline <- splinefun(xValues, pvgValues)
        for (i in pSmall) {
          zeroFun <- function(x) {
            pvgSpline(x) - pSort[i]
          }
          if (zeroFun(xTiny) >= 0) {
            qSort[i] <- xTiny
          } else {
            if (zeroFun(xSmall) <= 0) {
              qSort[i] <- xSmall
            } else {
              qSort[i] <- uniroot(zeroFun,
                                  interval = c(xTiny,xSmall),...)$root
            }
          }
        }
      }

      if (length(pLow) > 0) {
        xValues <- seq(xSmall, lowBreak, length = nInterpol)
        pvgValues <- pvg(q = xValues, param = param, small = small,
                         tiny = tiny, deriv = deriv,
                         subdivisions = subdivisions, accuracy = FALSE)
        pvgSpline <- splinefun(xValues, pvgValues)
        for (i in pLow) {
          zeroFun <- function (x) {
            pvgSpline(x) - pSort[i]
          }
          if (zeroFun(xSmall) >= 0) {
            qSort[i] <- xSmall
          } else {
            if (zeroFun(lowBreak) <= 0) {
              qSort[i] <- lowBreak
            } else {
              qSort[i] <- uniroot(zeroFun,
                                  interval = c(xSmall,lowBreak), ...)$root
            }
          }
        }
      }

      if (length(pLessEqMode) > 0) {
          xValues <- seq(lowBreak, modeDist, length = nInterpol)
          pvgValues <- pvg(q = xValues, param = param, small = small,
                           tiny = tiny, deriv = deriv,
                           subdivisions = subdivisions, accuracy = FALSE)
          pvgSpline <- splinefun(xValues, pvgValues)
          for (i in pLessEqMode) {
              zeroFun <- function(x) {
                  pvgSpline(x) - pSort[i]
              }
              if (zeroFun(lowBreak) >= 0) {
                  qSort[i] <- lowBreak
              }
              else {
                  if (zeroFun(modeDist) <= 0) {
                    qSort[i] <- modeDist
                  }
                  else {
                    qSort[i] <- uniroot(zeroFun,
                                        interval = c(lowBreak,modeDist),
                                        ...)$root
                  }
              }
          }
      }

      if (length(pGreatMode) > 0) {
        xValues <- seq(modeDist, highBreak, length = nInterpol)
        pvgValues <- pvg(q = xValues, param = param, small = small,
                         tiny = tiny, deriv = deriv,
                         subdivisions = subdivisions, accuracy = FALSE)
        pvgSpline <- splinefun(xValues, pvgValues)
        for (i in pGreatMode) {
          zeroFun <- function (x) {
            pvgSpline(x) - pSort[i]
          }
          if (zeroFun(modeDist) >= 0) {
            qSort[i] <- modeDist
          } else {
            if (zeroFun(highBreak) <= 0) {
              qSort[i] <- highBreak
            } else {
              qSort[i] <- uniroot(zeroFun,
                                  interval = c(modeDist,highBreak), ...)$root
            }
          }
        }
      }

      if (length(pHigh) > 0) {
        xValues <- seq(highBreak, xLarge, length = nInterpol)
        pvgValues <- pvg(q = xValues, param = param, small = small,
                         tiny = tiny, deriv = deriv,
                         subdivisions = subdivisions, accuracy = FALSE)
        pvgSpline <- splinefun(xValues, pvgValues)
        for (i in pHigh) {
          zeroFun <- function (x) {
            pvgSpline(x) - pSort[i]
          }
          if (zeroFun(highBreak) >= 0) {
            qSort[i] <- highBreak
          } else {
            if (zeroFun(xLarge) <= 0) {
              qSort[i] <- xLarge
            } else {
              qSort[i] <- uniroot(zeroFun,
                                  interval = c(highBreak,xLarge), ...)$root
            }
          }
        }
      }

      if (length(pLarge) > 0) {
        xValues <- seq(xLarge, xHuge, length = nInterpol)
        pvgValues <- pvg(q = xValues, param = param, small = small,
                         tiny = tiny, deriv = deriv,
                         subdivisions = subdivisions, accuracy = FALSE)
        pvgSpline <- splinefun(xValues, pvgValues)
        for (i in pLarge) {
          zeroFun <- function (x) {
            pvgSpline(x) - pSort[i]
          }
          if (zeroFun(xLarge) >= 0) {
            qSort[i] <- xLarge
          } else {
            if (zeroFun(xHuge) <= 0) {
              qSort[i] <- xHuge
            } else {
              qSort[i] <- uniroot(zeroFun,
                                  interval = c(xLarge,xHuge), ...)$root
            }
          }
        }
      }

      if (length(pHuge) > 0) {
        for (i in pHuge) {
          zeroFun <- function (x) {
            pvg(q = x, param = param) - pSort[i]
          }
          interval <- c(xHuge,xHuge + (xHuge - xLarge))
          while (zeroFun(interval[1]) * zeroFun(interval[2]) > 0) {
            interval[1] <- interval[1] + (xHuge - xLarge)
          }
            qSort[i] <- uniroot(zeroFun, interval)$root
        }
      }
      return(qSort[rank(p)])
    }
  if (log.p == TRUE) {
      stop("This function is not yet implemented")
    }
  }
  if (lower.tail == FALSE) {
    stop("This function is not yet implemented")
    if (log.p == FALSE) {
    }
    if (log.p == TRUE) {
    }
  }
}

rvg <- function (n, vgC = 0, sigma = 1, theta = 0, nu = 1,
                 param = c(vgC,sigma,theta,nu)) {

  ## check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  kkp1Param <- vgChangePars(1, 2, param = param)
  mu <- kkp1Param[3]
  kkp1Sigma <- kkp1Param[2]
  tau <- kkp1Param[4]
  kkp1Theta <- kkp1Param[1]
  kkp2Param <- vgChangePars(1, 3, param = param)
  kappa <- kkp2Param[3]
  rgamma1 <- rgamma(n, shape = tau, rate = 1)
  rgamma2 <- rgamma(n, shape = tau, rate = 1)

  X <- kkp1Theta + (kkp1Sigma/sqrt(2))*((1/kappa)*rgamma1 - kappa*rgamma2)
  return(X)
}

ddvg <- function (x,  vgC = 0, sigma = 1, theta = 0, nu = 1,
                  param = c(vgC,sigma,theta,nu), log = FALSE,
                  tolerance = .Machine$double.eps ^ 0.5, ...) {

  ## check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  if (log == TRUE) {
    stop ("This function is not yet implemented")
  } else {
    if (abs(nu - 2) < tolerance) {
      ddvg <- ifelse(abs(x - vgC) < tolerance, NA,
         exp(theta*(x - vgC)/sigma^2)*
          2^(1/2)*abs(x - vgC)^(-1/2*(-2 + nu)/nu)*
                     ((2*sigma^2 + theta^2*nu)/nu)^
          (1/4*(-2 + nu)/nu)*(theta*besselK(x = 1/sigma^2*abs(x - vgC)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2), nu = 1/2*(-2 + nu)/nu)-
          ((x - vgC)/abs(x - vgC))*besselK(x = 1/sigma^2*abs(x - vgC)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2), nu = 1/2*(-2 + 3*nu)/nu)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2))*nu^(-1/nu)/pi^(1/2)/sigma^3/
          gamma(1/nu))
    } else {
      if (nu < 2) {
        ddvg <- ifelse(abs(x - vgC) < tolerance, 0,
          exp(theta*(x - vgC)/sigma^2)*
          2^(1/2)*abs(x - vgC)^(-1/2*(-2 + nu)/nu)*
                       ((2*sigma^2 + theta^2*nu)/nu)^
          (1/4*(-2 + nu)/nu)*(theta*besselK(x = 1/sigma^2*abs(x - vgC)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2), nu = 1/2*(-2 + nu)/nu)-
          ((x - vgC)/abs(x - vgC))*besselK(x = 1/sigma^2*abs(x - vgC)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2), nu = 1/2*(-2 + 3*nu)/nu)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2))*nu^(-1/nu)/pi^(1/2)/sigma^3/
          gamma(1/nu))
      }
      if (nu > 2) {
        ddvg <- ifelse(abs(x - vgC) < tolerance, NA,
          exp(theta*(x - vgC)/sigma^2)*
          2^(1/2)*abs(x - vgC)^(-1/2*(-2 + nu)/nu)*
                       ((2*sigma^2 + theta^2*nu)/nu)^
          (1/4*(-2 + nu)/nu)*(theta*besselK(x = 1/sigma^2*abs(x - vgC)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2), nu = 1/2*(-2 + nu)/nu)-
          ((x - vgC)/abs(x - vgC))*besselK(x = 1/sigma^2*abs(x - vgC)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2), nu = 1/2*(-2 + 3*nu)/nu)*
          ((2*sigma^2 + theta^2*nu)/nu)^(1/2))*nu^(-1/nu)/pi^(1/2)/sigma^3/
          gamma(1/nu))
      }
    }
  }
  return(ddvg)
}

vgBreaks <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
                      param = c(vgC,sigma,theta,nu), small = 10^(-6),
                      tiny = 10^(-10), deriv = 0.3, ...) {

  #check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  if (nu < 2) {
  xTiny <- vgCalcRange(param = param, tol = tiny, density = TRUE)[1]
  xSmall <- vgCalcRange(param = param, tol = small, density = TRUE)[1]
  xLarge <- vgCalcRange(param = param, tol = small, density = TRUE)[2]
  xHuge <- vgCalcRange(param = param, tol = tiny, density = TRUE)[2]
  modeDist <- vgMode(param = param)
  xDeriv <- seq(xSmall, modeDist, length.out = 101)
  derivVals <- ddvg(x= xDeriv, param = param)
  maxDeriv <- max(derivVals)
  minDeriv <- min(derivVals)
  breakSize <- deriv*maxDeriv
  breakFun <- function (x) {
    ddvg(x, param = param) - breakSize
  }
  if ((maxDeriv < breakSize) || (derivVals[1] > breakSize)) {
    lowBreak <- xSmall
  } else {
    whichMaxDeriv <- which.max(derivVals)
    lowBreak <- uniroot(breakFun, c(xSmall,xDeriv[whichMaxDeriv]))$root
  }
  xDeriv <- seq(modeDist, xLarge, length.out = 101)
  derivVals <- -ddvg(x = xDeriv, param = param)
  maxDeriv <- max(derivVals)
  minDeriv <- min(derivVals)
  breakSize <- deriv*maxDeriv
  breakFun <- function (x) {
    -ddvg(x, vgC, sigma, theta, nu) - breakSize
  }
  if ((maxDeriv < breakSize) || (derivVals[101] > breakSize)) {
    highBreak <- xLarge
  } else {
    whichMaxDeriv <- which.max(derivVals)
    highBreak <- uniroot(breakFun, c(xDeriv[whichMaxDeriv],xLarge))$root
  }
  breaks <-c(xTiny, xSmall, lowBreak, highBreak, xLarge, xHuge, modeDist)
  breaks <- list(xTiny = breaks[1], xSmall = breaks[2], lowBreak = breaks[3],
                 highBreak = breaks[4], xLarge = breaks[5],
                 xHuge = breaks[6], modeDist = breaks[7])
  }
  if (nu >= 2) {

    xTiny <- vgCalcRange(param = param, tol = tiny, density = TRUE)[1]
    xSmall <- vgCalcRange(param = param, tol = small, density = TRUE)[1]
    xLarge <- vgCalcRange(param = param, tol = small, density = TRUE)[2]
    xHuge <- vgCalcRange(param = param, tol = tiny, density = TRUE)[2]
    modeDist <- vgMode(param = param)
    xDeriv <- seq(xSmall, modeDist - 0.165, length.out = 101)
    derivVals <- ddvg(x = xDeriv, param = param)
    maxDeriv <- max(derivVals)
    minDeriv <- min(derivVals)
    breakSize <- deriv*maxDeriv
    breakFun <- function (x) {
      ddvg(x, vgC, sigma, theta, nu) - breakSize
    }
    if ((maxDeriv < breakSize) || (derivVals[1] > breakSize)) {
      lowBreak <- xSmall
    } else {
      whichMaxDeriv <- which.max(derivVals)
      lowBreak <- uniroot(breakFun, c(xSmall,xDeriv[whichMaxDeriv]))$root
    }
    xDeriv <- seq(modeDist + 0.165, xLarge, length.out = 101)
    derivVals <- -ddvg(x = xDeriv, param = param)
    maxDeriv <- max(derivVals)
    minDeriv <- min(derivVals)
    breakSize <- deriv*maxDeriv
    breakFun <- function (x) {
      -ddvg(x, param = param) - breakSize
    }
    if ((maxDeriv < breakSize) || (derivVals[101] > breakSize)) {
      highBreak <- xLarge
    } else {
      whichMaxDeriv <- which.max(derivVals)
      highBreak <- uniroot(breakFun, c(xDeriv[whichMaxDeriv],xLarge))$root
    }
    breaks <-c(xTiny, xSmall, lowBreak, highBreak, xLarge, xHuge, modeDist)
    breaks <- list(xTiny = breaks[1], xSmall = breaks[2], lowBreak = breaks[3],
                   highBreak = breaks[4], xLarge = breaks[5],
                   xHuge = breaks[6], modeDist = breaks[7])
    }

  return(breaks)
}
