# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












qeunif <- function(p, min = 0, max = 1, Maxit.nr = 10, Tol.nr = 1.0e-6,
                   lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)   # 20150102 KaiH
  
  if (lower.tail) {
    if (log.arg)
      p <- exp(p)
  } else {
    p <- if (log.arg) -expm1(p) else 1 - p
  }


  ppp <- p
  vsmallno <- sqrt(.Machine$double.eps)
   smallno <- 0.10
  if (any(min >= max))
    stop("argument 'min' has values greater or equal ",
         "to argument 'max'")
  if (!is.Numeric( Tol.nr, length.arg = 1, positive = TRUE) ||
      Tol.nr > 0.10)
    stop("argument 'Tol.nr' is not a single positive value, ",
         "or is too large")
  nrok <- ppp >= vsmallno & ppp <= 1.0 - vsmallno & is.finite(ppp)

  eee <- qbeta(ppp, shape1 = 3, shape2 = 3)
  eee[ppp <        smallno] <- sqrt(ppp[ppp <  smallno])
  eee[ppp > 1.0 -  smallno] <- 1.0 - sqrt(1.0 - ppp[ppp > 1.0 -  smallno])


  for (iii in 1:Maxit.nr) {
    realdiff <- (peunif(eee[nrok]) - ppp[nrok]) / deunif(eee[nrok])
    eee[nrok] <- eee[nrok] - realdiff
    if (all(abs(realdiff) / (1.0 + abs(realdiff)) < Tol.nr )) break
    if (iii == Maxit.nr) warning("did not converge")
  }

  if (max(abs(peunif(eee[nrok]) - ppp[nrok])) > Tol.nr)
    warning("did not converge on the second check")

  eee[ppp <       vsmallno] <-       sqrt(      ppp[ppp <       vsmallno])
  eee[ppp > 1.0 - vsmallno] <- 1.0 - sqrt(1.0 - ppp[ppp > 1.0 - vsmallno])
  eee[ppp == 0] <- 0
  eee[ppp == 1] <- 1
  eee[ppp <  0] <- NA
  eee[ppp >  1] <- NA
  min + eee * (max - min)
}



peunif <- function(q, min = 0, max = 1,
                   lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (any(min >= max))
    stop("argument 'min' has values greater or equal to argument 'max'")

  eee <- (q - min) / (max - min)

  if (lower.tail) {
    if (log.p) {
      Gofy <- -log1p((1/eee - 1)^2)
      Gofy[eee < 0] <- -Inf
      Gofy[eee > 1] <- 0.0
    } else {
      Gofy <- eee^2 / (exp(2*log1p(-eee)) + eee^2)  # KaiH
      Gofy <- 1 / (1 + (1/eee - 1)^2)
      Gofy[eee < 0] <- 0.0
      Gofy[eee > 1] <- 1.0
    }
  } else {
    if (log.p) {
      Gofy <- 2*log1p(-eee) - log(exp(2*log1p(-eee)) + eee^2)
      Gofy[eee < 0] <- 0.0
      Gofy[eee > 1] <- -Inf
    } else {
      Gofy <- exp(2*log1p(-eee)) / (exp(2*log1p(-eee)) + eee^2)
      Gofy[eee < 0] <- 1
      Gofy[eee > 1] <- 0
    }
  }
  Gofy
}



deunif <- function(x, min = 0, max = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  if (any(min >= max))
    stop("argument 'min' has values greater or equal to argument 'max'")

  eee <- (x - min) / (max - min)

  if (log.arg) {
    ans <- log(2) + log(eee) + log1p(-eee) -
           2.0 * log(2*eee*(1-eee) - 1) - log(max - min)
    ans[eee <= 0.0] <- log(0.0)
    ans[eee >= 1.0] <- log(0.0)
  } else {
    gunif <- function(y)
        as.numeric(y >= 0 & y <= 1) * 2*y*(1-y) / (2*y*(1-y) - 1)^2
    ans <- gunif(eee) / (max - min)
    ans[is.infinite(x)] <- 0  # 20141209 KaiH
  }
  ans
}




reunif <- function(n, min = 0, max = 1) {

  qeunif(runif(n), min = min, max = max)
}







qenorm <- function(p, mean = 0, sd = 1, Maxit.nr = 10, Tol.nr = 1.0e-6,
                   lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  ppp <- p
  if (!is.Numeric( Tol.nr, length.arg = 1, positive = TRUE) ||
      Tol.nr > 0.10)
    stop("argument 'Tol.nr' is not a single ",
         "positive value, or is too large")
  nrok <- is.finite(ppp)

  eee <-  qnorm(ppp, sd = 2/3, lower.tail = lower.tail, log.p = log.p) 



  gnorm <- function(y) dnorm(y) / (y * (1-2*pnorm(y)) - 2*dnorm(y))^2


  for (iii in 1:Maxit.nr) {
    if (lower.tail) {
      realdiff <- if (log.p) {
        ln.ppp <- ppp
        (penorm(eee[nrok]) - exp(ln.ppp[nrok])) / gnorm(eee[nrok])
      } else {
        (penorm(eee[nrok]) - ppp[nrok]) / gnorm(eee[nrok])
      }
    } else {
      realdiff <- if (log.p) {
        ln.ppp <- ppp
        (penorm(eee[nrok]) + expm1(ln.ppp[nrok])) / gnorm(eee[nrok])
      } else {
        (penorm(eee[nrok]) + expm1(log(ppp[nrok]))) / gnorm(eee[nrok])
      }
    }
    eee[nrok] <- eee[nrok] - realdiff
    if (all(abs(realdiff) / (1.0 + abs(realdiff)) < Tol.nr ))
      break
    if (iii == Maxit.nr)
      warning("did not converge")
  }



  if (max(abs(penorm(eee[nrok]) - ppp[nrok])) > Tol.nr)
    warning("did not converge on the second check")


  if (lower.tail) {
    if (log.p) {
      eee[ln.ppp > 0] <- NaN
    } else {
      eee[ppp == 0] <- -Inf
      eee[ppp == 1] <-  Inf
      eee[ppp <  0] <- NaN
      eee[ppp >  1] <- NaN
    }
  } else {
    if (log.p) {
      eee[ln.ppp > 0] <- NaN
    } else { 
      eee[ppp == 0] <- Inf
      eee[ppp == 1] <- -Inf
      eee[ppp <  0] <- NaN
      eee[ppp >  1] <- NaN
    }
  }
  eee * ifelse(sd >= 0, sd, NaN) + mean
}



penorm <- function(q, mean = 0, sd = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  eee <- (q - mean) / sd
  tmp1 <- -dnorm(eee) - eee * pnorm(eee)

  if (lower.tail) {
    if (log.p) {
      Gofy <- log(tmp1 / (2 * tmp1 + eee))
      Gofy[eee <= -Inf] <- -Inf
      Gofy[eee >=  Inf] <- 0
    } else {
      Gofy <- tmp1 / (2 * tmp1 + eee)
      Gofy[eee <= -Inf] <- 0.0
      Gofy[eee >=  Inf] <- 1.0
    }
  } else {
    if (log.p) {
      Gofy <- log((tmp1 + eee) / (2 * tmp1 + eee))
      Gofy[eee <= -Inf] <- 0
      Gofy[eee >=  Inf] <- -Inf
    } else {
      Gofy <- (tmp1 + eee) / (2 * tmp1 + eee)
      Gofy[eee <= -Inf] <- 1
      Gofy[eee >=  Inf] <- 0
    }
  }
  Gofy
}



denorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  eee <- (x - mean) / sd
  if (log.arg) {
    ans <- dnorm(eee, log = TRUE) -
           2.0 * log(eee * (1-2*pnorm(eee)) - 2*dnorm(eee)) - log(sd)
  } else {
    gnorm <- function(y) dnorm(y) / (y * (1-2*pnorm(y)) - 2*dnorm(y))^2
    ans <- gnorm(eee) / sd
    ans[sd  <=  0.0] <- NaN
  }
  ans
}




renorm <- function(n, mean = 0, sd = 1) {

  qenorm(runif(n), mean = mean, sd = sd)
}







qeexp <- function(p, rate = 1, Maxit.nr = 10, Tol.nr = 1.0e-6,
                  lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)   # 20150102 KaiH
  
  if (lower.tail) {
    if (log.arg) p <- exp(p)
  } else {
    p <- if (log.arg) -expm1(p) else 1 - p
  }
  
  ppp <- p
  vsmallno <- sqrt(.Machine$double.eps)
  if (!is.Numeric( Tol.nr, length.arg = 1, positive = TRUE) ||
      Tol.nr > 0.10)
    stop("argument 'Tol.nr' is not a single positive value, or ",
         "is too large")
  nrok <- ppp >= vsmallno & is.finite(ppp)


  eee <- qf(1.0 * ppp, df1 =  4.0, df2 = 44)
  if ( any(rangex <- ppp < 0.8) )
      eee[rangex] <- qrayleigh(ppp[rangex], scale =  0.8)


  eee[ppp <       vsmallno] <- sqrt(ppp[ppp < vsmallno])


  for (iii in 1:Maxit.nr) {
    realdiff <- (peexp(eee[nrok]) - ppp[nrok]) / deexp(eee[nrok])
    eee[nrok] <- eee[nrok] - realdiff
    if (all(abs(realdiff) / (1.0 + abs(realdiff)) < Tol.nr )) break
    if (iii == Maxit.nr) warning("did not converge")
  }

  if (max(abs(peexp(eee[nrok]) - ppp[nrok])) > Tol.nr)
    warning("did not converge on the second check")

  eee[ppp < vsmallno] <- sqrt(ppp[ppp < vsmallno])
  eee[ppp == 0] <- 0
  eee[ppp == 1] <- Inf
  eee[ppp <  0] <- NaN
  eee[ppp >  1] <- NaN
  eee / rate
}



peexp <- function(q, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  eee <- q * rate
  tmp1 <- -expm1(-eee) - eee

  if (lower.tail) {
    if (log.p) {
      Gofy <- log(-tmp1) - log(expm1(-eee) + exp(-eee) + eee)
      Gofy[eee < 0 ] <- -Inf
      Gofy[eee == Inf] <- 0
    } else {
      Gofy <- tmp1 / (-expm1(-eee) - exp(-eee) - eee)
      Gofy[eee < 0] <- 0
      Gofy[eee == Inf] <- 1
    }
  } else {
    if (log.p) {
      Gofy <- -eee - log(expm1(-eee) + exp(-eee) + eee)
      Gofy[eee < 0] <- 0
      Gofy[eee == Inf] <- -Inf
    } else {
      Gofy <- exp(-eee)/(expm1(-eee) +exp(-eee) +eee)
      Gofy[eee < 0] <- 1
      Gofy[eee == Inf] <- 0
    }
  }
  Gofy
}



deexp <- function(x, rate = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  if (any(rate <= 0))
    stop("argument 'rate' must have positive values")

  eee <- x * rate

  if (log.arg) {
    ans <- log(eee) - eee + 2.0 * log((1-x) - 2 * exp(-x)) + log(rate)
    ans[is.infinite(x)] <- log(0)
  } else {
    gexp <- function(y)
      as.numeric(y >= 0) * y * exp(-y) / ((1-y) - 2 * exp(-y))^2
    ans <- gexp(eee) * rate
    ans[rate <=  0.0] <- NaN
    ans[is.infinite(x)] <- 0
  }
  ans
}



reexp <- function(n, rate = 1) {
  qeexp(runif(n), rate = rate)
}




dsc.t2 <- function(x, location = 0, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  zedd <- (x - location) / scale
  zedd[scale <= 0] <- NaN

  if (log.arg) {
    log(0.25) - 1.5 * log1p((zedd / 2)^2) - log(scale)
  } else {
    2 / (scale * (4 + zedd^2)^1.5)
  }
}





psc.t2 <- function(q, location = 0, scale = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  zedd <- (q - location) / scale
  zedd[scale <= 0] <- NaN

  if (lower.tail) {
    if (log.p) {
      ans <- log(0.5) + log1p(zedd / sqrt(4 + zedd^2))
      ans[q == -Inf] <- log(0)
      ans[q == Inf] <- log(1) 
      } else {
      ans <- 0.5 * (1 + zedd / sqrt(4 + zedd^2))
      ans[q == -Inf] <- 0 
      ans[q == Inf] <- 1
      }
  } else {
    if (log.p) {
      ans <- log(0.5) + log1p(-zedd / sqrt(4 + zedd^2))
      ans[q == -Inf] <- log(1)
      ans[q == Inf] <- log(0)
    } else {
      ans <- 0.5 * exp(log1p(-zedd / sqrt(4 + zedd^2)))
      ans[q == -Inf] <- 1 
      ans[q == Inf] <- 0 
    }
  }
  ans 
}






qsc.t2 <- function(p, location = 0, scale = 1,
                   lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- exp(0.5*(ln.p - log(-expm1(ln.p)))) -
             exp(0.5*(log(-expm1(ln.p)) - ln.p))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- exp(0.5*(log(p) - log1p(-p))) -
             exp(0.5*(log1p(-p) - log(p)))
      ans[p < 0] <- NaN
      ans[p == 0] <- -Inf
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- exp(0.5*(log(-expm1(ln.p)) - ln.p)) -
             exp(0.5*(ln.p - log(-expm1(ln.p))))
      ans[ln.p > 0] <- NaN
      ans
    } else { 
      ans <- exp(0.5*(log1p(-p) - log(p))) -
             exp(0.5*(log(p) - log1p(-p)))
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- -Inf
      ans[p > 1] <- NaN
    }
  }

  answer <- ans * scale + location
  answer[scale <= 0] <- NaN
  answer
}




rsc.t2 <- function(n, location = 0, scale = 1) {
  answer <- qsc.t2(runif(n)) * scale + location
  answer[scale <= 0] <- NaN
  answer
}







 sc.studentt2 <- function(percentile = 50,
                     llocation = "identitylink", lscale = "loge",
                     ilocation = NULL,   iscale = NULL,
                     imethod = 1,
                     zero = "scale") {

 

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  ilocat <- ilocation


  if (length(ilocat) &&
     (!is.Numeric(ilocat, length.arg = 1, positive = TRUE)))
      stop("bad input for argument 'ilocation'")
  if (length(iscale) && !is.Numeric(iscale))
    stop("bad input for argument 'iscale'")


  if (!is.Numeric(percentile, positive = TRUE) ||
      any(percentile >= 100))
    stop("bad input for argument 'percentile'")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Scaled Student t distribution with 2 degrees of freedom\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat, tag = FALSE), ", ",
            namesof("scale",    lscale, earg = escale, tag = FALSE), "\n\n",
            "Mean:     location\n",
            "Variance: infinite"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale"),
         llocation = .llocation ,
         lscale    = .lscale ,
         zero = .zero )
  }, list( .zero = zero, .llocation = llocation, .lscale = lscale ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
        namesof("location", .llocat, earg = .elocat, tag = FALSE),
        namesof("scale",    .lscale, earg = .escale, tag = FALSE))


    if (!length(etastart)) {
      locat.init <- if ( .imethod == 2) {
        weighted.mean(y, w)
      } else {
        median(y)
      }
      Scale.init <- if (length( .iscale )) .iscale else
        diff(quantile(y, prob = c(0.25, 0.75))) / (2 * 1.155) + 1.0e-5
      locat.init <- rep(locat.init, length = length(y))
      Scale.init <- rep(Scale.init, length = length(y))
      etastart <- cbind(theta2eta(locat.init, .llocat, earg = .elocat),
                        theta2eta(Scale.init, .lscale, earg = .escale))
    }
  }), list( .llocat = llocat, .lscale = lscale,
            .ilocat = ilocat, .iscale = iscale,
            .elocat = elocat, .escale = escale,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    Perce <- .percentile
    locat <- eta2theta(eta[, 1], link = .llocat, earg = .elocat)
    Scale <- eta2theta(eta[, 2], link = .lscale, earg = .escale)
    answer <- matrix(locat, nrow(eta), length(Perce))
    for (ii in 1:length(Perce))
      answer[, ii] <- qsc.t2(Perce[ii] / 100, loc = locat, sc = Scale)
    dimnames(answer) <- list(dimnames(eta)[[1]],
                             paste(as.character(Perce), "%", sep = ""))
    answer
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .percentile = percentile ))),
  last = eval(substitute(expression({
    misc$link <-    c("location" = .llocat, "scale" = .lscale)

    misc$earg <- list("location" = .elocat, "scale" = .escale)

    misc$expected <- TRUE
    misc$percentile <- .percentile
    misc$imethod <- .imethod
    misc$multipleResponses <- FALSE

      ncoly <- ncol(y)
      for (ii in 1:length( .percentile )) {
        y.use <- if (ncoly > 1) y[, ii] else y
        mu <- cbind(mu)
        extra$percentile[ii] <- 100 * weighted.mean(y.use <= mu[, ii], w)
      }
      names(extra$percentile) <- colnames(mu)
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .imethod = imethod, .percentile = percentile ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    locat <- eta2theta(eta[, 1], link = .llocat, earg = .elocat)
    Scale <- eta2theta(eta[, 2], link = .lscale, earg = .escale)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dsc.t2(x = y, location = locat, scale = Scale,
                                 log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale ))),
  vfamily = c("sc.studentt2"),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta[, 1], link = .llocat, earg = .elocat)
    Scale <- eta2theta(eta[, 2], link = .lscale, earg = .escale)

    dlocat.deta <- dtheta.deta(locat, link = .llocat, earg = .elocat)
    dscale.deta <- dtheta.deta(Scale, link = .lscale, earg = .escale)

    zedd <- (y - locat) / Scale

    dl.dlocat <- 3 * zedd   / (Scale * (4 + zedd^2))
    dl.dscale <- 3 * zedd^2 / (Scale * (4 + zedd^2)) - 1 / Scale

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta)
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  weight = eval(substitute(expression({
    ned2l.dlocat2 <- 0.3 / Scale^2
    ned2l.dscale2 <- 2.0 / (3 * Scale^2)

    wz <- matrix(-10, n, M)  # Diagonal EIM
    wz[, iam(1, 1, M = M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2, M = M)] <- ned2l.dscale2 * dscale.deta^2

    c(w) * wz
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))))
}










