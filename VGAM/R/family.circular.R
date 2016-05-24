# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.









dcard <- function(x, mu, rho, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(mu), length(rho))
  if (length(x)   != L) x   <- rep(x,   len = L)
  if (length(mu)  != L) mu  <- rep(mu,  len = L)
  if (length(rho) != L) rho <- rep(rho, len = L)

  logdensity <- rep(log(0), len = L)
  xok <- (x > 0) & (x < (2*pi))
  logdensity[xok] <- -log(2*pi) + log1p(2 * rho[xok] *
                      cos(x[xok]-mu[xok]))
  logdensity[mu  <=    0] <- NaN
  logdensity[mu  >= 2*pi] <- NaN
  logdensity[rho <= -0.5] <- NaN
  logdensity[rho >=  0.5] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



pcard <- function(q, mu, rho, lower.tail = TRUE, log.p = FALSE) {
  
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log((q + 2 * rho * (sin(q-mu) + sin(mu))) / (2*pi))
      ans[q <= 0 ] <- -Inf
      ans[q >= (2*pi)]  <- 0
    } else {
      ans <- (q + 2 * rho * (sin(q-mu) + sin(mu))) / (2*pi)
      ans[q <= 0] <- 0
      ans[q >= (2*pi)] <- 1
    }
  } else {
    if (log.p) {
      ans <- log1p(-(q + 2 * rho * (sin(q-mu) + sin(mu))) / (2*pi))
      ans[q <= 0] <- 0
      ans[q >= (2*pi)]  <- -Inf
    } else {
      ans <- (2*pi - q - 2 * rho * (sin(q-mu) + sin(mu))) / (2*pi)
      ans[q <= 0] <- 1
      ans[q >= (2*pi)] <- 0
    }
  } 
  ans[mu < 0 | mu > 2*pi] <- NaN  # A warning() may be a good idea here
  ans[abs(rho) > 0.5] <- NaN
  ans
}



qcard <- function(p, mu, rho, tolerance = 1.0e-7, maxits = 500,
                  lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(p) || any(p < 0) || any(p > 1))
    stop("'p' must be between 0 and 1")

  nn <- max(length(p), length(mu), length(rho))
  if (length(p)   != nn) p   <- rep(p,   len = nn)
  if (length(mu)  != nn) mu  <- rep(mu,  len = nn)
  if (length(rho) != nn) rho <- rep(rho, len = nn)


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      for (its in 1:maxits) {
        oldans <- 2 * pi * exp(ln.p)
        ans <- oldans - (oldans + 2 * rho * (sin(oldans-mu)+sin(mu)) -
                           2*pi*exp(ln.p)) / (1 + 2 * rho * cos(oldans - mu))
        index <- (ans < 0) | (ans > 2*pi)  # 20141216 KaiH  Remove ans == 0
        if (any(index)) {
          ans[index] <- runif (sum(index), 0, 2*pi)
        }
        if (max(abs(ans - oldans)) < tolerance)
          break
        if (its == maxits) {
          warning("did not converge")
          break
        }
        oldans <- ans
      }
    } else {
      for (its in 1:maxits) {
        oldans <- 2 * pi * p
        ans <- oldans - (oldans + 2 * rho * (sin(oldans-mu)+sin(mu)) -
               2*pi*p) / (1 + 2 * rho * cos(oldans - mu))
        index <- (ans < 0) | (ans > 2*pi)  # 20141216 KaiH  Remove ans == 0
        if (any(index)) {
          ans[index] <- runif(sum(index), 0, 2*pi)
        }
        if (max(abs(ans - oldans)) < tolerance)
          break
        if (its == maxits) {
          warning("did not converge")
          break
        }
        oldans <- ans
      }
    }
  } else {
    if (log.p) {
      ln.p <- p
      for (its in 1:maxits) {
        oldans <- - 2 * pi * expm1(ln.p)
        ans <- oldans - (oldans + 2 * rho * (sin(oldans-mu)+sin(mu)) +
               2*pi*expm1(ln.p)) / (1 + 2 * rho * cos(oldans - mu))
        index <- (ans < 0) | (ans > 2*pi)
        if (any(index)) {
          ans[index] <- runif (sum(index), 0, 2*pi)
        }
        if (max(abs(ans - oldans)) < tolerance)
          break
        if (its == maxits) {
          warning("did not converge")
          break
        }
        oldans <- ans 
       }
    } else { 
      for (its in 1:maxits) {
        oldans <- 2 * pi - 2 * pi * p
        ans <- oldans - (oldans + 2 * rho * (sin(oldans-mu)+sin(mu)) -
               2*pi + 2*pi*p) / (1 + 2 * rho * cos(oldans - mu))
        index <- (ans < 0) | (ans > 2*pi)  
        if (any(index)) {
          ans[index] <- runif (sum(index), 0, 2*pi)
        }
        if (max(abs(ans - oldans)) < tolerance)
          break
        if (its == maxits) {
          warning("did not converge")
          break
        }
        oldans <- ans
      }
    }
  }

  ans[mu < 0 | mu > 2*pi] <- NaN  # A warning() may be a good idea here
  ans[abs(rho) > 0.5] <- NaN
  ans
}



rcard <- function(n, mu, rho, ...) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  if (!is.Numeric(mu) || any(mu < 0) || any(mu > 2*pi))
    stop("argument 'mu' must be between 0 and 2*pi inclusive")
  if (!is.Numeric(rho) || max(abs(rho) > 0.5))
    stop("argument 'rho' must be between -0.5 and 0.5 inclusive")

  mu <- rep(mu, len = use.n)
  rho <- rep(rho, len = use.n)
  qcard(runif(use.n), mu = mu, rho = rho, ...)
}




cardioid.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 cardioid <- function(
     lmu  = extlogit(min = 0, max = 2*pi),
     lrho = extlogit(min = -0.5, max = 0.5),
     imu = NULL, irho = 0.3,
     nsimEIM = 100, zero = NULL) {

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")



  if (length(imu) && (!is.Numeric(imu, positive = TRUE) ||
      any(imu > 2*pi)))
    stop("bad input for argument 'imu'")
  if (!is.Numeric(irho) || max(abs(irho)) > 0.5)
    stop("bad input for argument 'irho'")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Cardioid distribution\n\n",
            "Links:    ",
            namesof("mu",  lmu,  earg = emu,  tag = FALSE), ", ", 
            namesof("rho", lrho, earg = erho, tag = FALSE), "\n",
            "Mean:     ",
            "pi + (rho/pi) *",
            "((2*pi-mu)*sin(2*pi-mu)+cos(2*pi-mu)-mu*sin(mu)-cos(mu))"),
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
         parameters.names = c("mu", "rho"),
         nsimEIM = .nsimEIM ,
         lmu  = .lmu  ,
         lrho = .lrho ,
         zero = .zero )
  }, list( .zero = zero, .lmu = lmu, .lrho = lrho,
           .nsimEIM = nsimEIM ))),


  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    if (any((y <= 0) | (y >=2*pi)))
      stop("the response must be in (0, 2*pi)")

    predictors.names <- c(
      namesof("mu",  .lmu  , earg = .emu  , tag = FALSE),
      namesof("rho", .lrho , earg = .erho , tag = FALSE))

    if (!length(etastart)) {
      rho.init <- rep(if (length( .irho )) .irho else 0.3, length = n)

      cardioid.Loglikfun <- function(mu, y, x, w, extraargs) {
        rho <- extraargs$irho
        sum(w * (-log(2*pi) + log1p(2*rho*cos(y-mu))))
      }
      mu.grid <- seq(0.1, 6.0, len = 19)
      mu.init <- if (length( .imu )) .imu else
          grid.search(mu.grid, objfun = cardioid.Loglikfun,
                      y = y,  x = x, w = w,
                      extraargs = list(irho = rho.init))
      mu.init <- rep(mu.init, length=length(y))
      etastart <-
        cbind(theta2eta( mu.init, .lmu  , earg = .emu  ),
              theta2eta(rho.init, .lrho , earg = .erho ))
    }
  }), list( .lmu = lmu, .lrho = lrho,
            .imu = imu, .irho = irho,
            .emu = emu, .erho = erho ))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    mu  <- eta2theta(eta[, 1], link = .lmu  , earg = .emu  )
    rho <- eta2theta(eta[, 2], link = .lrho , earg = .erho )
      pi + (rho/pi) *
      ((2*pi-mu)*sin(2*pi-mu) + cos(2*pi-mu) - mu*sin(mu) - cos(mu))
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmu , "rho" = .lrho )

    misc$earg <- list("mu" = .emu , "rho" = .erho )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    mu  <- eta2theta(eta[, 1], link = .lmu  , earg = .emu  )
    rho <- eta2theta(eta[, 2], link = .lrho , earg = .erho )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dcard(x = y, mu = mu, rho = rho, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho ))),
  vfamily = c("cardioid"),







  deriv = eval(substitute(expression({
    mu  <- eta2theta(eta[, 1], link = .lmu  , earg = .emu  )
    rho <- eta2theta(eta[, 2], link = .lrho , earg = .erho )

    dmu.deta  <- dtheta.deta(mu,  link = .lmu  , earg = .emu  )
    drho.deta <- dtheta.deta(rho, link = .lrho , earg = .erho )

    dl.dmu <-  2 * rho * sin(y-mu) / (1 + 2 * rho * cos(y-mu))
    dl.drho <- 2 * cos(y-mu) / (1 + 2 * rho * cos(y-mu))
    c(w) * cbind(dl.dmu  *  dmu.deta,
                 dl.drho * drho.deta)
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho, .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1   <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for (ii in 1:( .nsimEIM )) {
      ysim <- rcard(n, mu=mu, rho=rho)
      dl.dmu <-  2 * rho * sin(ysim-mu) / (1 + 2 * rho * cos(ysim-mu))
      dl.drho <- 2 * cos(ysim-mu) / (1 + 2 * rho * cos(ysim-mu))
      rm(ysim)
      temp3 <- cbind(dl.dmu, dl.drho)
      run.varcov <- ((ii-1) * run.varcov +
                    temp3[, ind1$row.index] *
                    temp3[, ind1$col.index]) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(run.varcov),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov

    dtheta.detas <- cbind(dmu.deta, drho.deta)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho, .nsimEIM = nsimEIM ))))
}



 vonmises <- function(llocation = extlogit(min = 0, max = 2*pi),
                      lscale  = "loge",
                      ilocation = NULL, iscale  = NULL,
                      imethod = 1, zero = NULL) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  ilocat <- ilocation


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")




  new("vglmff",
  blurb = c("Von Mises distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),
            "\n", "\n",
            "Mean:     location"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)


      predictors.names <-
        c(namesof("location", .llocat , earg = .elocat , tag = FALSE),
          namesof("scale",    .lscale , earg = .escale , tag = FALSE))

      if (!length(etastart)) {
        if ( .imethod == 1) {
          locat.init <- mean(y)
          rat10 <- sqrt((sum(w*cos(y )))^2 + sum(w*sin(y))^2) / sum(w)
          scale.init <- sqrt(1 - rat10)
        } else {
          locat.init <- median(y)
          scale.init <- sqrt(sum(w*abs(y - locat.init)) / sum(w))
        }

        locat.init <- if (length( .ilocat ))
                       rep( .ilocat , len = n) else
                       rep(locat.init, len = n)
        scale.init <- if (length( .iscale ))
                     rep( .iscale , len = n) else rep(1, len = n)
        etastart <- cbind(
            theta2eta(locat.init, .llocat , earg = .elocat ),
            theta2eta(scale.init, .lscale , earg = .escale ))
      }
      y <- y %% (2*pi)  # Coerce after initial values have been computed
  }), list( .imethod = imethod, .ilocat = ilocat,
            .escale = escale, .elocat = elocat,
            .lscale = lscale, .llocat = llocat,
            .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat ) %% (2*pi)
  }, list( .escale = escale, .lscale = lscale,
           .llocat = llocat, .elocat = elocat ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale )
    misc$earg <- list(location = .elocat , scale = .escale )



  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (Scale * cos(y - locat) -
                         log(mbesselI0(x = Scale)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .escale = escale, .lscale = lscale,
           .llocat = llocat, .elocat = elocat ))),
  vfamily = c("vonmises"),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )

    tmp6 <- mbesselI0(x = Scale, deriv = 2)
    dl.dlocat <- Scale * sin(y - locat)
    dl.dscale <- cos(y - locat) - tmp6[, 2] / tmp6[, 1]

    dlocat.deta <- dtheta.deta(locat, .llocat ,
                                 earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta)
  }), list( .escale = escale, .lscale = lscale,
            .llocat = llocat, .elocat = elocat ))),
  weight = eval(substitute(expression({
    ned2l.dlocat2 <- Scale * tmp6[, 2] / tmp6[, 1]
    ned2l.dscale2 <- tmp6[, 3] / tmp6[, 1] - (tmp6[, 2] / tmp6[, 1])^2

    wz <- matrix(0, nrow = n, ncol = 2)  # diagonal
    wz[, iam(1, 1, M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale2 * dscale.deta^2
    c(w) * wz
  }), list( .escale = escale, .elocat = elocat,
            .lscale = lscale, .llocat = llocat ))))
}



