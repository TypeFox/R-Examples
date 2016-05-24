# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












dkumar <- function(x, shape1, shape2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  N <- max(length(x), length(shape1), length(shape2))
  if (length(x)      != N) x      <- rep(x,      len = N)
  if (length(shape1) != N) shape1 <- rep(shape1, len = N)
  if (length(shape2) != N) shape2 <- rep(shape2, len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (0 <= x & x <= 1)
  logdensity[xok] <- log(shape1[xok]) + log(shape2[xok]) +
                     (shape1[xok] - 1) * log(x[xok]) +
                     (shape2[xok] - 1) * log1p(-x[xok]^shape1[xok])

  logdensity[shape1 <= 0] <- NaN
  logdensity[shape2 <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



rkumar <- function(n, shape1, shape2) {
  ans <- (1 - (runif(n))^(1/shape2))^(1/shape1)
  ans[(shape1 <= 0) | (shape2 <= 0)] <- NaN
  ans
}



qkumar <- function(p, shape1, shape2,
                   lower.tail = TRUE, log.p = FALSE) {



  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- (-expm1((1/shape2) * log(-expm1(ln.p))))^(1/shape1)
      ans[ln.p > 0] <- NaN
    } else {
      ans <- (-expm1((1/shape2) * log1p(-p)))^(1/shape1)
      ans[p < 0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- 1
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- (-expm1(ln.p / shape2))^(1/shape1)
      ans[ln.p > 0] <- NaN
      ans
    } else {
      ans <- (-expm1((1/shape2) * log(p)))^(1/shape1)
      ans[p < 0] <- NaN
      ans[p == 0] <- 1
      ans[p == 1] <- 0
      ans[p > 1] <- NaN
    }
  }
  ans[(shape1 <= 0) | (shape2 <= 0)] = NaN
  ans
}



pkumar <- function(q, shape1, shape2,
                   lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log(-expm1(shape2 * log1p(-q^shape1)))
      ans[q <= 0 ] <- -Inf
      ans[q >= 1] <- 0
    } else {
      ans <- -expm1(shape2 * log1p(-q^shape1))
      ans[q <= 0] <- 0
      ans[q >= 1] <- 1
    }
  } else {
    if (log.p) {
      ans <- shape2 * log1p(-q^shape1)
      ans[q <= 0] <- 0
      ans[q >= 1] <- -Inf
    } else {
      ans <- exp(shape2 * log1p(-q^shape1))
      ans[q <= 0] <- 1
      ans[q >= 1] <- 0
    }
  }

  ans[(shape1 <= 0) | (shape2 <= 0)] <- NaN
  ans
}










 kumar <-
  function(lshape1 = "loge", lshape2 = "loge",
           ishape1 = NULL,   ishape2 = NULL,
           grid.shape1 = c(0.4, 6.0), tol12 = 1.0e-4, zero = NULL) {
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")

  if (length(ishape1) &&
     (!is.Numeric(ishape1, length.arg = 1, positive = TRUE)))
    stop("bad input for argument 'ishape1'")
  if (length(ishape2) && !is.Numeric(ishape2))
    stop("bad input for argument 'ishape2'")

  if (!is.Numeric(tol12, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'tol12'")
  if (!is.Numeric(grid.shape1, length.arg = 2, positive = TRUE))
    stop("bad input for argument 'grid.shape1'")


  new("vglmff",
  blurb = c("Kumaraswamy distribution\n\n",
            "Links:    ", namesof("shape1", lshape1, eshape1, tag = FALSE), ", ",
                          namesof("shape2", lshape2, eshape2, tag = FALSE), "\n",
            "Mean:     shape2 * beta(1 + 1 / shape1, shape2)"),
 constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("shape1", "shape2"),
         lshape1 = .lshape1 ,
         lshape2 = .lshape2 ,
         zero = .zero )
  }, list( .zero = zero, .lshape1 = lshape1, .lshape2 = lshape2 ))),

  initialize = eval(substitute(expression({
    checklist <- w.y.check(w = w, y = y, Is.positive.y = TRUE,
                           ncol.w.max = Inf, ncol.y.max = Inf,
                           out.wy = TRUE, colsyperw = 1, maximize = TRUE)
    w <- checklist$w
    y <- checklist$y  # Now 'w' and 'y' have the same dimension.
    if (any((y <= 0) | (y >= 1)))
      stop("the response must be in (0, 1)")

    extra$ncoly <- ncoly <- ncol(y)
    extra$M1 <- M1 <- 2
    M <- M1 * ncoly
    mynames1 <- param.names("shape1", ncoly)
    mynames2 <- param.names("shape2", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lshape1 , earg = .eshape1 , tag = FALSE),
          namesof(mynames2, .lshape2 , earg = .eshape2 , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]

    if (!length(etastart)) {
      kumar.Loglikfun <- function(shape1, y, x, w, extraargs) {
        mediany <- colSums(y * w) / colSums(w)
        shape2 <- log(0.5) / log1p(-(mediany^shape1))
        sum(c(w) * dkumar(y, shape1 = shape1, shape2 = shape2, log = TRUE))
      }

      shape1.grid <- seq( .grid.shape1[1], .grid.shape1[2], len = 19)
      shape1.init <- if (length( .ishape1 )) .ishape1 else
        grid.search(shape1.grid, objfun = kumar.Loglikfun,
                    y = y,  x = x, w = w)
      shape1.init <- matrix(shape1.init, n, ncoly, byrow = TRUE)

      mediany <- colSums(y * w) / colSums(w)
      shape2.init <- if (length( .ishape2 )) .ishape2 else
        log(0.5) / log1p(-(mediany^shape1.init))
      shape2.init <- matrix(shape2.init, n, ncoly, byrow = TRUE)

      etastart <- cbind(theta2eta(shape1.init, .lshape1 , earg = .eshape1 ),
                        theta2eta(shape2.init, .lshape2 , earg = .eshape2 ))[,
                  interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .grid.shape1 = grid.shape1 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    shape2 * (base::beta(1 + 1/shape1, shape2))
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <- c(rep( .lshape1 , length = ncoly),
                   rep( .lshape2 , length = ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .eshape1
      misc$earg[[M1*ii  ]] <- .eshape2
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  loglikelihood = eval(substitute(
  function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dkumar(x = y, shape1, shape2, log = TRUE)
      if (summation) sum(ll.elts) else ll.elts
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("kumar"),
  simslot = eval(substitute(
  function(object, nsim) {
    eta <- predict(object)
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    rkumar(nsim * length(shape1), shape1 = shape1, shape2 = shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  deriv = eval(substitute(expression({
    shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2 , earg = .eshape2 )
    dshape1.deta <- dtheta.deta(shape1, link = .lshape1 , earg = .eshape1 )
    dshape2.deta <- dtheta.deta(shape2, link = .lshape2 , earg = .eshape2 )
    dl.dshape1 <- 1 / shape1 + log(y) - (shape2 - 1) * log(y) *
                  (y^shape1) / (1 - y^shape1)
    dl.dshape2 <- 1 / shape2 + log1p(-y^shape1)
    dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta,
                            dl.dshape2 * dshape2.deta)
    dl.deta[, interleave.VGAM(M, M1 = M1)]
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    ned2l.dshape11 <- (1 + (shape2 / (shape2 - 2)) *
      ((digamma(shape2) -  digamma(2))^2 -
      (trigamma(shape2) - trigamma(2)))) / shape1^2
    ned2l.dshape22 <- 1 / shape2^2
    ned2l.dshape12 <-
       (digamma(2) - digamma(1 + shape2)) / ((shape2 - 1) * shape1)

    index1 <- (abs(shape2 - 1) < .tol12 )  # Fix up singular point at shape2 == 1
    ned2l.dshape12[index1] <- -trigamma(2) / shape1[index1]
    index2 <- (abs(shape2 - 2) < .tol12 )  # Fix up singular point at shape2 == 2
    ned2l.dshape11[index2] <- (1 - 2 * psigamma(2, deriv = 2)) / shape1[index2]^2

    wz <- array(c(c(w) * ned2l.dshape11 * dshape1.deta^2,
                  c(w) * ned2l.dshape22 * dshape2.deta^2,
                  c(w) * ned2l.dshape12 * dshape1.deta * dshape2.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2, .tol12 = tol12 ))))
}






drice <- function(x, sigma, vee, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  N <- max(length(x), length(vee), length(sigma))
  if (length(x)      != N) x      <- rep(x,      len = N)
  if (length(vee)    != N) vee    <- rep(vee   , len = N)
  if (length(sigma ) != N) sigma  <- rep(sigma , len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (x > 0)
  x.abs <- abs(x[xok] * vee[xok] / sigma[xok]^2)
  logdensity[xok] <- log(x[xok]) - 2 * log(sigma[xok]) +
                     (-(x[xok]^2+vee[xok]^2)/(2*sigma[xok]^2)) +
                     log(besselI(x.abs, nu = 0, expon.scaled = TRUE)) +
                     x.abs
  logdensity[sigma <= 0] <- NaN
  logdensity[vee < 0] <- NaN

  logdensity[is.infinite(x)] <- -Inf  # 20141209 KaiH

  if (log.arg) logdensity else exp(logdensity)
}



rrice <- function(n, sigma, vee) {
  theta <- 1  # any number
  X <- rnorm(n, mean = vee * cos(theta), sd = sigma)
  Y <- rnorm(n, mean = vee * sin(theta), sd = sigma)
  sqrt(X^2 + Y^2)
}




marcumQ <- function(a, b, m = 1,
                    lower.tail = TRUE, log.p = FALSE, ... ) {
  pchisq(b^2, df = 2*m, ncp = a^2,
         lower.tail = lower.tail, log.p = log.p, ... )
}



price <- function(q, sigma, vee,
                  lower.tail = TRUE, log.p = FALSE, ...) {
  ans <- marcumQ(vee/sigma, q/sigma, m = 1,
                 lower.tail = lower.tail, log.p = log.p, ... )
  ans
}



qrice <- function(p, sigma, vee,
                  lower.tail = TRUE, log.p = FALSE, ... ) {
  sqrt(qchisq(p, df = 2, ncp = (vee/sigma)^2,
              lower.tail = lower.tail, log.p = log.p, ... )) * sigma
}









riceff.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 riceff <- function(lsigma = "loge", lvee = "loge",
                    isigma = NULL, ivee = NULL,
                    nsimEIM = 100, zero = NULL, nowarning = FALSE) {



  lvee     <- as.list(substitute(lvee))
  evee     <- link2list(lvee)
  lvee     <- attr(evee, "function.name")


  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")



  if (length(ivee) && !is.Numeric(ivee, positive = TRUE))
    stop("bad input for argument 'ivee'")
  if (length(isigma) && !is.Numeric(isigma, positive = TRUE))
    stop("bad input for argument 'isigma'")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Rice distribution\n\n",
            "Links:    ",
            namesof("sigma", lsigma, earg = esigma, tag = FALSE), ", ",
            namesof("vee",   lvee,   earg = evee,   tag = FALSE), "\n",
            "Mean:     ",
            "sigma*sqrt(pi/2)*exp(z/2)*((1-z)*",
            "besselI(-z/2, nu = 0) - z * besselI(-z/2, nu = 1)) ",
            "where z=-vee^2/(2*sigma^2)"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = FALSE,
         multipleResponses = FALSE,
         parameters.names = c("sigma", "vee"),
         nsimEIM = .nsimEIM,
         lsigma = .lsigma ,
         lvee = .lvee ,
         zero = .zero )
  }, list( .zero = zero, .lsigma = lsigma, .lvee = lvee,
           .nsimEIM = nsimEIM ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      c(namesof("sigma", .lsigma , earg = .esigma , tag = FALSE),
        namesof("vee",   .lvee   , earg = .evee   , tag = FALSE))
        



    if (!length(etastart)) {
      riceff.Loglikfun <- function(vee, y, x, w, extraargs) {
        sigma.init <- sd(rep(y, w))
        sum(c(w) * (log(y) - 2*log(sigma.init) +
                    log(besselI(y*vee/sigma.init^2, nu = 0)) -
                   (y^2 + vee^2) / (2*sigma.init^2)))
      }
    vee.grid <-
      seq(quantile(rep(y, w), probs = seq(0, 1, 0.2))["20%"],
          quantile(rep(y, w), probs = seq(0, 1, 0.2))["80%"], len = 11)
    vee.init <- if (length( .ivee )) .ivee else
      grid.search(vee.grid, objfun = riceff.Loglikfun, y = y,  x = x, w = w)
      vee.init <- rep(vee.init, length = length(y))
      sigma.init <- if (length( .isigma )) .isigma else
          sqrt(max((weighted.mean(y^2, w) - vee.init^2)/2, 0.001))
      sigma.init <- rep(sigma.init, length = length(y))

      etastart <-
        cbind(theta2eta(sigma.init, .lsigma , earg = .esigma ),
              theta2eta(vee.init,   .lvee ,   earg = .evee ))
    }
  }), list( .lvee = lvee, .lsigma = lsigma,
            .ivee = ivee, .isigma = isigma,
            .evee = evee, .esigma = esigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    vee   <- eta2theta(eta[, 1], link = .lvee ,   earg = .evee )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )
    temp9 <- -vee^2 / (2*sigma^2)


      sigma * sqrt(pi/2) *
      ((1-temp9) * besselI(-temp9/2, nu = 0, expon = TRUE) -
          temp9  * besselI(-temp9/2, nu = 1, expon = TRUE))
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),
  last = eval(substitute(expression({
    misc$link <-    c("sigma" = .lsigma , "vee" = .lvee )

    misc$earg <- list("sigma" = .esigma , "vee" = .evee )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * drice(x = y, sigma = sigma, vee = vee, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),
  vfamily = c("riceff"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )
    rrice(nsim * length(vee),
          vee = vee, sigma = sigma)
  }, list( .lvee = lvee, .lsigma = lsigma,
           .evee = evee, .esigma = esigma ))),



  deriv = eval(substitute(expression({
    sigma <- eta2theta(eta[, 1], link = .lsigma , earg = .esigma )
    vee   <- eta2theta(eta[, 2], link = .lvee   , earg = .evee )

    dvee.deta <- dtheta.deta(vee, link = .lvee , earg = .evee )
    dsigma.deta <- dtheta.deta(sigma, link = .lsigma , earg = .esigma )

    temp8 <- y * vee / sigma^2
    dl.dvee <- -vee/sigma^2 + (y/sigma^2) *
               besselI(temp8, nu = 1) / besselI(temp8, nu = 0)
    dl.dsigma <- -2/sigma + (y^2 + vee^2)/(sigma^3) -
                 (2 * temp8 / sigma) *
                 besselI(temp8, nu = 1) / besselI(temp8, nu = 0)

    c(w) * cbind(dl.dsigma * dsigma.deta,
                 dl.dvee   * dvee.deta)

  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.var <- run.cov <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rrice(n, vee = vee, sigma = sigma)
      temp8 <- ysim * vee / sigma^2
      dl.dvee <- -vee/sigma^2 + (ysim/sigma^2) *
                 besselI(temp8, nu = 1) / besselI(temp8, nu = 0)
      dl.dsigma <- -2/sigma + (ysim^2 + vee^2)/(sigma^3) -
                   (2 * temp8 / sigma) *
                   besselI(temp8, nu = 1) / besselI(temp8, nu = 0)

      rm(ysim)
      temp3 <- cbind(dl.dsigma, dl.dvee)
      run.var <- ((ii-1) * run.var + temp3^2) / ii
      run.cov <- ((ii-1) * run.cov + temp3[, 1] * temp3[, 2]) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var, run.cov)),
               n, dimm(M), byrow = TRUE) else cbind(run.var, run.cov)

    dtheta.detas <- cbind(dsigma.deta, dvee.deta)
    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list( .lvee = lvee, .lsigma = lsigma,
            .evee = evee, .esigma = esigma, .nsimEIM = nsimEIM ))))
}





dskellam <- function(x, mu1, mu2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(mu1), length(mu2))
  if (length(x)      != L) x      <- rep(x,      len = L)
  if (length(mu1)    != L) mu1    <- rep(mu1,    len = L)
  if (length(mu2)    != L) mu2    <- rep(mu2,    len = L)

  ok2 <- is.finite(mu1) & is.finite(mu2) & (mu1 >= 0) & (mu2 >= 0)
  ok3 <- (mu1 == 0) & (mu2 >  0)
  ok4 <- (mu1 >  0) & (mu2 == 0)
  ok5 <- (mu1 == 0) & (mu2 == 0)
    if (log.arg) {
      ans <- -mu1 - mu2 + 2 * sqrt(mu1*mu2) +
             0.5 * x * log(mu1) - 0.5 * x * log(mu2) +
             log(besselI(2 * sqrt(mu1*mu2),

                         nu = abs(x),

                         expon.scaled = TRUE))
      ans[ok3] <- dpois(x = -x[ok3], lambda = mu2[ok3], log = TRUE)
      ans[ok4] <- dpois(x = -x[ok4], lambda = mu1[ok4], log = TRUE)
      ans[ok5] <- dpois(x =  x[ok5], lambda = 0.0,      log = TRUE)
      ans[x != round(x)] = log(0.0)
    } else {
      ans <- (mu1/mu2)^(x/2) * exp(-mu1-mu2 + 2 * sqrt(mu1*mu2)) *
             besselI(2 * sqrt(mu1*mu2),

                     nu = abs(x),

                     expon.scaled = TRUE)
      ans[ok3] <- dpois(x = -x[ok3], lambda = mu2[ok3])
      ans[ok4] <- dpois(x = -x[ok4], lambda = mu1[ok4])
      ans[ok5] <- dpois(x =  x[ok5], lambda = 0.0)
      ans[x != round(x)] <- 0.0
    }
    ans[!ok2] <- NaN
    ans
}






rskellam <- function(n, mu1, mu2) {
  rpois(n, mu1) - rpois(n, mu2)
}



skellam.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 skellam <- function(lmu1 = "loge", lmu2 = "loge",
                     imu1 = NULL,   imu2 = NULL,
                     nsimEIM = 100, parallel = FALSE, zero = NULL) {

  lmu1 <- as.list(substitute(lmu1))
  emu1 <- link2list(lmu1)
  lmu1 <- attr(emu1, "function.name")

  lmu2 <- as.list(substitute(lmu2))
  emu2 <- link2list(lmu2)
  lmu2 <- attr(emu2, "function.name")


  if (length(imu1) &&
      !is.Numeric(imu1, positive = TRUE))
    stop("bad input for argument 'imu1'")
  if (length(imu2) &&
      !is.Numeric(imu2, positive = TRUE))
    stop("bad input for argument 'imu2'")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")

  new("vglmff",
  blurb = c("Skellam distribution\n\n",
         "Links:    ",
         namesof("mu1", lmu1, earg = emu1, tag = FALSE), ", ", 
         namesof("mu2", lmu2, earg = emu2, tag = FALSE), "\n",
         "Mean:     mu1-mu2", "\n",
         "Variance: mu1+mu2"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel , 
                           constraints = constraints,
                           apply.int = TRUE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = FALSE,
         multipleResponses = FALSE,
         parameters.names = c("mu1", "mu2"),
         nsimEIM = .nsimEIM,
         lmu1 = .lmu1 ,
         lmu2 = .lmu2 ,
         zero = .zero )
  }, list( .zero = zero, .lmu1 = lmu1, .lmu2 = lmu2,
           .nsimEIM = nsimEIM ))),
  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
      namesof("mu1", .lmu1, earg = .emu1, tag = FALSE),
      namesof("mu2", .lmu2, earg = .emu2, tag = FALSE))


    if (!length(etastart)) {
      junk <- lm.wfit(x = x, y = c(y), w = c(w))
      var.y.est <- sum(c(w) * junk$resid^2) / junk$df.residual
      mean.init <- weighted.mean(y, w)

      mu1.init <- max((var.y.est + mean.init) / 2, 0.01)
      mu2.init <- max((var.y.est - mean.init) / 2, 0.01)
      mu1.init <- rep(if (length( .imu1 )) .imu1 else mu1.init,
                      length = n)
      mu2.init <- rep(if (length( .imu2 )) .imu2 else mu2.init,
                      length = n)

      etastart <- cbind(theta2eta(mu1.init, .lmu1, earg = .emu1 ),
                        theta2eta(mu2.init, .lmu2, earg = .emu2 ))
      }
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .imu1 = imu1, .imu2 = imu2,
            .emu1 = emu1, .emu2 = emu2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
      mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
      mu1 - mu2
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2 ))),
  last = eval(substitute(expression({
      misc$link <-    c("mu1" = .lmu1, "mu2" = .lmu2)

      misc$earg <- list("mu1" = .emu1, "mu2" = .emu2 )

      misc$expected <- TRUE
      misc$nsimEIM <- .nsimEIM
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .emu1 = emu1, .emu2 = emu2,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      ll.elts <-
        if ( is.logical( .parallel ) &&
             length( .parallel ) == 1 &&
             .parallel )
          c(w) * log(besselI(2*mu1, nu = y, expon = TRUE)) else
          c(w) * (-mu1 - mu2 +
                  0.5 * y * log(mu1) -
                  0.5 * y * log(mu2) +
                  2 * sqrt(mu1*mu2) +  # Use this when expon = TRUE
                  log(besselI(2 * sqrt(mu1*mu2), nu = y, expon = TRUE)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2,
           .parallel = parallel ))),
  vfamily = c("skellam"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )
    rskellam(nsim * length(mu1), mu1, mu2)
  }, list( .lmu1 = lmu1, .lmu2 = lmu2,
           .emu1 = emu1, .emu2 = emu2,
           .parallel = parallel ))),





  deriv = eval(substitute(expression({
    mu1 <- eta2theta(eta[, 1], link = .lmu1, earg = .emu1 )
    mu2 <- eta2theta(eta[, 2], link = .lmu2, earg = .emu2 )

    dmu1.deta <- dtheta.deta(mu1, link = .lmu1, earg = .emu1 )
    dmu2.deta <- dtheta.deta(mu2, link = .lmu2, earg = .emu2 )

    temp8 <- 2 * sqrt(mu1*mu2)
    temp9 <-  besselI(temp8, nu = y  , expon = TRUE)
    temp7 <- (besselI(temp8, nu = y-1, expon = TRUE) +
              besselI(temp8, nu = y+1, expon = TRUE)) / 2
    temp6 <- temp7 / temp9

    dl.dmu1 <- -1 + 0.5 * y / mu1 + sqrt(mu2/mu1) * temp6
    dl.dmu2 <- -1 - 0.5 * y / mu2 + sqrt(mu1/mu2) * temp6

    c(w) * cbind(dl.dmu1 * dmu1.deta,
                 dl.dmu2 * dmu2.deta)
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .emu1 = emu1, .emu2 = emu2,
            .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.var <- run.cov <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rskellam(n, mu1=mu1, mu2=mu2)
      temp9 <-  besselI(temp8, nu = ysim,   expon = TRUE)
      temp7 <- (besselI(temp8, nu = ysim-1, expon = TRUE) +
                besselI(temp8, nu = ysim+1, expon = TRUE)) / 2
      temp6 <- temp7 / temp9
      dl.dmu1 <- -1 + 0.5 * ysim/mu1 + sqrt(mu2/mu1) * temp6
      dl.dmu2 <- -1 - 0.5 * ysim/mu2 + sqrt(mu1/mu2) * temp6
      rm(ysim)
      temp3 <- cbind(dl.dmu1, dl.dmu2)
      run.var <- ((ii-1) * run.var + temp3^2) / ii
      run.cov <- ((ii-1) * run.cov + temp3[, 1] * temp3[, 2]) / ii
    }
    wz <- if (intercept.only)
          matrix(colMeans(cbind(run.var, run.cov)),
                 n, dimm(M), byrow = TRUE) else
          cbind(run.var, run.cov)

    dtheta.detas <- cbind(dmu1.deta, dmu2.deta)
    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list( .lmu1 = lmu1, .lmu2 = lmu2,
            .emu1 = emu1, .emu2 = emu2,
            .nsimEIM = nsimEIM ))))
}




dyules <- function(x, rho, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if ( log.arg ) {
    ans <- log(rho) + lbeta(abs(x), rho+1)
    ans[(x != round(x)) | (x < 1)] <- log(0)
  } else {
    ans <- rho * beta(x, rho+1)
    ans[(x != round(x)) | (x < 1)] <- 0
  }
  ans[!is.finite(rho) | (rho <= 0) | (rho <= 0)] <- NA
  ans
}




ryules <- function(n, rho) {
  rgeom(n, prob = exp(-rexp(n, rate = rho))) + 1
}




pyules <- function(q, rho, log.p = FALSE) {


  tq <- trunc(q)
  ans <- 1 - tq * beta(abs(tq), rho+1)
  ans[q < 1] <- 0
  ans[is.infinite(q) & q > 0] <- 1  # 20141215 KaiH
  ans[(rho <= 0) | (rho <= 0)] <- NA
  if (log.p) log(ans) else ans
  ans
}




yulesimon.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 yulesimon <- function(link = "loge",
                       irho = NULL, nsimEIM = 200,
                       zero = NULL) {

  if (length(irho) &&
      !is.Numeric(irho, positive = TRUE))
    stop("argument 'irho' must be > 0")



  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")



  new("vglmff",
  blurb = c("Yule-Simon distribution f(y) = rho * beta(y, rho + 1), ",
            "rho > 0, y = 1, 2,..\n\n",
            "Link:    ",
            namesof("rho", link, earg = earg), "\n\n",
            "Mean:     rho / (rho - 1), provided rho > 1\n",
            "Variance: rho^2 / ((rho - 1)^2 * (rho - 2)), ",
            "provided rho > 2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         nsimEIM = .nsimEIM,
         parameters.names = c("rho"),
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)

    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("rho", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
      wmeany <- colSums(y * w) / colSums(w) + 1/8

      rho.init <- wmeany / (wmeany - 1)
      rho.init <- matrix(if (length( .irho )) .irho else
                         rho.init, n, M, byrow = TRUE)
      etastart <- theta2eta(rho.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .irho = irho ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- rho <- eta2theta(eta, .link , earg = .earg )
    ans[rho >  1] <- rho / (rho - 1)
    ans[rho <= 1] <- NA
    ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$irho <- .irho
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM,
            .irho = irho ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    rho <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dyules(x = y, rho = rho, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("yulesimon"),





  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    rho <- eta2theta(eta, .link , earg = .earg )
    ryules(nsim * length(rho), rho = rho)
  }, list( .link = link, .earg = earg ))),







  deriv = eval(substitute(expression({
    M1 <- 1
    rho <- eta2theta(eta, .link , earg = .earg )
    dl.drho <- 1/rho + digamma(1+rho) - digamma(1+rho+y)
    drho.deta <- dtheta.deta(rho, .link , earg = .earg )
    c(w) * dl.drho * drho.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({

    run.var <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- ryules(n, rho <- rho)
      dl.drho <- 1/rho + digamma(1+rho) - digamma(1+rho+ysim)
      rm(ysim)
      temp3 <- dl.drho
      run.var <- ((ii-1) * run.var + temp3^2) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, M, byrow = TRUE) else cbind(run.var)

    wz <- wz * drho.deta^2


    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}







dlind <- function(x, theta, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if ( log.arg ) {
    ans <- 2 * log(theta) + log1p(x) - theta * x - log1p(theta)
    ans[x < 0 | is.infinite(x)] <- log(0)  # 20141209 KaiH
  } else {
    ans <- theta^2 * (1 + x) * exp(-theta * x) / (1 + theta)
    ans[x < 0 | is.infinite(x)] <- 0  # 20141209 KaiH
  }
  ans[theta <= 0] <- NaN
  ans
}



plind <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ans <- log(-expm1(-theta * q + log1p(q / (1 + 1/theta))))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- -expm1(-theta * q + log1p(q / (1 + 1/theta)))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- -theta * q + log1p(q / (1 + 1/theta))
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- exp(-theta * q + log1p(q / (1 + 1/theta)))
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  }
  ans[theta <= 0] <- NaN
  ans
}






rlind <- function(n, theta) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE)) 
             stop("bad input for argument 'n'") else n



  ifelse(runif(use.n) < rep(1 / (1 + 1/theta), length = use.n),
         rexp(use.n, theta),
         rgamma(use.n, shape = 2, scale = 1 / theta))
}



 lindley <- function(link = "loge",
                     itheta = NULL, zero = NULL) {


  if (length(itheta) &&
      !is.Numeric(itheta, positive = TRUE))
    stop("argument 'itheta' must be > 0")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")




  new("vglmff",
  blurb = c("Lindley distribution f(y) = ",
            "theta^2 * (1 + y) * exp(-theta * y) / (1 + theta), ",  
            "theta > 0, y > 0,\n\n",
            "Link:    ",
            namesof("theta", link, earg = earg), "\n\n",
            "Mean:     (theta + 2) / (theta * (theta + 1))\n",
            "Variance: (theta^2 + 4 * theta + 2) / (theta * (theta + 1))^2"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("theta"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)

    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("theta", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
      wmeany <- colSums(y * w) / colSums(w) + 1/8


      theta.init <- 1 / (wmeany + 1)
      theta.init <- matrix(if (length( .itheta )) .itheta else
                           theta.init, n, M, byrow = TRUE)
      etastart <- theta2eta(theta.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .itheta = itheta ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta <- eta2theta(eta, .link , earg = .earg )
    (theta + 2) / (theta * (theta + 1))
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$itheta <- .itheta
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg,
            .itheta = itheta ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    theta <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlind(x = y, theta = theta, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("lindley"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    theta <- eta2theta(eta, .link , earg = .earg )
    rlind(nsim * length(theta), theta = theta)
  }, list( .link = link, .earg = earg ))),



  deriv = eval(substitute(expression({
    M1 <- 1
    theta <- eta2theta(eta, .link , earg = .earg )

    dl.dtheta <- 2 / theta - 1 / (1 + theta) - y

    DTHETA.DETA <- dtheta.deta(theta, .link , earg = .earg )

    c(w) * dl.dtheta * DTHETA.DETA
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({

    ned2l.dtheta2 <- (theta^2 + 4 * theta + 2) / (theta * (1 + theta))^2

    c(w) * ned2l.dtheta2 * DTHETA.DETA^2
  }), list( .zero = zero ))))
}






dpoislindley <- function(x, theta, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if ( log.arg ) {
    ans <- 2 * log(theta) + log(theta + 2 + x) -
           (x+3) * log1p(theta)
    ans[(x != round(x)) | (x < 0)] <- log(0)
  } else {
    ans <- theta^2 * (theta + 2 + x) / (theta + 1)^(x+3)
    ans[(x != round(x)) | (x < 0)] <- 0
  }
  ans[ # !is.finite(theta) |
     (theta <= 0)] <- NA
  ans
}


if (FALSE)
rpoislindley <- function(n, theta) {
}


if (FALSE)
ppoislindley <- function(q, theta) {
}



if (FALSE)
poislindley.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


if (FALSE)
 poissonlindley <-
  function(link = "loge",
           itheta = NULL, nsimEIM = 200,
           zero = NULL) {

  stop("not working since rpoislindley() not written")



  if (length(itheta) &&
      !is.Numeric(itheta, positive = TRUE))
    stop("argument 'itheta' must be > 0")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")



  new("vglmff",
  blurb = c("Poisson-Lindley distribution f(y) = ",
            "theta^2 * (theta + 2 + y) / (theta + 1)^(y+3), ",  
            "theta > 0, y = 0, 1, 2,..\n\n",
            "Link:    ",
            namesof("theta", link, earg = earg), "\n\n",
            "Mean:     (theta + 2) / (theta * (theta + 1)),\n",
            "Variance: (theta^3 + 4 * theta^2 + 6 * theta + 2) / ",
                      "(theta * (theta + 1))^2, "
            ),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         nsimEIM = .nsimEIM ,
         parameters.names = c("theta"),
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)

    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1  <- param.names("theta", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg , tag = FALSE)

    if (!length(etastart)) {
      wmeany <- colSums(y * w) / colSums(w) + 1/8

      MOM <- (sqrt((wmeany - 1)^2 + 8 * wmeany) -
              wmeany + 1) / (2 * wmeany)
      MOM[MOM < 0.01] <- 0.01


      theta.init <- MOM
      theta.init <- matrix(if (length( .itheta )) .itheta else
                           theta.init, n, M, byrow = TRUE)
      etastart <- theta2eta(theta.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .itheta = itheta ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta <- eta2theta(eta, .link , earg = .earg )
    (theta + 2) / (theta * (theta + 1))
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$itheta <- .itheta
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM,
            .itheta = itheta ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    theta = eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dpoislindley(x = y, theta = theta, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("poissonlindley"),
  deriv = eval(substitute(expression({
    M1 <- 1
    theta <- eta2theta(eta, .link , earg = .earg )

    dl.dtheta <- 2 / theta + 1 / (y + 2 + theta) - (y + 3) / (theta + 1)

    DTHETA.DETA <- dtheta.deta(theta, .link , earg = .earg )

    c(w) * dl.dtheta * DTHETA.DETA
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({

    run.var <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- rpoislindley(n, theta = theta)
      dl.dtheta <- 2 / theta + 1 / (ysim + 2 + theta) -
                   (ysim + 3) / (theta + 1)
      rm(ysim)
      temp3 <- dl.dtheta
      run.var <- ((ii-1) * run.var + temp3^2) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(cbind(run.var)),
               n, M, byrow = TRUE) else cbind(run.var)

    wz <- wz * DTHETA.DETA^2


    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}








dslash <- function(x, mu = 0, sigma = 1, log = FALSE,
                   smallno = .Machine$double.eps * 1000) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (!is.Numeric(sigma) || any(sigma <= 0))
    stop("argument 'sigma' must be positive")
  L <- max(length(x), length(mu), length(sigma))
  if (length(x)     != L) x     <- rep(x,     len = L)
  if (length(mu)    != L) mu    <- rep(mu,    len = L)
  if (length(sigma) != L) sigma <- rep(sigma, len = L)

  zedd <- (x-mu)/sigma
  if (log.arg) {
    ifelse(abs(zedd) < smallno,
           -log(2*sigma*sqrt(2*pi)),
           log1p(-exp(-zedd^2/2)) - log(sqrt(2*pi)*sigma*zedd^2))
  } else {
    ifelse(abs(zedd) < smallno,
           1/(2*sigma*sqrt(2*pi)),
           -expm1(-zedd^2/2)/(sqrt(2*pi)*sigma*zedd^2))
  }
}




pslash <- function(q, mu = 0, sigma = 1, very.negative = -10000,
                   lower.tail = TRUE, log.p = FALSE) {
  if (any(is.na(q)))
    stop("argument 'q' must have non-missing values")
  if (!is.Numeric(mu))
    stop("argument 'mu' must have finite and non-missing values")
  if (!is.Numeric(sigma, positive = TRUE))
    stop("argument 'sigma' must have positive finite non-missing values")
  if (!is.Numeric(very.negative, length.arg = 1) ||
     (very.negative >= 0))
    stop("argument 'very.negative' must be quite negative")

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  L <- max(length(q), length(mu), length(sigma))
  if (length(q)     != L) q     <- rep(q,     len = L)
  if (length(mu)    != L) mu    <- rep(mu,    len = L)
  if (length(sigma) != L) sigma <- rep(sigma, len = L)

  zedd <- (q - mu)/sigma
  ans <- as.numeric(q * NA)
  extreme.q <- FALSE
  for (ii in 1:L) {
    use.trick <- (-abs(zedd[ii]) <= very.negative)
    if (use.trick) {
      ans[ii] <- ifelse(zedd[ii] < 0, 0.0, 1.0)
      extreme.q <- TRUE
    } else
    if ((zedd[ii] >= very.negative) &&
         zedd[ii] <= 0.0) {
      temp2 <- integrate(dslash, lower = q[ii], upper = mu[ii],
                         mu = mu[ii], sigma = sigma[ii])
      if (temp2$message != "OK")
        warning("integrate() failed on 'temp2'")
      ans[ii] <- 0.5 - temp2$value
    } else {
      temp1 <- integrate(dslash, lower = mu[ii], upper =  q[ii],
                         mu = mu[ii], sigma = sigma[ii])
      if (temp1$message != "OK")
        warning("integrate() failed")
      ans[ii] <- 0.5 + temp1$value
    }
  }
  if (extreme.q)
    warning("returning 0 or 1 values for extreme values of argument 'q'")

  if (lower.tail) {
    if (log.p) log(ans) else ans
  } else {
    if (log.p) log1p(-ans) else -expm1(log(ans))
  }
}




rslash <- function (n, mu = 0, sigma = 1) {
  rnorm(n = n, mean = mu, sd = sigma) / runif(n = n)
}



slash.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 slash <- function(lmu = "identitylink", lsigma = "loge",
                   imu = NULL, isigma = NULL,
                   iprobs = c(0.1, 0.9),
                   nsimEIM = 250, zero = NULL,
                   smallno = .Machine$double.eps * 1000) {

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")


  if (length(isigma) &&
      !is.Numeric(isigma, positive = TRUE))
    stop("argument 'isigma' must be > 0")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
    stop("argument 'nsimEIM' should be an integer greater than 50")

  if (!is.Numeric(iprobs, positive = TRUE) ||
      max(iprobs) >= 1 ||
      length(iprobs) != 2)
    stop("bad input for argument 'iprobs'")
  if (!is.Numeric(smallno, positive = TRUE) ||
      smallno > 0.1)
    stop("bad input for argument 'smallno'")


  new("vglmff",
  blurb = c("Slash distribution\n\n",
         "Links:    ",
         namesof("mu",    lmu,    earg = emu,    tag = FALSE), ", ",
         namesof("sigma", lsigma, earg = esigma, tag = FALSE), "\n",
         paste(
         "1-exp(-(((y-mu)/sigma)^2)/2))/(sqrt(2*pi)*",
         "sigma*((y-mu)/sigma)^2)",
         "\ty!=mu",
         "\n1/(2*sigma*sqrt(2*pi))",
         "\t\t\t\t\t\t\ty=mu\n")),

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
         parameters.names = c("mu", "sigma"),
         lmu    = .lmu ,
         lsigma = .lsigma ,
         zero = .zero )
  }, list( .zero = zero, .lmu = lmu, .lsigma = lsigma ))),


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
        namesof("mu",    .lmu ,    earg = .emu,    tag = FALSE),
        namesof("sigma", .lsigma , earg = .esigma, tag = FALSE))


    if (!length(etastart)) {

      slash.Loglikfun <- function(mu, y, x, w, extraargs) {
          sigma <- if (is.Numeric(.isigma)) .isigma else
            max(0.01,
               ((quantile(rep(y, w), prob = 0.75)/2)-mu)/qnorm(0.75))
          zedd <- (y-mu)/sigma
          sum(c(w) * ifelse(abs(zedd)<.smallno,
                         -log(2*sigma*sqrt(2*pi)),
                         log1p(-exp(-zedd^2/2)) -
                         log(sqrt(2*pi) * sigma * zedd^2)))
      }
      iprobs <- .iprobs
      mu.grid <- quantile(rep(y, w), probs=iprobs)
      mu.grid <- seq(mu.grid[1], mu.grid[2], length=100)
      mu.init <- if (length( .imu )) .imu else
                 grid.search(mu.grid, objfun = slash.Loglikfun,
                             y = y,  x = x, w = w)
      sigma.init <- if (is.Numeric(.isigma)) .isigma else
        max(0.01,
           ((quantile(rep(y, w), prob = 0.75)/2) -
                      mu.init) / qnorm(0.75))
      mu.init <- rep(mu.init, length = length(y))
      etastart <- matrix(0, n, 2)
      etastart[, 1] <- theta2eta(mu.init, .lmu , earg = .emu )
      etastart[, 2] <- theta2eta(sigma.init, .lsigma , earg = .esigma )
    }
  }), list( .lmu = lmu, .lsigma = lsigma,
            .imu = imu, .isigma = isigma,
            .emu = emu, .esigma = esigma,
            .iprobs = iprobs, .smallno = smallno))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      NA * eta2theta(eta[, 1], link = .lmu , earg = .emu )
  }, list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmu , "sigma" = .lsigma )

    misc$earg <- list("mu" = .emu , "sigma" = .esigma )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )
    zedd <- (y - mu) / sigma
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dslash(x = y, mu = mu, sigma = sigma, log = TRUE,
                               smallno = .smallno)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lsigma = lsigma,
           .emu = emu, .esigma = esigma, .smallno = smallno ))),
  vfamily = c("slash"),





  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )
    rslash(nsim * length(sigma), mu = mu, sigma = sigma)
  }, list( .lmu = lmu, .lsigma = lsigma,
           .emu = emu, .esigma = esigma, .smallno = smallno ))),




  deriv = eval(substitute(expression({
    mu    <- eta2theta(eta[, 1], link = .lmu    , earg = .emu    )
    sigma <- eta2theta(eta[, 2], link = .lsigma , earg = .esigma )

    dmu.deta    <- dtheta.deta(mu,    link = .lmu    , earg = .emu    )
    dsigma.deta <- dtheta.deta(sigma, link = .lsigma , earg = .esigma )

    zedd <- (y - mu) / sigma
    d3 <- deriv3(~ w * log(1 - exp(-(((y - mu) / sigma)^2) / 2)) -
                 log(sqrt(2 * pi) * sigma * ((y - mu) / sigma)^2),
                 c("mu", "sigma"))
    eval.d3 <- eval(d3)
    dl.dthetas <-  attr(eval.d3, "gradient")
    dl.dmu    <- dl.dthetas[, 1]
    dl.dsigma <- dl.dthetas[, 2]
    ind0 <- (abs(zedd) < .smallno)
    dl.dmu[ind0] <- 0
    dl.dsigma[ind0] <- -1 / sigma[ind0]
    c(w) * cbind(dl.dmu * dmu.deta, dl.dsigma * dsigma.deta)
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma, .smallno = smallno ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    sd3 <- deriv3(~ w * log(1 - exp(-(((ysim - mu) / sigma)^2) / 2))-
                  log(sqrt(2 * pi) * sigma * ((ysim - mu) / sigma)^2),
                  c("mu", "sigma"))
    for (ii in 1:( .nsimEIM )) {
      ysim <- rslash(n, mu = mu, sigma = sigma)
      seval.d3 <- eval(sd3)

      dl.dthetas <-  attr(seval.d3, "gradient")
      dl.dmu    <- dl.dthetas[, 1]
      dl.dsigma <- dl.dthetas[, 2]

      temp3 <- cbind(dl.dmu, dl.dsigma)
      run.varcov <- run.varcov + temp3[, ind1$row] * temp3[, ind1$col]
    }
    run.varcov <- run.varcov / .nsimEIM
    wz <- if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = FALSE),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov
    dthetas.detas <- cbind(dmu.deta, dsigma.deta)
    wz <- wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
    c(w) * wz
  }), list( .lmu = lmu, .lsigma = lsigma,
            .emu = emu, .esigma = esigma,
            .nsimEIM = nsimEIM, .smallno = smallno ))))
}




dnefghs <- function(x, tau, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(tau))
  if (length(x)   != N) x   <- rep(x,   len = N)
  if (length(tau) != N) tau <- rep(tau, len = N)

  logdensity <- log(sin(pi*tau)) + (1-tau)*x - log(pi) - log1pexp(x)
  logdensity[tau < 0] <- NaN
  logdensity[tau > 1] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



 nefghs <- function(link = "logit",
                    itau = NULL, imethod = 1) {

  if (length(itau) &&
      !is.Numeric(itau, positive = TRUE) ||
      any(itau >= 1))
    stop("argument 'itau' must be in (0, 1)")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
       imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Natural exponential family generalized hyperbolic ",
            "secant distribution\n",
            "f(y) = sin(pi*tau)*exp((1-tau)*y)/(pi*(1+exp(y))\n\n",
            "Link:    ",
            namesof("tau", link, earg = earg), "\n\n",
            "Mean:     pi / tan(pi * tau)\n"),
  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
      namesof("tau", .link , earg = .earg , tag = FALSE) 


    if (!length(etastart)) {
      wmeany <- if ( .imethod == 1) weighted.mean(y, w) else
                median(rep(y, w))
      if (abs(wmeany) < 0.01)
        wmeany <- 0.01
      tau.init <- atan(pi / wmeany) / pi + 0.5
      tau.init[tau.init < 0.03] <- 0.03
      tau.init[tau.init > 0.97] <- 0.97
      tau.init <- rep(if (length( .itau )) .itau else tau.init,
                      len = n)
      etastart <- theta2eta(tau.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg,
            .itau = itau,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    tau <- eta2theta(eta, .link , earg = .earg )
    pi / tan(pi * tau)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(tau = .link)

    misc$earg <- list(tau = .earg )

    misc$expected <- TRUE
    misc$imethod <- .imethod
  }), list( .link = link, .earg = earg,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    tau <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnefghs(x = y, tau = tau, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("nefghs"),
  deriv = eval(substitute(expression({
    tau <- eta2theta(eta, .link , earg = .earg )
    dl.dtau <- pi / tan(pi * tau) - y
    dtau.deta <- dtheta.deta(tau, .link , earg = .earg )
    w * dl.dtau * dtau.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.dtau2 <- (pi / sin(pi * tau))^2
    wz <- ned2l.dtau2 * dtau.deta^2
    c(w) * wz
  }), list( .link = link ))))
}




dlogF <- function(x, shape1, shape2, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  logdensity <- shape1*x - lbeta(shape1, shape2) -
                (shape1 + shape2) * log1pexp(x)

  logdensity[is.infinite(x)] <- -Inf  # 20141209 KaiH

  if (log.arg) logdensity else exp(logdensity)
}




 logF <- function(lshape1 = "loge", lshape2 = "loge",
                  ishape1 = NULL, ishape2 = 1,
                  imethod = 1) {

  if (length(ishape1) &&
      !is.Numeric(ishape1, positive = TRUE))
    stop("argument 'ishape1' must be positive")
  if ( # length(ishape2) &&
      !is.Numeric(ishape2, positive = TRUE))
    stop("argument 'ishape2' must be positive")


  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")


  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
      stop("argument 'imethod' must be 1 or 2")

  new("vglmff",
  blurb = c("log F distribution\n",
            "f(y) = exp(-shape2 * y) / (beta(shape1, shape2) * ",
            "(1 + exp(-y))^(shape1 + shape2))\n\n",
            "Link:    ",
            namesof("shape1", lshape1, earg = eshape1), ", ",
            namesof("shape2", lshape2, earg = eshape2), "\n\n",
            "Mean:     digamma(shape1) - digamma(shape2)"),
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
      namesof("shape1", .lshape1 , earg = .eshape1 , tag = FALSE),
      namesof("shape2", .lshape2 , earg = .eshape2 , tag = FALSE))


    if (!length(etastart)) {
      wmeany <- if ( .imethod == 1) weighted.mean(y, w) else
                median(rep(y, w))


      shape1.init <- shape2.init <- rep( .ishape2 , len = n)
      shape1.init <- if (length( .ishape1))
                            rep( .ishape1, len = n) else {
                index1 <- (y > wmeany)
                shape1.init[ index1] <- shape2.init[ index1] + 1/1
                shape1.init[!index1] <- shape2.init[!index1] - 1/1
                shape1.init <- pmax(shape1.init, 1/8)
                shape1.init
              }
      etastart <-
          cbind(theta2eta(shape1.init, .lshape1 , earg = .eshape1 ),
                theta2eta(shape2.init, .lshape2 , earg = .eshape2 ))
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )
    digamma(shape1) - digamma(shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  last = eval(substitute(expression({
    misc$link <-    c(shape1 = .lshape1 , shape2 = .lshape2)

    misc$earg <- list(shape1 = .eshape1 , shape2 = .eshape2 )

    misc$expected <- TRUE
    misc$imethod <- .imethod
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlogF(x = y, shape1 = shape1,
                              shape2 = shape2, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .eshape1 = eshape1, .eshape2 = eshape2 ))),
  vfamily = c("logF"),












  deriv = eval(substitute(expression({
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .eshape1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .eshape2 )

    tmp888 <- digamma(shape1 + shape2) - log1pexp(-y)
    dl.dshape1 <- tmp888 - digamma(shape1)
    dl.dshape2 <- tmp888 - digamma(shape2) - y

    dshape1.deta <- dtheta.deta(shape1, .lshape1 , earg = .eshape1 )
    dshape2.deta <- dtheta.deta(shape2, .lshape2 , earg = .eshape2 )

    c(w) * cbind(dl.dshape1 * dshape1.deta,
                 dl.dshape2 * dshape2.deta)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))),
  weight = eval(substitute(expression({
    tmp888 <- trigamma(shape1 + shape2)
    ned2l.dshape12 <- trigamma(shape1) - tmp888
    ned2l.dshape22 <- trigamma(shape2) - tmp888
    ned2l.dshape1shape2 <- -tmp888

    wz <- matrix(0, n, dimm(M))
    wz[,iam(1, 1, M = M)] <- ned2l.dshape12 * dshape1.deta^2
    wz[,iam(2, 2, M = M)] <- ned2l.dshape22 * dshape2.deta^2
    wz[,iam(1, 2, M = M)] <- ned2l.dshape1shape2 * dshape1.deta *
                                                   dshape2.deta
    c(w) * wz
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .eshape1 = eshape1, .eshape2 = eshape2 ))))
}







dbenf <- function(x, ndigits = 1, log = FALSE) {
  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  ans <- x * NA
  indexTF <- is.finite(x) & (x >= lowerlimit)

  ans[indexTF] <- log10(1 + 1/x[indexTF])
  ans[!is.na(x) & !is.nan(x) &
     ((x < lowerlimit) |
      (x > upperlimit) |
      (x != round(x)))] <- 0.0
  if (log.arg) log(ans) else ans
}



rbenf <- function(n, ndigits = 1) {
  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE)) 
             stop("bad input for argument 'n'") else n
  myrunif <- runif(use.n)

  ans <- rep(lowerlimit, length = use.n)
  for (ii in (lowerlimit+1):upperlimit) {
      indexTF <- (pbenf(ii-1, ndigits = ndigits) < myrunif) &
                 (myrunif <= pbenf(ii, ndigits = ndigits))
      ans[indexTF] <- ii
  }
  ans
}



pbenf <- function(q, ndigits = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)

  ans <- q * NA
  floorq <- floor(q)
  indexTF <- is.finite(q) & (floorq >= lowerlimit)

  if (ndigits == 1) {
    if (lower.tail) {
      if (log.p) {
        ans[indexTF] <- log(log10(1 + floorq[indexTF]))
        ans[q <  lowerlimit ] <- -Inf
        ans[q >= upperlimit] <- 0
      } else {
        ans[indexTF] <- log10(1 + floorq[indexTF])
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- 1
      }
    } else {
      if (log.p) {
        ans[indexTF] <- log1p(-log10(1 + floorq[indexTF]))
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- -Inf
      } else {
        ans[indexTF] <- log10(10 / (1 + floorq[indexTF]))
        ans[q <  lowerlimit] <- 1
        ans[q >= upperlimit] <- 0
      }
    }
  } else {
    if (lower.tail) {
      if (log.p) {
        ans[indexTF] <- log(log10((1 + floorq[indexTF])/10))
        ans[q <  lowerlimit ] <- -Inf
        ans[q >= upperlimit] <- 0
      } else {
        ans[indexTF] <- log10((1 + floorq[indexTF])/10)
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- 1
     }
    } else {
      if (log.p) {
        ans[indexTF] <- log(log10(100/(1 + floorq[indexTF])))
        ans[q <  lowerlimit] <- 0
        ans[q >= upperlimit] <- -Inf
      } else {
        ans[indexTF] <- log10(100/(1 + floorq[indexTF]))
        ans[q <  lowerlimit] <- 1
        ans[q >= upperlimit] <- 0
      }
    }
  }
  ans
}



if (FALSE)
qbenf <- function(p, ndigits = 1) {

  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")
  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  bad <- !is.na(p) & !is.nan(p) & ((p < 0) | (p > 1))
  if (any(bad))
    stop("bad input for argument 'p'")

  ans <- rep(lowerlimit, length = length(p))
  for (ii in (lowerlimit+1):upperlimit) {
    indexTF <- is.finite(p) &
              (pbenf(ii-1, ndigits = ndigits) < p) &
              (p <= pbenf(ii, ndigits = ndigits))
    ans[indexTF] <- ii
  }

  ans[ is.na(p) |  is.nan(p)] <- NA
  ans[!is.na(p) & !is.nan(p) & (p == 0)] <- lowerlimit
  ans[!is.na(p) & !is.nan(p) & (p == 1)] <- upperlimit
  ans
}




qbenf <- function(p, ndigits = 1,
                  lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(ndigits, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      ndigits > 2)
    stop("argument 'ndigits' must be 1 or 2")

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (log.p) { 
    bad <- ((p > 0) | is.na(p) | is.nan(p))
  } else {
    bad <- ((p < 0) | (p > 1) | is.na(p) | is.nan(p))
  }
  if (any(bad))
    stop("bad input for argument 'p'")

  lowerlimit <- ifelse(ndigits == 1, 1, 10)
  upperlimit <- ifelse(ndigits == 1, 9, 99)
  ans <- rep(lowerlimit, length = length(p))

  if (lower.tail) {
    for (ii in (lowerlimit+1):upperlimit) {
      indexTF <- is.finite(p) &
                 (pbenf(ii-1, ndigits = ndigits,
                        lower.tail = lower.tail, log.p = log.p) < p) &
              (p <= pbenf(ii, ndigits = ndigits,
                          lower.tail = lower.tail, log.p = log.p))
      ans[indexTF] <- ii
    }
  } else {  ## when lower.tail = F, pbenf(ii-1) >= p & pben(ii) < p
    for (ii in (lowerlimit+1):upperlimit) {
      indexTF <- is.finite(p) &
                 (pbenf(ii-1, ndigits = ndigits,
                        lower.tail = lower.tail, log.p = log.p) >= p) &
                 (p > pbenf(ii, ndigits = ndigits,
                            lower.tail = lower.tail, log.p = log.p))
      ans[indexTF] <- ii
    }
  }

  if (lower.tail) {
    if (log.p) {
      ans[p > 0] <- NaN
      ans[p == -Inf] <- lowerlimit
    } else {
      ans[p < 0] <- NaN
      ans[p == 0] <- lowerlimit
      ans[p == 1] <- upperlimit
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ans[p > 0] <- NaN
      ans[p == -Inf] <- upperlimit
    } else {
      ans[p < 0] <- NaN
      ans[p == 0] <- upperlimit
      ans[p == 1] <- lowerlimit
      ans[p > 1] <- NaN
    }
  }
  ans
}


















 truncgeometric <-
  function(upper.limit = Inf,  # lower.limit = 1,  # Inclusive
           link = "logit", expected = TRUE,
           imethod = 1, iprob = NULL, zero = NULL) {

  if (is.finite(upper.limit) &&
      !is.Numeric(upper.limit, integer.valued = TRUE,
                  positive = TRUE))
    stop("bad input for argument 'upper.limit'")

  if (any(upper.limit < 0))
    stop("bad input for argument 'upper.limit'")



  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")



  uu.ll <- min(upper.limit)


  new("vglmff",
  blurb = c("Truncated geometric distribution ",
            "(P[Y=y] =\n",
            "     ",
            "prob * (1 - prob)^y / [1-(1-prob)^",
             uu.ll+1, "], y = 0,1,...,",
             uu.ll, ")\n",
            "Link:     ",
            namesof("prob", link, earg = earg), "\n",
            "Mean:     mu = 1 / prob - 1 ",
            ifelse(is.finite(upper.limit),
                   paste("- (", upper.limit+1, ") * (1 - prob)^",
                         upper.limit+1, " / (1 - ",
                         "(1 - prob)^", upper.limit+1, ")", sep = ""),
                         "")),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = .expected ,
         imethod = .imethod ,
         multipleResponses = TRUE,
         parameters.names = c("prob"),
         upper.limit = .upper.limit ,
         zero = .zero )
  }, list( .zero = zero,
           .expected = expected,
           .imethod = imethod,
           .upper.limit = upper.limit ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly
    extra$upper.limit <- matrix( .upper.limit , n, ncoly, byrow = TRUE)

    if (any(y > extra$upper.limit))
      stop("some response values greater than argument 'upper.limit'")


    mynames1 <- param.names("prob", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 2)
                      1 / (1 + y + 1/16) else
                  if ( .imethod == 3)
                      1 / (1 + apply(y, 2, median) + 1/16) else
                      1 / (1 + colSums(y * w) / colSums(w) + 1/16)

      if (!is.matrix(prob.init))
        prob.init <- matrix(prob.init, n, M, byrow = TRUE)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, M, byrow = TRUE)


        etastart <- theta2eta(prob.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg,
            .upper.limit = upper.limit,
            .imethod = imethod, .iprob = iprob ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob <- eta2theta(eta, .link , earg = .earg )
    QQQ <- 1 - prob
    upper.limit <- extra$upper.limit
    tmp1 <- QQQ^(upper.limit+1)
    answer <- 1 / prob - 1 - (upper.limit+1) * tmp1 / (1 - tmp1)
    answer[!is.finite(answer)] <- 1 / prob[!is.finite(answer)] - 1
    answer
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep( .link , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .earg
    }

    misc$M1 <- M1
    misc$multipleResponses <- TRUE
    misc$expected <- .expected
    misc$imethod <- .imethod
    misc$iprob <- .iprob
  }), list( .link = link, .earg = earg,
            .iprob = iprob,
            .upper.limit = upper.limit,
            .expected = expected, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    prob <- eta2theta(eta, .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      upper.limit <- extra$upper.limit
      ll.elts <- c(w) * (dgeom(x = y, prob = prob, log = TRUE) -
                         log1p(-(1.0 - prob)^(1 + upper.limit)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("truncgeometric"),
  deriv = eval(substitute(expression({
    prob <- eta2theta(eta, .link , earg = .earg )
    sss <- upper.limit <- extra$upper.limit  # Is a matrix

    QQQ <- 1 - prob
    tmp1 <- QQQ^(upper.limit + 1)
    dl.dprob <- 1 / prob  + (0 - y) / (1 - prob) -
                (1+upper.limit) * QQQ^(upper.limit - 0) / (1 - tmp1)
    dl.dprob[!is.finite(upper.limit)] <-  1 / prob[!is.finite(upper.limit)] +
      (0 - y[!is.finite(upper.limit)]) / (1 - prob[!is.finite(upper.limit)])


    dprobdeta <- dtheta.deta(prob, .link , earg = .earg )
    c(w) * cbind(dl.dprob * dprobdeta)
  }), list( .link = link, .earg = earg,
            .upper.limit = upper.limit,
            .expected = expected ))),
  weight = eval(substitute(expression({

    eim.oim.fun <- function(mu.y, sss)
      ifelse(is.finite(sss),
             1/prob^2 + (0 + mu.y) / QQQ^2 - (1+sss) *
             ((sss-0) * QQQ^(sss-1) / (1 - tmp1) +
             (1+sss) * QQQ^(2*sss) / (1 - tmp1)^2),
             1 / (prob^2 * (1 - prob)))


    ned2l.dprob2 <- if ( .expected ) {
      eim.oim.fun(mu, sss)
    } else {
      eim.oim.fun(y, sss)
    }
    wz <- ned2l.dprob2 * dprobdeta^2
    if ( !( .expected ))
      wz <- wz - dl.dprob * d2theta.deta2(prob, .link , earg = .earg )
    c(w) * wz
  }), list( .link = link, .earg = earg,
            .expected = expected ))))
}




