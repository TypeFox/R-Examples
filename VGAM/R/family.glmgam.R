# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












 binomialff <-
  function(link = "logit",
           dispersion = 1,
           multiple.responses = FALSE, onedpar = !multiple.responses,
           parallel = FALSE,  # apply.parint = FALSE,
           zero = NULL,
           bred = FALSE,
           earg.link = FALSE) {


 if (!is.logical(bred) || length(bred) > 1)
   stop("argument 'bred' must be a single logical")



  apply.parint <- FALSE
  estimated.dispersion <- dispersion == 0





  if (earg.link) {
    earg <- link
  } else {
    link <- as.list(substitute(link))
    earg <- link2list(link)
  }
  link <- attr(earg, "function.name")


  ans <-
  new("vglmff",
  blurb = if (multiple.responses) c("Multiple binomial model\n\n", 
         "Link:     ", namesof("mu[,j]", link, earg = earg), "\n",
         "Variance: mu[,j]*(1-mu[,j])") else
         c("Binomial model\n\n", 
         "Link:     ", namesof("prob", link, earg = earg), "\n",
         "Variance: mu * (1 - mu)"),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x, 
                           bool = .parallel , 
                           constraints = constraints,
                           apply.int = .apply.parint )

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero,
            .parallel = parallel, .apply.parint = apply.parint ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         bred = .bred ,
         expected = TRUE,
         parameters.names = c("prob"),  # new.name
         zero = .zero )
  }, list( .zero = zero,
           .bred = bred ))),

  initialize = eval(substitute(expression({
    assign("CQO.FastAlgorithm",
           ( .link == "logit" || .link == "cloglog"),
           envir = VGAMenv)
    assign("modelno", if ( .link == "logit") 1 else
                      if ( .link == "cloglog") 4 else NULL,
           envir = VGAMenv)



    old.name <- "mu"
    new.name <- "prob"



    if ( .multiple.responses ) {
      temp5 <-
      w.y.check(w = w, y = y,
                Is.nonnegative.y = TRUE,
                ncol.w.max = Inf,
                ncol.y.max = Inf,
                out.wy = TRUE,
                colsyperw = 1,
                maximize = TRUE)
      w <- temp5$w
      y <- temp5$y


      y.counts <- y
      y <- y / w





      M <- ncol(y)

  if (FALSE)
      if (!all(y == 0 | y == 1))
        stop("response must contain 0s and 1s only")



      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2)) {
        paste("E[", dn2, "]", sep = "") 
      } else {
        paste(new.name, 1:M, sep = "") 
      }
      predictors.names <-
          namesof(if (M > 1) dn2 else new.name,
                  .link , earg = .earg , short = TRUE)

      if (!length(mustart) && !length(etastart))
        mustart <- matrix(colMeans(y.counts), nrow = nrow(y), ncol = ncol(y),
                         byrow = TRUE) /
                   matrix(colMeans(w), nrow = nrow(w), ncol = ncol(w),
                         byrow = TRUE)



      extra$multiple.responses <- TRUE

    } else {

      if (!all(w == 1))
          extra$orig.w <- w


      NCOL <- function (x) if (is.array(x) && length(dim(x)) > 1 ||
                          is.data.frame(x)) ncol(x) else as.integer(1)
      if (NCOL(y) == 1) {
        if (is.factor(y))
          y <- (y != levels(y)[1])
        nvec <- rep(1, n)
        y[w == 0] <- 0
        if (!all(y == 0 | y == 1))
          stop("response values 'y' must be 0 or 1")
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + w * y) / (1 + w)


        no.successes <- y
        if (min(y) < 0)
          stop("Negative data not allowed!")
        if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
          stop("Number of successes must be integer-valued")
      } else if (NCOL(y) == 2) {
        if (min(y) < 0)
          stop("Negative data not allowed!")
        if (any(abs(y - round(y)) > 1.0e-8))
          stop("Count data must be integer-valued")
        y <- round(y)
        nvec <- y[, 1] + y[, 2]
        y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
        w <- w * nvec
        if (!length(mustart) && !length(etastart))
          mustart <- (0.5 + nvec * y) / (1 + nvec)
        } else {
          stop("for the binomialff family, response 'y' must be a ",
               "vector of 0 and 1's\n",
               "or a factor (first level = fail, other levels = success),\n",
               "or a 2-column matrix where col 1 is the no. of ",
               "successes and col 2 is the no. of failures")
        }
        predictors.names <-
          namesof(new.name, .link , earg = .earg , short = TRUE)
    }


    if ( .bred ) {
      if ( !control$save.weights ) {
       save.weights <- control$save.weights <- TRUE
      }
    }



    }), list( .link = link, .multiple.responses = multiple.responses,
              .earg = earg, .bred = bred ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu <-  eta2theta(eta, link = .link , earg = .earg )


    colnames(mu) <- NULL

    mu
  }, list( .link = link, .earg = earg  ))),

  last = eval(substitute(expression({
    if (exists("CQO.FastAlgorithm", envir = VGAMenv))
      rm("CQO.FastAlgorithm", envir = VGAMenv)
    if (exists("modelno", envir = VGAMenv))
      rm("modelno", envir = VGAMenv)

    dpar <- .dispersion
    if (!dpar) {
      temp87 <- (y-mu)^2 * wz / (
                dtheta.deta(mu, link = .link ,
                            earg = .earg )^2)  # w cancel
      if (.multiple.responses && ! .onedpar ) {
        dpar <- rep(NA_real_, len = M)
        temp87 <- cbind(temp87)
        nrow.mu <- if (is.matrix(mu)) nrow(mu) else length(mu)
        for (ii in 1:M)
          dpar[ii] <- sum(temp87[, ii]) / (nrow.mu - ncol(x))
        if (is.matrix(y) && length(dimnames(y)[[2]]) == length(dpar))
          names(dpar) <- dimnames(y)[[2]]
      } else {
        dpar <- sum(temp87) / (length(mu) - ncol(x))
      }
    }

    misc$multiple.responses <- .multiple.responses
    misc$dispersion <- dpar
    misc$default.dispersion <- 1
    misc$estimated.dispersion <- .estimated.dispersion
    misc$bred <- .bred
    misc$expected <- TRUE

    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else new.name  # Was old.name=="mu"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

  }), list( .dispersion = dispersion,
            .estimated.dispersion = estimated.dispersion,
            .onedpar = onedpar, .multiple.responses = multiple.responses,
            .bred = bred,
            .link = link, .earg = earg))),

  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, .link , earg = .earg )
  }, list( .link = link, .earg = earg))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      c(w) * (y / mu - (1-y) / (1-mu))
    } else {

      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                 y * w  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                round(w)


      smallno <- 1.0e6 * .Machine$double.eps
      smallno <- sqrt(.Machine$double.eps)
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <- if ( .multiple.responses ) {
        c(w) * (    ycounts  * log(   mu) +
               (1 - ycounts) * log1p(-mu))
      } else {
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
        dbinom(x = ycounts, size = nvec, prob = mu, log = TRUE)
      }
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .multiple.responses = multiple.responses ))),

  vfamily = c("binomialff", "VGAMcategorical"),



  simslot = function (object, nsim) {

    ftd <- fitted(object)


    if (ncol(ftd) > 1)
      stop("simulate() does not work with more than one response")



    n <- length(ftd)
    ntot <- n * nsim
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")



    if (any(pwts %% 1 != 0))
      stop("cannot simulate from non-integer prior.weights")
    if (length(m <- object@model) > 0) {
      y <- model.response(m)
      if (is.factor(y)) {
        yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
                     labels = levels(y))
        split(yy, rep(seq_len(nsim), each = n))
      } else if (is.matrix(y) && ncol(y) == 2) {
        yy <- vector("list", nsim)
        for (i in seq_len(nsim)) {
          Y <- rbinom(n, size = pwts, prob = ftd)
          YY <- cbind(Y, pwts - Y)
          colnames(YY) <- colnames(y)
          yy[[i]] <- YY
        }
        yy
      } else {
        rbinom(ntot, size = pwts, prob = ftd)/pwts
      }
    } else {

      rbinom(ntot, size = c(pwts), prob = c(ftd))/c(pwts)
    }
  },




  deriv = eval(substitute(expression({
    yBRED <- if ( .bred ) {
      Hvector <- hatvaluesbasic(X.vlm = X.vlm.save,
                                diagWm = c(t(w * mu)))  # Handles M>1

      varY <- mu * (1 - mu) / w  # Is a matrix if M>1. Seems the most correct.
      d1.ADJ <-   dtheta.deta(mu, .link , earg = .earg )

      temp.earg <- .earg
      temp.earg$inverse <- FALSE
      temp.earg$inverse <- TRUE
      d2.ADJ <- d2theta.deta2(mu, .link , earg = temp.earg )



      yBRED <- y + matrix(Hvector, n, M, byrow = TRUE) *
                   varY * d2.ADJ / (2 * d1.ADJ^2)
      yBRED
    } else {
      y
    }


    answer <- if ( .link == "logit") {
      c(w) * (yBRED - mu)
    } else if ( .link == "cloglog") {
      mu.use <- mu
      smallno <- 100 * .Machine$double.eps
      mu.use[mu.use <       smallno] <-       smallno
      mu.use[mu.use > 1.0 - smallno] <- 1.0 - smallno
      -c(w) * (yBRED - mu) * log1p(-mu.use) / mu.use
    } else {
      c(w) * dtheta.deta(mu, link = .link , earg = .earg ) *
             (yBRED / mu - 1.0) / (1.0 - mu)
    }

    answer
  }), list( .link = link, .earg = earg, .bred = bred))),

  weight = eval(substitute(expression({
    tmp100 <- mu * (1.0 - mu)

    tmp200 <- if ( .link == "logit") {
      cbind(c(w) * tmp100)
    } else if ( .link == "cloglog") {
      cbind(c(w) * (1.0 - mu.use) * (log1p(-mu.use))^2 / mu.use)
    } else {
      cbind(c(w) * dtheta.deta(mu, link = .link ,
                               earg = .earg )^2 / tmp100)
    }
    for (ii in 1:M) {
      index500 <- !is.finite(tmp200[, ii]) |
                  (abs(tmp200[, ii]) < .Machine$double.eps)
      if (any(index500)) {  # Diagonal 0s are bad
        tmp200[index500, ii] <- .Machine$double.eps
      }
    }
    tmp200
  }), list( .link = link, .earg = earg))))






    ans@deviance <- 
      if (multiple.responses)
        function(mu, y, w, residuals = FALSE, eta, extra = NULL,
                 summation = TRUE) {
      Deviance.categorical.data.vgam(mu  = mu,
                                     y   = y,
                                     w   = w, residuals = residuals,
                                     eta = eta, extra = extra,
                                     summation = summation)
        } else
        function(mu, y, w, residuals = FALSE, eta, extra = NULL,
                 summation = TRUE) {
      Deviance.categorical.data.vgam(mu  = cbind(mu, 1-mu),
                                     y   = cbind(y , 1-y),
                                     w   = w, residuals = residuals,
                                     eta = eta, extra = extra,
                                     summation = summation)
        }


  ans
}



 gammaff <- function(link = "nreciprocal", dispersion = 0) {
  estimated.dispersion <- dispersion == 0


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Gamma distribution\n\n",
            "Link:     ", namesof("mu", link, earg = earg), "\n",
            "Variance: mu^2 / k"),
  deviance =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    devi <- -2 * c(w) * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
    if (residuals) {
      sign(y - mu) * sqrt(abs(devi) * w)
    } else {
      dev.elts <- c(w) * devi
      if (summation) {
        sum(dev.elts)
      } else {
        dev.elts
      }
    }
  },
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("mu"),
         dispersion = .dispersion )
  }, list( .dispersion = dispersion ))),
  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    mustart <- y + 0.167 * (y == 0)

    M <- if (is.matrix(y)) ncol(y) else 1
    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) {
      paste("E[", dn2, "]", sep = "") 
    } else {
      paste("mu", 1:M, sep = "") 
    }

    predictors.names <-
      namesof(if (M > 1) dn2 else "mu", .link ,
              earg = .earg , short = TRUE)

    if (!length(etastart))
      etastart <- theta2eta(mustart, link = .link , earg = .earg )
  }), list( .link = link, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, link = .link , earg = .earg )
  }, list( .link = link, .earg = earg))),
  last = eval(substitute(expression({
    dpar <- .dispersion
    if (!dpar) {
      if (M == 1) {
        temp <- c(w) * dmu.deta^2
        dpar <- sum(c(w) * (y-mu)^2 * wz / temp) / (length(mu) - ncol(x))
      } else {
        dpar <- rep(0, len = M)
        for (spp in 1:M) {
          temp <- c(w) * dmu.deta[, spp]^2
          dpar[spp] <- sum(c(w) * (y[,spp]-mu[, spp])^2 * wz[, spp]/temp) / (
                       length(mu[,spp]) - ncol(x))
        }
      }
    }
    misc$dispersion <- dpar
    misc$default.dispersion <- 0
    misc$estimated.dispersion <- .estimated.dispersion

    misc$link <- rep( .link , length = M)
    names(misc$link) <- param.names("mu", M)

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .dispersion = dispersion, .earg = earg,
            .estimated.dispersion = estimated.dispersion,
            .link = link ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link = .link , earg = .earg )
  }, list( .link = link, .earg = earg))),
  vfamily = "gammaff",
  deriv = eval(substitute(expression({
    M1 <- 1
    ncoly <- ncol(as.matrix(y))

    dl.dmu <- (y-mu) / mu^2
    dmu.deta <- dtheta.deta(theta = mu, link = .link , earg = .earg )
    c(w) * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
    d2l.dmu2 <- 1 / mu^2
    wz <- dmu.deta^2 * d2l.dmu2
    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .link = link, .earg = earg))))
}




 inverse.gaussianff <- function(link = "natural.ig",
                                dispersion = 0) {
  estimated.dispersion <- dispersion == 0
  warning("@deviance() not finished")
  warning("needs checking, but I'm sure it works")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Inverse Gaussian distribution\n\n",
            "Link:     ", namesof("mu", link, earg = earg), "\n",
            "Variance: mu^3 / k"),

  deviance =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    pow <- 3  # Use Quasi()$deviance with pow==3
    devy  <- y^(2-pow) / (1-pow) - y^(2-pow) / (2-pow)
    devmu <- y * mu^(1-pow) / (1-pow) - mu^(2-pow) / (2-pow)
    devi <- 2 * (devy - devmu)
    if (residuals) {
      sign(y - mu) * sqrt(abs(devi) * w)
    } else {
      dev.elts <- c(w) * devi
      if (summation) {
        sum(dev.elts)
      } else {
        dev.elts
      }
    }
  },

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("mu"),
         dispersion = .dispersion )
  }, list( .earg = earg , .dispersion = dispersion ))),
  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    mu <- y + 0.167 * (y == 0)



    M <- if (is.matrix(y)) ncol(y) else 1
    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) {
      paste("E[", dn2, "]", sep = "") 
    } else {
      paste("mu", 1:M, sep = "") 
    }

    predictors.names <-
      namesof(if (M > 1) dn2 else "mu", .link , .earg , short = TRUE)


    if (!length(etastart))
      etastart <- theta2eta(mu, link = .link , .earg )
  }), list( .link = link, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, link = .link , earg = .earg )
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    dpar <- .dispersion
    if (!dpar) {
      temp <- c(w) * dmu.deta^2
      dpar <- sum( c(w) * (y-mu)^2 * wz / temp ) / (length(mu) - ncol(x))
    }
    misc$dispersion <- dpar
    misc$default.dispersion <- 0
    misc$estimated.dispersion <- .estimated.dispersion

    misc$link <- rep( .link , length = M)
    names(misc$link) <- param.names("mu", M)

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .dispersion = dispersion,
            .estimated.dispersion = estimated.dispersion,
            .link = link, .earg = earg ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
      theta2eta(mu, link = .link, earg = .earg )
  }, list( .link = link, .earg = earg ))),
  vfamily = "inverse.gaussianff",
  deriv = eval(substitute(expression({
    M1 <- 1
    ncoly <- ncol(as.matrix(y))

    dl.dmu <- (y - mu) / mu^3
    dmu.deta <- dtheta.deta(theta = mu, link = .link , earg = .earg )
    c(w) * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    d2l.dmu2 <- 1 / mu^3
    wz <- dmu.deta^2 * d2l.dmu2
    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .link = link, .earg = earg ))))
}




dinv.gaussian <- function(x, mu, lambda, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(mu), length(lambda))
  x      <- rep(x,      len = LLL);
  mu     <- rep(mu,     len = LLL);
  lambda <- rep(lambda, len = LLL)
  logdensity <- rep(log(0), len = LLL)

  xok <- (x > 0)
  logdensity[xok] = 0.5 * log(lambda[xok] / (2 * pi * x[xok]^3)) -
                    lambda[xok] *
                    (x[xok]-mu[xok])^2 / (2*mu[xok]^2 * x[xok])
  logdensity[mu     <= 0] <- NaN
  logdensity[lambda <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


pinv.gaussian <- function(q, mu, lambda) {
  if (any(mu  <= 0))
    stop("mu must be positive")
  if (any(lambda  <= 0))
    stop("lambda must be positive")

  LLL <- max(length(q), length(mu), length(lambda))
  q      <- rep(q,      len = LLL)
  mu     <- rep(mu,     len = LLL)
  lambda <- rep(lambda, len = LLL)
  ans <- q

  ans[q <= 0] <- 0
  bb <- q > 0
  ans[bb] <- pnorm( sqrt(lambda[bb]/q[bb]) * (q[bb]/mu[bb] - 1)) +
             exp(2*lambda[bb]/mu[bb]) *
             pnorm(-sqrt(lambda[bb]/q[bb]) * (q[bb]/mu[bb] + 1))
  ans
}


rinv.gaussian <- function(n, mu, lambda) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  mu     <- rep(mu,     len = use.n);
  lambda <- rep(lambda, len = use.n)

  u <- runif(use.n)
  Z <- rnorm(use.n)^2 # rchisq(use.n, df = 1)
  phi <- lambda / mu
  y1 <- 1 - 0.5 * (sqrt(Z^2 + 4*phi*Z) - Z) / phi
  ans <- mu * ifelse((1+y1)*u > 1, 1/y1, y1)
  ans[mu     <= 0] <- NaN
  ans[lambda <= 0] <- NaN
  ans
}











 inv.gaussianff <- function(lmu = "loge", llambda = "loge",
                            imethod = 1,  ilambda = NULL,
                            parallel = FALSE,
                            ishrinkage = 0.99,
                            zero = NULL) {


  apply.parint <- FALSE


  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")



  if (is.logical(parallel) && parallel && length(zero))
    stop("set 'zero = NULL' if 'parallel = TRUE'")




  new("vglmff",
  blurb = c("Inverse Gaussian distribution\n\n",
            "f(y) = sqrt(lambda/(2*pi*y^3)) * ",
            "exp(-lambda * (y - mu)^2 / (2 * mu^2 * y)); y, mu & lambda > 0",
            "Link:     ", namesof("mu",     lmu,     earg = emu), ", ",
                          namesof("lambda", llambda, earg = elambda), "\n",
            "Mean:     ", "mu\n",
            "Variance: mu^3 / lambda"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x, 
                           bool = .parallel , 
                           constraints = constraints,
                           apply.int = .apply.parint )

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero,
            .parallel = parallel, .apply.parint = apply.parint ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         parameters.names = c("mu", "lambda"),
         expected = TRUE,
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
    M1 <- 2
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly



    mynames1 <- param.names("mu",     ncoly)
    mynames2 <- param.names("lambda", ncoly)
    predictors.names <-
      c(namesof(mynames1, .lmu ,     earg = .emu ,     short = TRUE),
        namesof(mynames2, .llambda , earg = .elambda , short = TRUE))[
          interleave.VGAM(M, M1 = M1)]




    if (!length(etastart)) {
      init.mu <-
        if ( .imethod == 2) {
          mediany <- apply(y, 2, median)
          matrix(1.1 * mediany + 1/8, n, ncoly, byrow = TRUE)
        } else if ( .imethod == 3) {
          use.this <- colSums(y * w) / colSums(w)  # weighted.mean(y, w)
          (1 - .ishrinkage ) * y  + .ishrinkage * use.this
        } else {
          matrix(colSums(y * w) / colSums(w) + 1/8,
                 n, ncoly, byrow = TRUE)
        }

      variancey <- apply(y, 2, var)
      init.la <- matrix(if (length( .ilambda )) .ilambda else
                        (init.mu^3) / (0.10 + variancey),
                        n, ncoly, byrow = TRUE)

      etastart <- cbind(
          theta2eta(init.mu, link = .lmu , earg = .emu ),
          theta2eta(init.la, link = .llambda , earg = .elambda ))[,
          interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda,
            .ishrinkage = ishrinkage,
            .imethod = imethod, .ilambda = ilambda ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, c(TRUE, FALSE)], link = .lmu , earg = .emu )
  }, list( .lmu = lmu, .emu = emu, .elambda = elambda ))),

  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lmu ,     length = ncoly),
        rep( .llambda , length = ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .emu
      misc$earg[[M1*ii  ]] <- .elambda
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$ishrinkage <- .ishrinkage 
    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
    misc$parallel <- .parallel
    misc$apply.parint <- .apply.parint
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda,
            .parallel = parallel, .apply.parint = apply.parint,
            .ishrinkage = ishrinkage,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    mymu   <- eta2theta(eta[, c(TRUE, FALSE)],
                        link = .lmu , earg = .emu )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)],
                        link = .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dinv.gaussian(x = y, mu = mymu,
                                      lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .llambda = llambda,
           .emu = emu, .elambda = elambda ))),

  vfamily = "inv.gaussianff",

  deriv = eval(substitute(expression({
    M1 <- 2
    mymu   <- eta2theta(eta[, c(TRUE, FALSE)],
                        link = .lmu ,     earg = .emu )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)],
                        link = .llambda , earg = .elambda )

    dmu.deta <- dtheta.deta(theta = mymu , link = .lmu , earg = .emu )
    dlambda.deta <- dtheta.deta(theta = lambda, link = .llambda ,
                                earg = .elambda )

    dl.dmu <- lambda * (y - mymu) / mymu^3
    dl.dlambda <- 0.5 / lambda - (y - mymu)^2 / (2 * mymu^2 * y)
    myderiv <- c(w) * cbind(dl.dmu * dmu.deta,
                            dl.dlambda * dlambda.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda ))),

  weight = eval(substitute(expression({

    ned2l.dmu2 <- lambda / mymu^3
    ned2l.dlambda2 <- 0.5 / (lambda^2)

    wz <- cbind(dmu.deta^2 * ned2l.dmu2,
                dlambda.deta^2 * ned2l.dlambda2)[,
                interleave.VGAM(M, M1 = M1)]

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)
  }), list( .lmu = lmu, .llambda = llambda,
            .emu = emu, .elambda = elambda ))))
}




 poissonff <- function(link = "loge",
                       dispersion = 1, onedpar = FALSE,
                       imu = NULL, imethod = 1,
                       parallel = FALSE, zero = NULL,
                       bred = FALSE,
                       earg.link = FALSE) {



  if (!is.logical(bred) || length(bred) > 1)
    stop("argument 'bred' must be a single logical")

  estimated.dispersion <- (dispersion == 0)


  if (earg.link) {
    earg <- link
  } else {
    link <- as.list(substitute(link))
    earg <- link2list(link)
  }
  link <- attr(earg, "function.name")




  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (length(imu) &&
      !is.Numeric(imu, positive = TRUE))
    stop("bad input for argument 'imu'")


  new("vglmff",
  blurb = c("Poisson distribution\n\n",
            "Link:     ", namesof("lambda", link, earg = earg), "\n",
            "Variance: lambda"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x, 
                           bool = .parallel , 
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel = parallel, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("lambda"),
         bred = .bred ,
         zero = .zero )
  }, list( .zero = zero,
           .bred = bred ))),


  deviance =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    nz <- (y > 0)
    devi <-  -(y - mu)
    devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
    if (residuals) {
      sign(y - mu) * sqrt(2 * abs(devi) * c(w))
    } else {
      dev.elts <- 2 * c(w) * devi
      if (summation) {
        sum(dev.elts)
      } else {
        dev.elts
      }
    }
  },

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    M <- ncoly <- ncol(y)

    assign("CQO.FastAlgorithm", ( .link == "loge"), envir = VGAMenv)



    old.name <- "mu"
    new.name <- "lambda"
    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) {
      paste("E[", dn2, "]", sep = "") 
    } else {
      paste(new.name, 1:M, sep = "") 
    }
    predictors.names <-
      namesof(if (M > 1) dn2 else new.name, # was "mu" == old.name
              .link ,
              earg = .earg , short = TRUE)


    if ( .bred ) {
      if ( !control$save.weights ) {
       save.weights <- control$save.weights <- TRUE
      }
    }



    if (!length(etastart)) {
      mu.init <- pmax(y, 1/8)
      for (iii in 1:ncol(y)) {
        if ( .imethod == 2) {
          mu.init[, iii] <- weighted.mean(y[, iii], w[, iii]) + 1/8
        } else if ( .imethod == 3) {
          mu.init[, iii] <- median(y[, iii]) + 1/8
        }
      }
      if (length( .imu ))
        mu.init <- matrix( .imu , n, ncoly, byrow = TRUE)
      etastart <- theta2eta(mu.init, link = .link , earg = .earg )
    }
  }), list( .link = link, .estimated.dispersion = estimated.dispersion,
            .bred = bred,
            .imethod = imethod, .imu = imu, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu <- eta2theta(eta, link = .link , earg = .earg )
    mu
  }, list( .link = link, .earg = earg))),

  last = eval(substitute(expression({
    if (exists("CQO.FastAlgorithm", envir = VGAMenv))
      rm("CQO.FastAlgorithm", envir = VGAMenv)
    dpar <- .dispersion
    if (!dpar) {
      temp87 <- (y-mu)^2 *
          wz / (dtheta.deta(mu, link = .link , earg = .earg )^2)  # w cancel
      if (M > 1 && ! .onedpar ) {
        dpar <- rep(NA_real_, length = M)
        temp87 <- cbind(temp87)
        nrow.mu <- if (is.matrix(mu)) nrow(mu) else length(mu)
        for (ii in 1:M)
          dpar[ii] <- sum(temp87[, ii]) / (nrow.mu - ncol(x))
        if (is.matrix(y) && length(dimnames(y)[[2]]) == length(dpar))
          names(dpar) <- dimnames(y)[[2]]
      } else {
        dpar <- sum(temp87) / (length(mu) - ncol(x))
      }
    }
    misc$dispersion <- dpar
    misc$default.dispersion <- 1
    misc$estimated.dispersion <- .estimated.dispersion

    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$multipleResponses <- TRUE
    misc$bred <- .bred


    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else new.name  # Was old.name=="mu"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

  }), list( .dispersion = dispersion, .imethod = imethod,
            .estimated.dispersion = estimated.dispersion,
            .bred = bred,
            .onedpar = onedpar, .link = link, .earg = earg))),

  linkfun = eval(substitute( function(mu, extra = NULL) {
    theta2eta(mu, link = .link , earg = .earg )
  }, list( .link = link, .earg = earg))),

  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      c(w) * (y / mu - 1)
    } else {
      ll.elts <- c(w) * dpois(x = y, lambda = mu, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = "poissonff",




  simslot =
    function(object, nsim) {


    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    rpois(nsim * length(ftd), ftd)
  },



  deriv = eval(substitute(expression({
    yBRED <- if ( .bred ) {
      Hvector <- hatvaluesbasic(X.vlm = X.vlm.save,
                                diagWm = c(t(c(w) * mu)))  # Handles M>1


      varY <- mu # Is a matrix if M>1.
      d1.BRED <-   dtheta.deta(mu, .link , earg = .earg )
      d2.BRED <- d2theta.deta2(mu, .link , earg = .earg )
      y + matrix(Hvector, n, M, byrow = TRUE) *
                 varY * d2.BRED / (2 * d1.BRED^2)
    } else {
      y
    }


    answer <- if ( .link == "loge" && (any(mu < .Machine$double.eps))) {
      c(w) * (yBRED - mu)
    } else {
      lambda <- mu
      dl.dlambda <- (yBRED - lambda) / lambda
      dlambda.deta <- dtheta.deta(theta = lambda,
                                  link = .link , earg = .earg )
      c(w) * dl.dlambda * dlambda.deta
    }

    answer
  }), list( .link = link, .earg = earg, .bred = bred))),

  weight = eval(substitute(expression({
    if ( .link == "loge" && (any(mu < .Machine$double.eps))) {
      tmp600 <- mu
      tmp600[tmp600 < .Machine$double.eps] <- .Machine$double.eps
      c(w) * tmp600
    } else {
      ned2l.dlambda2 <- 1 / lambda
      ned2lambda.deta2 <- d2theta.deta2(theta = lambda,
                                        link = .link , earg = .earg )
      c(w) * dlambda.deta^2 * ned2l.dlambda2
    }
  }), list( .link = link, .earg = earg))))
}




 quasibinomialff <-
  function(
           link = "logit",
           multiple.responses = FALSE, onedpar = !multiple.responses,
           parallel = FALSE, zero = NULL) {


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  dispersion <- 0 # Estimated; this is the only difference with binomialff()
  ans <- binomialff(link = earg, earg.link = TRUE,
                    dispersion = dispersion,
                    multiple.responses = multiple.responses,
                    onedpar = onedpar,
                    parallel = parallel, zero = zero)
  ans@vfamily <- "quasibinomialff"
  ans@infos <- eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = .multiple.responses ,
         parameters.names = c("prob"),
         zero = .zero )
  }, list( .zero = zero,
           .multiple.responses = multiple.responses )))

  ans
}






 quasipoissonff <- function(link = "loge", onedpar = FALSE,
                            parallel = FALSE, zero = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  dispersion <- 0 # Estimated; this is the only difference with poissonff()
  ans <- poissonff(link = earg, earg.link = TRUE,
                   dispersion = dispersion, onedpar = onedpar,
                   parallel = parallel, zero = zero)
  ans@vfamily <- "quasipoissonff"
  ans@infos <- eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = TRUE,
         parameters.names = c("lambda"),
         zero = .zero )
  }, list( .zero = zero )))

  ans
}





 double.exppoisson <-
  function(lmean = "loge",
           ldispersion = "logit",
           idispersion = 0.8,
           zero = NULL) {

  if (!is.Numeric(idispersion, positive = TRUE))
    stop("bad input for 'idispersion'")


  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  ldisp <- as.list(substitute(ldispersion))
  edisp <- link2list(ldisp)
  ldisp <- attr(edisp, "function.name")

  idisp <- idispersion


  new("vglmff",
  blurb = c("Double exponential Poisson distribution\n\n",
            "Link:     ",
            namesof("mean",       lmean,       earg = emean), ", ",
            namesof("dispersion", ldisp, earg = edisp), "\n",
            "Mean:     ", "mean\n",
            "Variance: mean / dispersion"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         parameters.names = c("mean", "dispersion"),
         lmean       = .lmean ,
         ldispersion = .ldispersion ,
         zero = .zero )
  }, list( .lmean       = lmean,
           .ldispersion = ldispersion,
           .zero = zero ))),


  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1)


    M <- if (is.matrix(y)) ncol(y) else 1
    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) {
        paste("E[", dn2, "]", sep = "") 
    } else {
        "mu"
    }
    predictors.names <-
      c(namesof(dn2,          link = .lmean, earg = .emean, short = TRUE),
        namesof("dispersion", link = .ldisp, earg = .edisp, short = TRUE))

    init.mu <- pmax(y, 1/8)
    tmp2 <- rep( .idisp , length.out = n)

    if (!length(etastart))
      etastart <-
        cbind(theta2eta(init.mu, link = .lmean , earg = .emean ),
              theta2eta(tmp2,    link = .ldisp , earg = .edisp ))
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp,
            .idisp = idisp ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], link = .lmean, earg = .emean)
  }, list( .lmean = lmean, .emean = emean,
           .ldisp = ldisp, .edisp = edisp ))),
  last = eval(substitute(expression({
    misc$link <-    c(mean = .lmean , dispersion = .ldisp )

    misc$earg <- list(mean = .emean , dispersion = .edisp )

    misc$expected <- TRUE
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
      lambda <- eta2theta(eta[, 1], link = .lmean,
                          earg = .emean )
      Disper <- eta2theta(eta[, 2], link = .ldisp,
                          earg = .edisp )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (0.5 * log(Disper) +
                         Disper*(y-lambda) + Disper*y*log(lambda))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmean = lmean, .emean = emean,
           .ldisp = ldisp, .edisp = edisp ))),
  vfamily = "double.exppoisson",
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], link = .lmean, earg = .emean)
    Disper <- eta2theta(eta[, 2], link = .ldisp,
                        earg = .edisp)

    dl.dlambda <- Disper * (y / lambda - 1)
    dl.dDisper <- y * log(lambda) + y - lambda + 0.5 / Disper

    dlambda.deta <- dtheta.deta(theta = lambda, link = .lmean,
                                earg = .emean)
    dDisper.deta <- dtheta.deta(theta = Disper, link = .ldisp,
                                earg = .edisp)

    c(w) * cbind(dl.dlambda * dlambda.deta,
                 dl.dDisper * dDisper.deta)
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, nrow = n, ncol = 2)  # diagonal
    usethis.lambda <- pmax(lambda, .Machine$double.eps / 10000)
    wz[, iam(1, 1, M)] <- (Disper / usethis.lambda) * dlambda.deta^2
    wz[, iam(2, 2, M)] <- (0.5 / Disper^2) * dDisper.deta^2
    c(w) * wz
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp,
            .edisp = edisp ))))
}



 double.expbinomial <-
  function(lmean = "logit", ldispersion = "logit",
           idispersion = 0.25, zero = "dispersion") {

  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  ldisp <- as.list(substitute(ldispersion))
  edisp <- link2list(ldisp)
  ldisp <- attr(edisp, "function.name")
  idisp <- idispersion


  if (!is.Numeric(idispersion, positive = TRUE))
      stop("bad input for 'idispersion'")


  new("vglmff",
  blurb = c("Double Exponential Binomial distribution\n\n",
            "Link:     ",
            namesof("mean",       lmean, earg = emean), ", ",
            namesof("dispersion", ldisp, earg = edisp), "\n",
            "Mean:     ", "mean\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = NA,
         parameters.names = c("mean", "dispersion"),
         lmean = .lmean ,
         ldisp = .ldisp ,
         multipleResponses = FALSE,
         zero = .zero )
  }, list( .lmean = lmean,
           .zero = zero,
           .ldisp = ldisp ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w


    if (ncol(cbind(w)) != 1)
      stop("'weights' must be a vector or a one-column matrix")

        NCOL <- function (x)
            if (is.array(x) && length(dim(x)) > 1 ||
            is.data.frame(x)) ncol(x) else as.integer(1)

        if (NCOL(y) == 1) {


          if (is.factor(y)) y <- (y != levels(y)[1])
          nvec <- rep(1, n)
          y[w == 0] <- 0
          if (!all(y == 0 | y == 1))
            stop("response values 'y' must be 0 or 1")
          init.mu <- (0.5 + w * y) / (1 + w)


          no.successes <- y
          if (min(y) < 0)
            stop("Negative data not allowed!")
          if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
            stop("Number of successes must be integer-valued")
        } else if (NCOL(y) == 2) {
            if (min(y) < 0)
              stop("Negative data not allowed!")
            if (any(abs(y - round(y)) > 1.0e-8))
              stop("Count data must be integer-valued")
            y <- round(y)
            nvec <- y[, 1] + y[, 2]
            y <- ifelse(nvec > 0, y[, 1] / nvec, 0)
            w <- w * nvec
            init.mu <- (0.5 + nvec * y) / (1 + nvec)
        } else
            stop("for the double.expbinomial family, response 'y' must be a ",
                 "vector of 0 and 1's\n",
                     "or a factor (first level = fail, ",
                     "other levels = success),\n",
                     "or a 2-column matrix where col 1 is the no. of ",
                     "successes and col 2 is the no. of failures")

    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) paste("E[", dn2, "]", sep = "") else "mu"

    predictors.names <-
    c(namesof(dn2,          .lmean , earg = .emean , short = TRUE),
      namesof("dispersion", .ldisp , earg = .edisp , short = TRUE))

    tmp2 <- rep( .idisp , len = n)

    if (!length(etastart))
      etastart <- cbind(theta2eta(init.mu, .lmean, earg = .emean),
                        theta2eta(tmp2,    .ldisp, earg = .edisp))
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp,
            .idisp = idisp ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], link = .lmean , earg = .emean )
  }, list( .lmean = lmean, .emean = emean,
           .ldisp = ldisp, .edisp = edisp ))),
  last = eval(substitute(expression({
    misc$link <-    c("mean" = .lmean, "dispersion" = .ldisp)

    misc$earg <- list( mean  = .emean,  dispersion  = .edisp)

    misc$expected <- TRUE
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    prob   <- eta2theta(eta[, 1], link = .lmean, earg = .emean)
    Disper <- eta2theta(eta[, 2], link = .ldisp, earg = .edisp)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {



      temp1 <- y * log(ifelse(y > 0, y, 1))  # y*log(y)
      temp2 <- (1.0-y) * log1p(ifelse(y < 1, -y, 0))  # (1-y)*log(1-y)


      ll.elts <-
         (0.5 * log(Disper) + w * (y * Disper * log(prob) +
         (1-y) * Disper * log1p(-prob) +
         temp1 * (1-Disper) + temp2 * (1 - Disper)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmean = lmean, .emean = emean,
           .ldisp = ldisp, .edisp = edisp ))),
  vfamily = "double.expbinomial",
  deriv = eval(substitute(expression({
    prob   <- eta2theta(eta[, 1], link = .lmean, earg = .emean)
    Disper <- eta2theta(eta[, 2], link = .ldisp, earg = .edisp)
    temp1 <- y * log(ifelse(y > 0, y, 1))  # y*log(y)
    temp2 <- (1.0-y) * log1p(ifelse(y < 1, -y, 0))  # (1-y)*log(1-y)
    temp3 <- prob * (1.0-prob)
    temp3 <- pmax(temp3, .Machine$double.eps * 10000)

    dl.dprob <- w * Disper * (y - prob) / temp3
    dl.dDisper <- 0.5 / Disper + w * (y * log(prob) + 
                 (1-y)*log1p(-prob) - temp1 - temp2)

    dprob.deta   <- dtheta.deta(theta = prob,   .lmean, earg = .emean)
    dDisper.deta <- dtheta.deta(theta = Disper, .ldisp, earg = .edisp)

    cbind(dl.dprob   * dprob.deta,
          dl.dDisper * dDisper.deta)
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, nrow = n, ncol = 2)  # diagonal
    wz[, iam(1, 1, M)] <- w * (Disper / temp3) * dprob.deta^2
    wz[, iam(2, 2, M)] <- (0.5 / Disper^2) * dDisper.deta^2
    wz
  }), list( .lmean = lmean, .emean = emean,
            .ldisp = ldisp, .edisp = edisp ))))
}










 augbinomial <- function(link = "logit", multiple.responses = FALSE,
                        parallel = TRUE) {

    if (!is.logical(parallel) ||
        length(parallel) != 1 ||
        !parallel)
      warning("Argument 'parallel' should be assigned 'TRUE' only")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = if (multiple.responses)
          c("Augmented multivariate binomial model\n\n", 
         "Link:     ",
         namesof("mu.1[,j]", link, earg = earg), ", ",
         namesof("mu.2[,j]", link, earg = earg),
         "\n",
         "Variance: mu[,j]*(1-mu[,j])") else
         c("Augmented binomial model\n\n", 
         "Link:     ",
         namesof("mu.1[,j]", link, earg = earg), ", ",
         namesof("mu.2[,j]", link, earg = earg),
         "\n",
         "Variance: mu*(1-mu)"),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL,
                      summation = TRUE) {
      Deviance.categorical.data.vgam(mu = cbind(mu, 1-mu), y=cbind(y, 1-y),
                                     w = w, residuals = residuals,
                                     eta = eta, extra = extra,
                                     summation = summation)
  },
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         parameters.names = c("mu.1[,j]", "mu.2[,j]"),
         parallel = .parallel)
  }, list( .parallel = parallel ))),
  initialize = eval(substitute(expression({

    M1 = 2

    if ( .multiple.responses ) {
        y = as.matrix(y)
        M = M1 * ncol(y)
        if (!all(y == 0 | y == 1))
            stop("response must contain 0's and 1's only")
        dn2 = if (is.matrix(y)) dimnames(y)[[2]] else NULL
        dn2 = if (length(dn2)) {
            paste("E[", dn2, "]", sep = "") 
        } else {
            paste("mu", 1:M, sep = "") 
        }
        predictors.names <-
          c(namesof(if (M > 1) dn2 else
                    "mu.1", .link , earg = .earg , short = TRUE),
            namesof(if (M > 1) dn2 else
                    "mu.2", .link , earg = .earg , short = TRUE))
        NOS = M / M1
        predictors.names <-
        predictors.names[interleave.VGAM(M1 * NOS, M1 = M1)]


        if (!length(mustart) && !length(etastart))
          mustart = (0.5 + w * y) / (1 + w)
    } else {

      dn2 = c("mu1.", "mu2.")
      M = M1



        if (!all(w == 1))
          extra$orig.w = w


        NCOL = function (x) if (is.array(x) && length(dim(x)) > 1 ||
                          is.data.frame(x)) ncol(x) else as.integer(1)
        if (NCOL(y) == 1) {
            if (is.factor(y)) y = (y != levels(y)[1])
            nvec = rep(1, n)
            y[w == 0] <- 0
            if (!all(y == 0 | y == 1))
                stop("response values 'y' must be 0 or 1")
            if (!length(mustart) && !length(etastart))
              mustart = (0.5 + w * y) / (1 + w)


              no.successes = y
              if (min(y) < 0)
                  stop("Negative data not allowed!")
              if (any(abs(no.successes - round(no.successes)) > 1.0e-8))
                  stop("Number of successes must be integer-valued")
          } else if (NCOL(y) == 2) {
              if (min(y) < 0)
                  stop("Negative data not allowed!")
              if (any(abs(y - round(y)) > 1.0e-8))
                  stop("Count data must be integer-valued")
              y = round(y)
              nvec = y[, 1] + y[, 2]
              y = ifelse(nvec > 0, y[, 1] / nvec, 0)
              w = w * nvec
              if (!length(mustart) && !length(etastart))
                mustart = (0.5 + nvec * y) / (1 + nvec)
          } else {
              stop("for the binomialff family, response 'y' must be a ",
                   "vector of 0 and 1's\n",
                   "or a factor (first level = fail, ",
                                 "other levels = success),\n",
                   "or a 2-column matrix where col 1 is the no. of ",
                   "successes and col 2 is the no. of failures")
          }
          predictors.names <-
            c(namesof("mu.1", .link , earg = .earg , short = TRUE),
              namesof("mu.2", .link , earg = .earg , short = TRUE))
      }
  }), list( .link = link,
            .multiple.responses = multiple.responses, .earg = earg))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Mdiv2  =  ncol(eta) / 2
    index1 =  2*(1:Mdiv2) - 1
    mu =  eta2theta(eta[, index1],
                    link = .link , earg = .earg )
    mu
  }, list( .link = link, .earg = earg  ))),
  last = eval(substitute(expression({
    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else "mu"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$parallel <- .parallel
    misc$expected <- TRUE
    misc$multiple.responses <- .multiple.responses
  }), list( .link = link,
            .multiple.responses = multiple.responses, .earg = earg,
            .parallel = parallel ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    usualanswer = theta2eta(mu, .link , earg = .earg )
    kronecker(usualanswer, matrix(1, 1, 2))
  }, list( .link = link, .earg = earg))),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      c(w) * (y / mu - (1-y) / (1-mu))
    } else {
      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                y * c(w)  # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                round(w)

      smallno <- 1.0e6 * .Machine$double.eps
      smallno <- sqrt(.Machine$double.eps)
      if (max(abs(ycounts - round(ycounts))) > smallno)
          warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dbinom(x = ycounts, size = nvec, prob = mu, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("augbinomial", "VGAMcategorical"),
  deriv = eval(substitute(expression({
    M1 <- 2
    Mdiv2 <-  M / 2

    NOS <- M / M1

    Konst1 <- 1  # Works with this
    deriv1 <- Konst1 * w *
      if ( .link == "logit") {
          y * (1 - mu)
      } else  {
          stop("this is not programmed in yet")
          dtheta.deta(mu, link = .link , earg = .earg ) *
          (y / mu - 1.0) / (1.0 - mu)
      }
    deriv2 = Konst1 * w *
      if ( .link == "logit") {
         -(1 - y) * mu
      } else  {
          stop("this is not programmed in yet")
          dtheta.deta(mu, link = .link , earg = .earg ) *
          (y / mu - 1.0) / (1.0 - mu)
      }

    myderiv = (cbind(deriv1,
                     deriv2))[, interleave.VGAM(M1 * NOS, M1 = M1)]
    myderiv
  }), list( .link = link, .earg = earg))),
  weight = eval(substitute(expression({
      tmp100 <- mu * (1.0 - mu)

      tmp200 <- if ( .link == "logit") {
          cbind(w * tmp100)
        } else {
          cbind(w * dtheta.deta(mu, link = .link , earg = .earg )^2 / tmp100)
        }

      wk.wt1 <- (Konst1^2) * tmp200 * (1 - mu)
      wk.wt2 <- (Konst1^2) * tmp200 *      mu




    my.wk.wt <- cbind(wk.wt1, wk.wt2)
    my.wk.wt <- my.wk.wt[, interleave.VGAM(M1 * NOS, M1 = M1)]
    my.wk.wt
  }), list( .link = link, .earg = earg))))
}





