# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












process.binomial2.data.VGAM <- expression({



  if (!all(w == 1))
    extra$orig.w <- w


  if (!is.matrix(y)) {
    yf <- as.factor(y)
    lev <- levels(yf)
    llev <- length(lev)
    if (llev != 4)
        stop("response must have 4 levels")
    nn <- length(yf)
    y <- matrix(0, nn, llev)
    y[cbind(1:nn, as.vector(unclass(yf)))] <- 1
    colnamesy <- paste(lev, ":", c("00", "01", "10", "11"), sep = "")
    dimnames(y) <- list(names(yf), colnamesy)
    input.type <- 1
  } else if (ncol(y) == 2) {
    if (!all(y == 0 | y == 1))
      stop("response must contains 0's and 1's only")
    col.index <- y[, 2] + 2*y[, 1] + 1    # 1:4
    nn <- nrow(y)
    y <- matrix(0, nn, 4)
    y[cbind(1:nn, col.index)] <- 1
    dimnames(y) <- list(dimnames(y)[[1]],
                        c("00", "01", "10", "11"))
    input.type <- 2
  } else if (ncol(y) == 4) {
    input.type <- 3
  } else
    stop("response unrecognized")



  nvec <- rowSums(y)

  w <- w * nvec
  y <- y / nvec             # Convert to proportions

  if (length(mustart) + length(etastart) == 0) {
    mu <- y + (1 / ncol(y) - y) / nvec
    dimnames(mu) <- dimnames(y)
    mustart <- mu
  }
})






betabinomial.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 betabinomial <- function(lmu = "logit",
                          lrho = "logit",
                          irho = NULL,
                          imethod = 1, ishrinkage = 0.95,
                          nsimEIM = NULL, zero = "rho") {
  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2, 3 or 4")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
      ishrinkage < 0 ||
      ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (!is.null(nsimEIM)) {
    if (!is.Numeric(nsimEIM, length.arg = 1,
                    integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 10)
      warning("'nsimEIM' should be an integer greater than 10, say")
  }

  new("vglmff",
  blurb = c("Beta-binomial model\n",
            "Links:      ",
            namesof("mu",  lmu,  earg = emu), ", ",
            namesof("rho", lrho, earg = erho), "\n",
            "Mean:       mu", "\n",
            "Variance:   mu*(1-mu)*(1+(w-1)*rho)/w"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu", "rho"),
         imethod  = .imethod ,
         ishrinkage  = .ishrinkage ,
         nsimEIM  = .nsimEIM ,
         lmu  = .lmu ,
         lrho = .lrho ,
         zero = .zero )
  }, list( .lmu = lmu, .lrho = lrho,
           .imethod = imethod, .ishrinkage = ishrinkage,
           .zero = zero,
           .nsimEIM = nsimEIM ))),


  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w

    if (is.null( .nsimEIM )) {
      save.weights <- control$save.weights <- FALSE
    }

    mustart.orig <- mustart
    eval(binomialff()@initialize)  # Note: n,w,y,mustart is changed 
    if (length(mustart.orig))
      mustart <- mustart.orig  # Retain it if inputted


    ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
               y * w  # Convert proportions to counts
    if (max(abs(ycounts - round(ycounts))) > 1.0e-6)
      warning("the response (as counts) does not appear to ",
              "be integer-valued. Am rounding to integer values.")
    ycounts <- round(ycounts)  # Make sure it is an integer
    predictors.names <-
      c(namesof("mu",  .lmu ,  earg = .emu  , tag = FALSE),
        namesof("rho", .lrho , earg = .erho , tag = FALSE))

    if (!length(etastart)) {
      betabinomial.Loglikfun <- function(rhoval, y, x, w, extraargs) {
        shape1 <-    extraargs$mustart  * (1-rhoval) / rhoval
        shape2 <- (1-extraargs$mustart) * (1-rhoval) / rhoval
        ycounts <- extraargs$ycounts  # Ought to be integer-valued
        nvec <- extraargs$nvec
        sum(dbetabinom.ab(x = ycounts, size = nvec, shape1 = shape1,
                          shape2 = shape2, log = TRUE))
      }
      rho.grid <- seq(0.05, 0.95, len = 25)  # rvar =
      mustart.use <- if (length(mustart.orig)) {
        mustart.orig
      } else if ( .imethod == 1) {
        rep(weighted.mean(y, w), len = n)
      } else if ( .imethod == 2) {
        .ishrinkage * weighted.mean(y, w) + (1 - .ishrinkage ) * y
      } else if ( .imethod == 3) {
        y.matrix <- cbind(y)
        mat.temp <- matrix(colMeans(y.matrix), nrow(y.matrix),
                           ncol(y.matrix), byrow = TRUE)
        0.5 * mustart + 0.5 * mat.temp
      } else {
        mustart
      }
      try.this <- grid.search(rho.grid, objfun = betabinomial.Loglikfun,
                              y = y,  x = x, w = w,
                              extraargs = list(
                              ycounts = ycounts,
                              nvec = if (is.numeric(extra$orig.w))
                                     round(w / extra$orig.w) else round(w),
                              mustart = mustart.use))
      init.rho <- if (is.Numeric( .irho ))
                    rep( .irho , length = n) else
                    rep(try.this, length = n)
      etastart <-
        cbind(theta2eta(mustart.use,  .lmu ,  earg = .emu ),
              theta2eta(init.rho,     .lrho , earg = .erho ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho,
            .imethod = imethod, .ishrinkage = ishrinkage,
            .nsimEIM = nsimEIM, .irho = irho ))),
  linkinv = eval(substitute(function(eta, extra = NULL)
    eta2theta(eta[, 1], .lmu , earg = .emu ), 
  list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
    misc$link <-    c(mu = .lmu , rho = .lrho)

    misc$earg <- list(mu = .emu , rho = .erho )

    misc$zero <- .zero
    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$rho <- 1 / (shape1 + shape2 + 1)
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho,
            .nsimEIM = nsimEIM, .zero = zero ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
               y * w  # Convert proportions to counts

    mymu <- eta2theta(eta[, 1], .lmu ,  earg = .emu )
    rho  <- eta2theta(eta[, 2], .lrho , earg = .erho )
    smallno <- 1.0e4 * .Machine$double.eps

    if (max(abs(ycounts - round(ycounts))) > smallno)
      warning("converting 'ycounts' to integer in @loglikelihood")
    ycounts <- round(ycounts)

    rho  <- pmax(rho,     smallno)
    rho  <- pmin(rho, 1 - smallno)
    shape1 <-      mymu  * (1 - rho) / rho
    shape2 <- (1 - mymu) * (1 - rho) / rho

    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
            round(w)

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dbetabinom.ab(x = ycounts, size = nvec, shape1 = shape1,
                       shape2 = shape2, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho  ))),
  vfamily = c("betabinomial"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    w <- pwts
    eta <- predict(object)
    extra <- object@extra

    mymu <- eta2theta(eta[, 1], .lmu ,  earg = .emu )
    rho  <- eta2theta(eta[, 2], .lrho , earg = .erho )
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)

    rbetabinom(nsim * length(rho), size = nvec,
               prob = mymu, rho  = rho)
  }, list( .lmu = lmu, .lrho = lrho,
           .emu = emu, .erho = erho  ))),





  deriv = eval(substitute(expression({
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
    ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
              y * w  # Convert proportions to counts

    ycounts <- round(ycounts)
    mymu <- eta2theta(eta[, 1], .lmu ,  earg = .emu )
    rho  <- eta2theta(eta[, 2], .lrho , earg = .erho )
    smallno <- 100 * .Machine$double.eps
    rho  <- pmax(rho, smallno)
    rho  <- pmin(rho, 1-smallno)

    shape1 <-      mymu  * (1 - rho) / rho
    shape2 <- (1 - mymu) * (1 - rho) / rho
    dshape1.dmu <-  (1 - rho) / rho
    dshape2.dmu <- -(1 - rho) / rho
    dshape1.drho <-       -mymu  / rho^2
    dshape2.drho <-  -(1 - mymu) / rho^2

    dmu.deta  <- dtheta.deta(mymu, .lmu  , earg = .emu )
    drho.deta <- dtheta.deta(rho,  .lrho , earg = .erho )

    dl.dmu <- dshape1.dmu * (digamma(shape1+ycounts) -
              digamma(shape2+nvec-ycounts) -
              digamma(shape1) + digamma(shape2))
    dl.drho <- (-1/rho^2) * (mymu * digamma(shape1 + ycounts) +
               (1 - mymu) * digamma(shape2 + nvec - ycounts) -
               digamma(shape1 + shape2 + nvec) -
               mymu * digamma(shape1) -
               (1 - mymu)*digamma(shape2) + digamma(shape1+shape2))

    c(if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
    cbind(dl.dmu * dmu.deta, dl.drho * drho.deta)
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho  ))),
  weight = eval(substitute(expression({
    if (is.null( .nsimEIM )) {
      wz <- matrix(NA_real_, n, dimm(M))  #3=dimm(2)
      wz11 <- -(expected.betabin.ab(nvec, shape1, shape2, TRUE) -
               trigamma(shape1+shape2+nvec) -
               trigamma(shape1) + trigamma(shape1+shape2))
      wz22 <- -(expected.betabin.ab(nvec, shape1, shape2, FALSE) -
               trigamma(shape1+shape2+nvec) -
               trigamma(shape2) + trigamma(shape1+shape2))
      wz21 <- -(trigamma(shape1+shape2) - trigamma(shape1+shape2+nvec))

      wz[, iam(1, 1, M)] <- dmu.deta^2 * (wz11 * dshape1.dmu^2 +
                                      wz22 * dshape2.dmu^2 +
                         2 * wz21 * dshape1.dmu * dshape2.dmu)
      wz[, iam(2, 2, M)] <- drho.deta^2 * (wz11 * dshape1.drho^2 +
                                       wz22 * dshape2.drho^2 +
                         2 * wz21 * dshape1.drho * dshape2.drho)
      wz[, iam(2, 1, M)] <- dmu.deta * drho.deta *
                  (dshape1.dmu*(wz11*dshape1.drho + wz21*dshape2.drho) +
                  dshape2.dmu*(wz21*dshape1.drho + wz22*dshape2.drho))

      wz * (if (is.numeric(extra$orig.w)) extra$orig.w else 1)
    } else {
      run.varcov <- 0
      ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
      dthetas.detas <- cbind(dmu.deta, drho.deta)

      for (ii in 1:( .nsimEIM )) {
        ysim <- rbetabinom.ab(n = n, size = nvec,
                              shape1 = shape1,
                              shape2 = shape2)
        dl.dmu <- dshape1.dmu * (digamma(shape1+ysim) -
                 digamma(shape2+nvec-ysim) -
                 digamma(shape1) + digamma(shape2))
        dl.drho <- (-1/rho^2) * (mymu * digamma(shape1+ysim) +
                  (1-mymu) * digamma(shape2+nvec-ysim) -
                  digamma(shape1+shape2+nvec) - 
                  mymu * digamma(shape1) -
                  (1-mymu)*digamma(shape2) + digamma(shape1+shape2))


        temp3 <- cbind(dl.dmu, dl.drho)  # n x M matrix
        run.varcov <- run.varcov +
                      temp3[, ind1$row.index] * temp3[, ind1$col.index]
      }
      run.varcov <- run.varcov / .nsimEIM


      wz <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz <- wz * dthetas.detas[, ind1$row] *
                 dthetas.detas[, ind1$col]
      wz * (if (is.numeric(extra$orig.w)) extra$orig.w else 1)
    }
  }), list( .lmu = lmu, .lrho = lrho,
            .emu = emu, .erho = erho, 
            .nsimEIM = nsimEIM ))))
}






dbinom2.or <-
  function(mu1,
           mu2 = if (exchangeable) mu1 else
                 stop("'mu2' not specified"),
           oratio = 1,
           exchangeable = FALSE,
           tol = 0.001,
           colnames = c("00", "01", "10", "11"),
           ErrorCheck = TRUE) {
  if (ErrorCheck) {
    if (!is.Numeric(mu1, positive = TRUE) || max(mu1) >= 1)
      stop("bad input for argument 'mu1'") 
    if (!is.Numeric(mu2, positive = TRUE) || max(mu2) >= 1)
      stop("bad input for argument 'mu2'") 
    if (!is.Numeric(oratio, positive = TRUE))
      stop("bad input for argument 'oratio'") 
    if (!is.Numeric(tol, positive = TRUE, length.arg = 1) ||
        tol > 0.1)
      stop("bad input for argument 'tol'") 
    if (exchangeable && max(abs(mu1 - mu2)) > 0.00001)
      stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ")
  }

  L <- max(length(mu1), length(mu2), length(oratio))
  if (length(oratio) != L) oratio <- rep(oratio, len = L)
  if (length(mu1   ) != L) mu1    <- rep(mu1,    len = L)
  if (length(mu2   ) != L) mu2    <- rep(mu2,    len = L)

  a.temp <- 1 + (mu1+mu2)*(oratio-1)
  b.temp <- -4 * oratio * (oratio-1) * mu1 * mu2
  temp <- sqrt(a.temp^2 + b.temp)
  p11 <- ifelse(abs(oratio-1) < tol,
                mu1*mu2,
               (a.temp-temp)/(2*(oratio-1)))
  p01 <- mu2 - p11
  p10 <- mu1 - p11
  p00 <- 1 - p11 - p01 - p10
  matrix(c(p00, p01, p10, p11), L, 4, dimnames = list(NULL, colnames))
}




rbinom2.or <-
  function(n, mu1,
           mu2 = if (exchangeable) mu1 else
                   stop("argument 'mu2' not specified"),
           oratio = 1,
           exchangeable = FALSE,
           tol = 0.001,
           twoCols = TRUE,
           colnames = if (twoCols) c("y1", "y2") else
                         c("00", "01", "10", "11"),
           ErrorCheck = TRUE) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (ErrorCheck) {
    if (!is.Numeric(mu1, positive = TRUE) || max(mu1) >= 1)
      stop("bad input for argument 'mu1'") 
    if (!is.Numeric(mu2, positive = TRUE) || max(mu2) >= 1)
      stop("bad input for argument 'mu2'") 
    if (!is.Numeric(oratio, positive = TRUE))
      stop("bad input for argument 'oratio'") 
    if (!is.Numeric(tol, positive = TRUE, length.arg = 1) ||
        tol > 0.1)
      stop("bad input for argument 'tol'") 
    if (exchangeable && max(abs(mu1 - mu2)) > 0.00001)
      stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ")
  }

  dmat <- dbinom2.or(mu1 = mu1, mu2 = mu2, oratio = oratio,
                     exchangeable = exchangeable,
                     tol = tol, ErrorCheck = ErrorCheck)

  answer <- matrix(0, use.n, 2,
                   dimnames = list(NULL,
                                   if (twoCols) colnames else NULL))
  yy <- runif(use.n)
  cs1 <- dmat[, "00"] + dmat[, "01"]
  cs2 <- cs1 + dmat[, "10"]
  index <- (dmat[, "00"] < yy) & (yy <= cs1)
  answer[index, 2] <- 1
  index <- (cs1 < yy) & (yy <= cs2)
  answer[index, 1] <- 1
  index <- (yy > cs2)
  answer[index,] <- 1
  if (twoCols) {
    answer
  } else {
    answer4 <- matrix(0, use.n, 4, dimnames = list(NULL, colnames))
    answer4[cbind(1:use.n, 1 + 2*answer[, 1] + answer[, 2])] <- 1
    answer4
  }
}




 binom2.or <- function(lmu = "logit", lmu1 = lmu, lmu2 = lmu,
                       loratio = "loge",
                       imu1 = NULL, imu2 = NULL, ioratio = NULL,
                       zero = "oratio",
                       exchangeable = FALSE,
                       tol = 0.001,
                       more.robust = FALSE) {

  lmu1 <- lmu1
  lmu2 <- lmu2


  lmu1 <- as.list(substitute(lmu1))
  emu1 <- link2list(lmu1)
  lmu1 <- attr(emu1, "function.name")

  lmu2 <- as.list(substitute(lmu2))
  emu2 <- link2list(lmu2)
  lmu2 <- attr(emu2, "function.name")


  loratio <- as.list(substitute(loratio))
  eoratio <- link2list(loratio)
  loratio <- attr(eoratio, "function.name")



  if (!is.logical(exchangeable))
    warning("argument 'exchangeable' should be a single logical") 

  if (is.logical(exchangeable) && exchangeable &&
     ((lmu1 != lmu2) || !all.equal(emu1, emu2)))
    warning("exchangeable = TRUE but marginal links are not equal") 



  if (!is.Numeric(tol, positive = TRUE, length.arg = 1) ||
      tol > 0.1)
    stop("bad input for argument 'tol'") 


  new("vglmff",
  blurb = c("Bivariate binomial regression with an odds ratio\n",
            "Links:    ",
            namesof("mu1", lmu1, earg = emu1), ", ",
            namesof("mu2", lmu2, earg = emu2), "; ",
            namesof("oratio", loratio, earg = eoratio)),
  constraints = eval(substitute(expression({
    cm.intercept.default <- diag(3)
    constraints <- cm.VGAM(matrix(c(1, 1, 0, 0, 0, 1), 3, 2), x = x,
                           bool = .exchangeable ,
                           constraints = constraints,
                           apply.int = TRUE,
                           cm.default           = cm.intercept.default,
                           cm.intercept.default = cm.intercept.default)
      constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                  predictors.names = predictors.names,
                                  M1 = 3)
  }), list( .exchangeable = exchangeable, .zero = zero ))),
  deviance = Deviance.categorical.data.vgam,

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu1", "mu2", "oratio"),
         exchangeable = .exchangeable ,
         lmu1 = .lmu1 ,
         lmu2 = .lmu2 ,
         loratio = .loratio ,
         zero = .zero )
  }, list( .lmu1 = lmu1,
           .lmu2 = lmu2,
           .loratio = loratio,
           .zero = zero,
           .exchangeable = exchangeable
         ))),


  initialize = eval(substitute(expression({
    mustart.orig <- mustart
    eval(process.binomial2.data.VGAM)
    if (length(mustart.orig))
      mustart <- mustart.orig  # Retain it if inputted


    predictors.names <-
       c(namesof("mu1",     .lmu1,    earg = .emu1,    short = TRUE),
         namesof("mu2",     .lmu2,    earg = .emu2,    short = TRUE),
         namesof("oratio",  .loratio, earg = .eoratio, short = TRUE))


    if (!length(etastart)) {
        pmargin <- cbind(mustart[, 3] + mustart[, 4],
                         mustart[, 2] + mustart[, 4])
        ioratio <- if (length( .ioratio)) rep( .ioratio , len = n) else
                   mustart[, 4] * mustart[, 1] / (mustart[, 2] *
                                                  mustart[, 3])
        if (length( .imu1 )) pmargin[, 1] <- .imu1
        if (length( .imu2 )) pmargin[, 2] <- .imu2
        etastart <- cbind(theta2eta(pmargin[, 1], .lmu1, earg = .emu1),
                          theta2eta(pmargin[, 2], .lmu2, earg = .emu2),
                          theta2eta(ioratio, .loratio, earg = .eoratio))
    }
  }), list( .lmu1 = lmu1, .lmu2 = lmu2, .loratio = loratio,
            .emu1 = emu1, .emu2 = emu2, .eoratio = eoratio,
            .imu1 = imu1, .imu2 = imu2, .ioratio = ioratio ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pmargin <- cbind(eta2theta(eta[, 1], .lmu1, earg = .emu1),
                     eta2theta(eta[, 2], .lmu2, earg = .emu2))
    oratio <- eta2theta(eta[, 3], .loratio, earg = .eoratio)
    a.temp <- 1 + (pmargin[, 1]+pmargin[, 2])*(oratio-1)
    b.temp <- -4 * oratio * (oratio-1) * pmargin[, 1] * pmargin[, 2]
    temp <- sqrt(a.temp^2 + b.temp)
    pj4 <- ifelse(abs(oratio-1) < .tol, pmargin[, 1]*pmargin[, 2],
                 (a.temp-temp)/(2*(oratio-1)))
    pj2 <- pmargin[, 2] - pj4
    pj3 <- pmargin[, 1] - pj4
    cbind("00" = 1-pj4-pj2-pj3,
          "01" = pj2,
          "10" = pj3,
          "11" = pj4)
  }, list( .lmu1 = lmu1, .lmu2 = lmu2, .loratio = loratio,
           .emu1 = emu1, .emu2 = emu2, .eoratio = eoratio,
           .tol = tol ))),
  last = eval(substitute(expression({
    misc$link <-    c(mu1 = .lmu1 , mu2 = .lmu2 , oratio = .loratio )

    misc$earg <- list(mu1 = .emu1 , mu2 = .emu2 , oratio = .eoratio )

    misc$tol <- .tol
    misc$expected <- TRUE
  }), list( .lmu1 = lmu1, .lmu2 = lmu2, .loratio = loratio,
            .emu1 = emu1, .emu2 = emu2, .eoratio = eoratio,
            .tol = tol ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    pmargin <- cbind(mu[, 3]+mu[, 4], mu[, 2]+mu[, 4])
    oratio <- mu[, 4]*mu[, 1] / (mu[, 2]*mu[, 3])
    cbind(theta2eta(pmargin[, 1], .lmu1 , earg = .emu1),
          theta2eta(pmargin[, 2], .lmu2 , earg = .emu2), 
          theta2eta(oratio,      .loratio, earg = .eoratio))
  }, list( .lmu1 = lmu1, .lmu2 = lmu2, .loratio = loratio,
           .emu1 = emu1, .emu2 = emu2, .eoratio = eoratio ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      if ( .more.robust) {
        vsmallno <-  1.0e4 * .Machine$double.xmin
        mu[mu < vsmallno] <- vsmallno
      }

      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                 y * w # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                 round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .more.robust = more.robust ))),
  vfamily = c("binom2.or", "binom2"),








  deriv = eval(substitute(expression({
    smallno <- 1.0e4 * .Machine$double.eps
    mu.use <- mu
    mu.use[mu.use <     smallno] <-     smallno
    mu.use[mu.use > 1 - smallno] <- 1 - smallno
    pmargin <- cbind(mu.use[, 3] + mu.use[, 4],
                     mu.use[, 2] + mu.use[, 4])
    pmargin[, 1] <- pmax(    smallno, pmargin[, 1])
    pmargin[, 1] <- pmin(1 - smallno, pmargin[, 1])
    pmargin[, 2] <- pmax(    smallno, pmargin[, 2])
    pmargin[, 2] <- pmin(1 - smallno, pmargin[, 2])

    oratio <- mu.use[, 4]*mu.use[, 1] / (mu.use[, 2]*mu.use[, 3])
    use.oratio <- pmax(smallno, oratio)
    a.temp <- 1 + (pmargin[, 1]+pmargin[, 2])*(oratio-1)
    b.temp <- -4 * oratio * (oratio-1) * pmargin[, 1] * pmargin[, 2]
    temp9 <- sqrt(a.temp^2 + b.temp)

    coeff12 <- -0.5 + (2*oratio*pmargin - a.temp) / (2*temp9)
    dl.dmu1 <- coeff12[, 2] * (y[, 1]/mu.use[, 1]-y[, 3]/mu.use[, 3]) -
       (1+coeff12[, 2]) * (y[, 2]/mu.use[, 2]-y[, 4]/mu.use[, 4])
    
    dl.dmu2 <- coeff12[, 1] * (y[, 1]/mu.use[, 1]-y[, 2]/mu.use[, 2]) -
       (1+coeff12[, 1]) * (y[, 3]/mu.use[, 3]-y[, 4]/mu.use[, 4])

    coeff3 <- (y[, 1]/mu.use[, 1] - y[, 2]/mu.use[, 2] -
              y[, 3]/mu.use[, 3] + y[, 4]/mu.use[, 4])
    Vab <- pmax(smallno, 1 / (1/mu.use[, 1] + 1/mu.use[, 2] +
                             1/mu.use[, 3] + 1/mu.use[, 4]))
    dp11.doratio <- Vab / use.oratio
    dl.doratio <- coeff3 * dp11.doratio

    c(w) * cbind(dl.dmu1 * dtheta.deta(pmargin[, 1], .lmu1, earg = .emu1),
                 dl.dmu2 * dtheta.deta(pmargin[, 2], .lmu2, earg = .emu2),
              dl.doratio * dtheta.deta(oratio, .loratio, earg = .eoratio))
  }), list( .lmu1 = lmu1, .lmu2 = lmu2, .loratio = loratio,
            .emu1 = emu1, .emu2 = emu2, .eoratio = eoratio ))),
  weight = eval(substitute(expression({
    Deltapi <- mu.use[, 3]*mu.use[, 2] - mu.use[, 4]*mu.use[, 1]
    myDelta  <- pmax(smallno, mu.use[, 1] * mu.use[, 2] *
                             mu.use[, 3] * mu.use[, 4])
    pqmargin <- pmargin * (1-pmargin)
    pqmargin[pqmargin < smallno] <- smallno

    wz <- matrix(0, n, 4)
    wz[, iam(1, 1, M)] <- (pqmargin[, 2] * Vab / myDelta) *
                      dtheta.deta(pmargin[, 1], .lmu1, earg = .emu1)^2
    wz[, iam(2, 2, M)] <- (pqmargin[, 1] * Vab / myDelta) *
                      dtheta.deta(pmargin[, 2], .lmu2, earg = .emu2)^2
    wz[, iam(3, 3, M)] <- (Vab / use.oratio^2) *
                 dtheta.deta(use.oratio, .loratio, earg = .eoratio)^2
    wz[, iam(1, 2, M)] <- (Vab * Deltapi / myDelta) *
                      dtheta.deta(pmargin[, 1], .lmu1, earg = .emu1) *
                      dtheta.deta(pmargin[, 2], .lmu2, earg = .emu2)
    c(w) * wz
  }), list( .lmu1 = lmu1, .lmu2 = lmu2, .loratio = loratio,
            .emu1 = emu1, .emu2 = emu2, .eoratio = eoratio ))))
}  # binom2.or










setClass("binom2",         contains = "vglmff")
setClass("binom2.or",      contains = "binom2")





setMethod("summaryvglmS4VGAM",  signature(VGAMff = "binom2.or"),
  function(object,
           VGAMff,
           ...) {

  cfit <- coef.vlm(object, matrix = TRUE)
  if (rownames(cfit)[1] == "(Intercept)" &&
      all(cfit[-1, 3] == 0)) {
    object@post$oratio <- eta2theta(cfit[1, 3],
                                    link = object@misc$link[3],
                                    earg = object@misc$earg[[3]])
  }

  object@post
})


setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "binom2.or"),
  function(object,
           VGAMff,
           ...) {
 if (length(object@post$oratio) == 1 &&
      is.numeric(object@post$oratio)) {
    cat("\nOdds ratio: ", round(object@post$oratio, digits = 4), "\n")
  }
})





















dbinom2.rho <-
  function(mu1,
           mu2 = if (exchangeable) mu1 else stop("'mu2' not specified"),
           rho = 0,
           exchangeable = FALSE,
           colnames = c("00", "01", "10", "11"),
           ErrorCheck = TRUE) {
  if (ErrorCheck) {
    if (!is.Numeric(mu1, positive = TRUE) || max(mu1) >= 1)
      stop("bad input for argument 'mu1'") 
    if (!is.Numeric(mu2, positive = TRUE) || max(mu2) >= 1)
      stop("bad input for argument 'mu2'") 
    if (!is.Numeric(rho) || min(rho) <= -1 || max(rho) >= 1)
      stop("bad input for argument 'rho'") 
    if (exchangeable && max(abs(mu1 - mu2)) > 0.00001)
      stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ")
  }

  nn <- max(length(mu1), length(mu2), length(rho))
  rho <- rep(rho, len = nn)
  mu1 <- rep(mu1, len = nn)
  mu2 <- rep(mu2, len = nn)
  eta1 <- qnorm(mu1)
  eta2 <- qnorm(mu2)
  p11 <- pbinorm(eta1, eta2, cov12 = rho)
  p01 <- mu2 - p11
  p10 <- mu1 - p11
  p00 <- 1.0 - p01 - p10 - p11
  matrix(c(p00, p01, p10, p11), nn, 4,
         dimnames = list(NULL, colnames))
}



rbinom2.rho <-
  function(n, mu1,
           mu2 = if (exchangeable) mu1 else
                   stop("argument 'mu2' not specified"),
           rho = 0,
           exchangeable = FALSE,
           twoCols = TRUE,
           colnames = if (twoCols) c("y1", "y2") else
                      c("00", "01", "10", "11"),
           ErrorCheck = TRUE) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (ErrorCheck) {
    if (!is.Numeric(mu1, positive = TRUE) ||
        max(mu1) >= 1)
      stop("bad input for argument 'mu1'") 
    if (!is.Numeric(mu2, positive = TRUE) ||
        max(mu2) >= 1)
      stop("bad input for argument 'mu2'") 
    if (!is.Numeric(rho) || min(rho) <= -1 ||
        max(rho) >= 1)
      stop("bad input for argument 'rho'") 


    if (exchangeable &&
        max(abs(mu1 - mu2)) > 0.00001)
      stop("argument 'exchangeable' is TRUE but 'mu1' and 'mu2' differ")
  }

  dmat <- dbinom2.rho(mu1 = mu1, mu2 = mu2, rho = rho,
                      exchangeable = exchangeable,
                      ErrorCheck = ErrorCheck)

  answer <- matrix(0, use.n, 2,
                   dimnames = list(NULL,
                                   if (twoCols) colnames else NULL))
  yy <- runif(use.n)
  cs1 <- dmat[, "00"] + dmat[, "01"]
  cs2 <- cs1 + dmat[, "10"]
  index <- (dmat[, "00"] < yy) & (yy <= cs1)
  answer[index, 2] <- 1
  index <- (cs1 < yy) & (yy <= cs2)
  answer[index, 1] <- 1
  index <- (yy > cs2)
  answer[index,] <- 1
  if (twoCols) {
    answer
  } else {
    answer4 <- matrix(0, use.n, 4, dimnames = list(NULL, colnames))
    answer4[cbind(1:use.n, 1 + 2*answer[, 1] + answer[, 2])] <- 1
    answer4
  }
}





binom2.rho.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 binom2.rho <-
  function(lmu = "probit",  # added 20120817, order swapped 20151128
           lrho = "rhobit",
           imu1 = NULL, imu2 = NULL, irho = NULL,
           imethod = 1,
           zero =  "rho",   # 3
           exchangeable = FALSE,
           grho = seq(-0.95, 0.95, by = 0.05),
           nsimEIM = NULL) {



  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  if (lmu != "probit")
    warning("argument 'lmu' should be 'probit'")

    lmu12 <- "probit" # But emu may contain some arguments.
    emu12 <- emu # list()



  if (is.Numeric(nsimEIM)) {
    if (!is.Numeric(nsimEIM, length.arg = 1,
                    integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 100)
      warning("'nsimEIM' should be an integer greater than 100")
  }

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")




  new("vglmff",
  blurb = c("Bivariate probit model\n",
            "Links:    ",
            namesof("mu1", lmu12, earg = emu12), ", ",
            namesof("mu2", lmu12, earg = emu12), ", ",
            namesof("rho", lrho,  earg = erho)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(c(1, 1, 0, 0, 0, 1), 3, 2), x = x,
                           bool = .exchangeable ,
                           constraints = constraints,
                           apply.int = TRUE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .exchangeable = exchangeable, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu1", "mu2", "rho"),
         lmu1 = .lmu12,
         lmu2 = .lmu12,
         lrho = .lrho ,
         zero = .zero )
  }, list( .lmu12 = lmu12, .lrho = lrho,
           .zero = zero ))),


  initialize = eval(substitute(expression({
    mustart.orig <- mustart
    eval(process.binomial2.data.VGAM)

    if (length(mustart.orig))
      mustart <- mustart.orig  # Retain it if inputted

    predictors.names <- c(
        namesof("mu1", .lmu12 , earg = .emu12 , short = TRUE),
        namesof("mu2", .lmu12 , earg = .emu12 , short = TRUE),
        namesof("rho", .lrho ,  earg = .erho,  short = TRUE))

    if (is.null( .nsimEIM )) {
      save.weights <- control$save.weights <- FALSE
    }


    ycounts <- if (is.numeric(extra$orig.w)) y * c(w) / extra$orig.w else
               y * c(w)  # Convert proportions to counts
    if (max(abs(ycounts - round(ycounts))) > 1.0e-6)
       warning("the response (as counts) does not appear to ",
               "be integer-valued. Am rounding to integer values.")
    ycounts <- round(ycounts)  # Make sure it is an integer
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                                          round(w)


    if (is.null(etastart)) {
      if (length(mustart.orig)) {
        mu1.init <- mustart.orig[, 3] + mustart.orig[, 4]
        mu2.init <- mustart.orig[, 2] + mustart.orig[, 4]
      } else if ( .imethod == 1) {
          glm1.fit <- glm(cbind(ycounts[, 3] + ycounts[, 4],
                                ycounts[, 1] + ycounts[, 2]) ~ x - 1,
                          fam = binomial("probit"))
          glm2.fit <- glm(cbind(ycounts[, 2] + ycounts[, 4],
                                ycounts[, 1] + ycounts[, 3]) ~ x - 1,
                          fam = binomial("probit"))
          mu1.init <- fitted(glm1.fit)
          mu2.init <- fitted(glm2.fit)
      } else if ( .imethod == 2) {
          mu1.init <- if (is.Numeric( .imu1 ))
                      rep( .imu1 , length = n) else
                      mu[, 3] + mu[, 4]
          mu2.init <- if (is.Numeric( .imu2 ))
                      rep( .imu2 , length = n) else
                      mu[, 2] + mu[, 4]
      } else {
        stop("bad value for argument 'imethod'")
      }



      binom2.rho.Loglikfun <-
          function(rhoval, y, x, w, extraargs) {
          init.mu1 <-    extraargs$initmu1
          init.mu2 <-    extraargs$initmu2
          ycounts  <-    extraargs$ycounts
          nvec     <-    extraargs$nvec
          eta1 <- qnorm(init.mu1)
          eta2 <- qnorm(init.mu2)
          p11 <- pbinorm(eta1, eta2, cov12 = rhoval)
          p01 <- pmin(init.mu2 - p11, init.mu2)
          p10 <- pmin(init.mu1 - p11, init.mu1)
          p00 <- 1.0 - p01 - p10 - p11
          mumat <- abs(cbind("00" = p00,
                             "01" = p01,
                             "10" = p10,
                             "11" = p11))
          mumat <- mumat / rowSums(mumat)
          mumat[mumat < 1.0e-100] <- 1.0e-100

      sum((if (is.numeric(extraargs$orig.w)) extraargs$orig.w else 1) *
          dmultinomial(x = ycounts, size = nvec, prob = mumat,
                       log = TRUE, dochecking = FALSE))
        }
        rho.grid <- .grho # seq(-0.95, 0.95, len = 31)
        try.this <- grid.search(rho.grid, objfun = binom2.rho.Loglikfun,
                                y = y, x = x, w = w, extraargs = list(
                                orig.w = extra$orig.w,
                                ycounts = ycounts,
                                initmu1 = mu1.init,
                                initmu2 = mu2.init,
                                nvec = nvec ))


      rho.init <- if (is.Numeric( .irho ))
                   rep( .irho , len = n) else {
          try.this
      }

      etastart <- cbind(theta2eta(mu1.init, .lmu12 , earg = .emu12 ),
                        theta2eta(mu2.init, .lmu12 , earg = .emu12 ),
                        theta2eta(rho.init, .lrho ,  earg = .erho ))
      mustart <- NULL # Since etastart has been computed.
    }
  }), list( .lmu12 = lmu12, .lrho = lrho,
            .emu12 = emu12, .erho = erho, 
                            .grho = grho,
                            .irho = irho,
            .imethod = imethod, .nsimEIM = nsimEIM,
            .imu1 = imu1, .imu2 = imu2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pmargin <- cbind(eta2theta(eta[, 1], .lmu12 , earg = .emu12 ),
                     eta2theta(eta[, 2], .lmu12 , earg = .emu12 ))
    rho <- eta2theta(eta[, 3], .lrho , earg = .erho )
    p11 <- pbinorm(eta[, 1], eta[, 2], cov12 = rho)
    p01 <- pmin(pmargin[, 2] - p11, pmargin[, 2])
    p10 <- pmin(pmargin[, 1] - p11, pmargin[, 1])
    p00 <- 1.0 - p01 - p10 - p11
    ansmat <- abs(cbind("00" = p00,
                        "01" = p01,
                        "10" = p10,
                        "11" = p11))
    ansmat / rowSums(ansmat)
  }, list( .lmu12 = lmu12, .lrho = lrho,
           .emu12 = emu12, .erho = erho ))),
  last = eval(substitute(expression({
    misc$link <-    c(mu1 = .lmu12 , mu2 = .lmu12 , rho = .lrho )

    misc$earg <- list(mu1 = .emu12 , mu2 = .emu12 , rho = .erho )

    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lmu12 = lmu12, .lrho = lrho, .nsimEIM = nsimEIM,
            .emu12 = emu12, .erho = erho ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      ycounts <- if (is.numeric(extra$orig.w))
                 y * c(w) / extra$orig.w else
                 y * c(w)  # Convert proportions to counts

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
        warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                round(w)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
        dmultinomial(x = ycounts, size = nvec, prob = mu,
                     log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .erho = erho ))),
  vfamily = c("binom2.rho", "binom2"),






  deriv = eval(substitute(expression({
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
            round(w)
    ycounts <- if (is.numeric(extra$orig.w)) y * c(w) / extra$orig.w else
               y * c(w)  # Convert proportions to counts

    pmargin <- cbind(eta2theta(eta[, 1], .lmu12 , earg = .emu12 ),
                     eta2theta(eta[, 2], .lmu12 , earg = .emu12 ))
    rhovec <- eta2theta(eta[, 3], .lrho , earg = .erho )
    p11 <- pbinorm(eta[, 1], eta[, 2], cov12 = rhovec)
    p01 <- pmargin[, 2] - p11
    p10 <- pmargin[, 1] - p11
    p00 <- 1 - p01 - p10 - p11

    ABmat <- (eta[, 1:2] -
              rhovec * eta[, 2:1]) /  sqrt(pmax(1e5 * .Machine$double.eps,
                                                1.0 - rhovec^2))
    PhiA <- pnorm(ABmat[, 1])
    PhiB <- pnorm(ABmat[, 2])
    onemPhiA <- pnorm(ABmat[, 1], lower.tail = FALSE)
    onemPhiB <- pnorm(ABmat[, 2], lower.tail = FALSE)

    smallno <- 1000 * .Machine$double.eps
    p00[p00 < smallno] <- smallno
    p01[p01 < smallno] <- smallno
    p10[p10 < smallno] <- smallno
    p11[p11 < smallno] <- smallno

    dprob00 <- dbinorm(eta[, 1], eta[, 2], cov12 = rhovec)
    dl.dprob1 <-     PhiB * (ycounts[, 4]/p11 - ycounts[, 2]/p01) +
                 onemPhiB * (ycounts[, 3]/p10 - ycounts[, 1]/p00)
    dl.dprob2 <-     PhiA * (ycounts[, 4]/p11 - ycounts[, 3]/p10) +
                 onemPhiA * (ycounts[, 2]/p01 - ycounts[, 1]/p00)
    dl.drho <- (ycounts[, 4]/p11 - ycounts[, 3]/p10 -
                ycounts[, 2]/p01 + ycounts[, 1]/p00) * dprob00

    dprob1.deta <- dtheta.deta(pmargin[, 1], .lmu12 , earg = .emu12 )
    dprob2.deta <- dtheta.deta(pmargin[, 2], .lmu12 , earg = .emu12 )
    drho.deta <- dtheta.deta(rhovec, .lrho , earg = .erho )
    dthetas.detas <- cbind(dprob1.deta, dprob2.deta, drho.deta)

    (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
    cbind(dl.dprob1,
          dl.dprob2,
          dl.drho) * dthetas.detas
  }), list( .lmu12 = lmu12, .lrho = lrho,
            .emu12 = emu12, .erho = erho ))),
  weight = eval(substitute(expression({
    if (is.null( .nsimEIM )) {
      ned2l.dprob1prob1 <-      PhiB^2 * (1/p11 + 1/p01) +
                            onemPhiB^2 * (1/p10 + 1/p00)
      ned2l.dprob2prob2 <-      PhiA^2 * (1/p11 + 1/p10) +
                            onemPhiA^2 * (1/p01 + 1/p00)
      ned2l.dprob1prob2 <-      PhiA * (    PhiB/p11 - onemPhiB/p10) +
                            onemPhiA * (onemPhiB/p00 -     PhiB/p01)
      ned2l.dprob1rho <-     (PhiB * (1/p11 + 1/p01) -
                          onemPhiB * (1/p10 + 1/p00)) * dprob00
      ned2l.dprob2rho <-     (PhiA * (1/p11 + 1/p10) -
                          onemPhiA * (1/p01 + 1/p00)) * dprob00
      ned2l.drho2 <-  (1/p11 + 1/p01 + 1/p10 + 1/p00) * dprob00^2

      wz <- matrix(0, n, dimm(M))  # 6=dimm(M)
      wz[, iam(1, 1, M)] <- ned2l.dprob1prob1 * dprob1.deta^2
      wz[, iam(2, 2, M)] <- ned2l.dprob2prob2 * dprob2.deta^2
      wz[, iam(3, 3, M)] <- ned2l.drho2 * drho.deta^2
      wz[, iam(1, 2, M)] <- ned2l.dprob1prob2 * dprob1.deta * dprob2.deta
      wz[, iam(2, 3, M)] <- ned2l.dprob2rho * dprob2.deta * drho.deta
      wz[, iam(1, 3, M)] <- ned2l.dprob1rho * dprob1.deta * drho.deta
    } else {
      run.varcov <- 0
      ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
      for (ii in 1:( .nsimEIM )) {
        ysim <- rbinom2.rho(n, mu1 = pmargin[, 1], mu2 = pmargin[, 2],
                            twoCols = FALSE, rho = rhovec)
        dl.dprob1 <-     PhiB * (ysim[, 4]/p11 - ysim[, 2]/p01) +
                     onemPhiB * (ysim[, 3]/p10 - ysim[, 1]/p00)
        dl.dprob2 <-     PhiA * (ysim[, 4]/p11 - ysim[, 3]/p10) +
                     onemPhiA * (ysim[, 2]/p01 - ysim[, 1]/p00)
        dl.drho <- (ysim[, 4]/p11 - ysim[, 3]/p10 -
                    ysim[, 2]/p01 + ysim[, 1]/p00) * dprob00

        rm(ysim)
        temp3 <- cbind(dl.dprob1, dl.dprob2, dl.drho)
          run.varcov <- ((ii-1) * run.varcov +
                 temp3[, ind1$row.index] * temp3[, ind1$col.index]) / ii
      }
      wz <- if (intercept.only)
        matrix(colMeans(run.varcov),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz <- wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
    }
    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}





dnorm2 <- function(x, y, rho = 0, log = FALSE) {


  warning("decommissioning dnorm2() soon; use ",
          "dbinorm(..., cov12 = rho) instead")

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  logdnorm2 <-
    (-0.5*(x * (x - 2*y*rho) + y^2) / (1.0 - rho^2)) - log(2 * pi) -
      0.5 * log1p(-rho^2)

  if (log.arg) {
    logdnorm2
  } else {
    exp(logdnorm2)
  }
}




 pbinorm <-
 function(q1, q2,
          mean1 = 0, mean2 = 0,
          var1 = 1, var2 = 1,
          cov12 = 0) {



  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  rho <- cov12 / (sd1 * sd2)

  if (any(is.na(q1)    | is.na(q2)    |
          is.na(sd1)   | is.na(sd2)   |
          is.na(mean1) | is.na(mean2) | is.na(rho)))
    stop("no NAs allowed in arguments or variables 'q1', 'q2', 'mean1', ",
         "'mean2', 'sd1', 'sd2', 'cov12'")
  if (min(rho) < -1 || max(rho) > +1)
    stop("correlation 'rho' is out of range")


  if (length(mean1) > 1 && length(mean2) == 1 &&
      length(var1) == 1 && length(var2)  == 1 && length(cov12) == 1)
    warning("the call to pnorm2() seems based on the old version ",
            "of the arguments")

  LLL <- max(length(q1), length(q2),
             length(mean1), length(mean2),
             length(sd1), length(sd2),
             length(rho))
  if (length(q1)    != LLL) q1    <- rep(q1,     len = LLL)
  if (length(q2)    != LLL) q2    <- rep(q2,     len = LLL)
  if (length(mean1) != LLL) mean1 <- rep(mean1,  len = LLL)
  if (length(mean2) != LLL) mean2 <- rep(mean2,  len = LLL)
  if (length(sd1)   != LLL) sd1   <- rep(sd1,    len = LLL)
  if (length(sd2)   != LLL) sd2   <- rep(sd2,    len = LLL)
  if (length(rho)   != LLL) rho   <- rep(rho,    len = LLL)

  Zedd1 <- Z1 <- (q1 - mean1) / sd1
  Zedd2 <- Z2 <- (q2 - mean2) / sd2

  is.inf1.neg <- is.infinite(Z1) & Z1 < 0  # -Inf
  is.inf1.pos <- is.infinite(Z1) & Z1 > 0  # +Inf
  is.inf2.neg <- is.infinite(Z2) & Z2 < 0  # -Inf
  is.inf2.pos <- is.infinite(Z2) & Z2 > 0  # +Inf
  Zedd1[is.inf1.neg] <- 0
  Zedd1[is.inf1.pos] <- 0
  Zedd2[is.inf2.neg] <- 0
  Zedd2[is.inf2.pos] <- 0

  ans <- Zedd1
  singler <- ifelse(length(rho) == 1, 1, 0)
  answer <- .C("pnorm2",
       ah = as.double(-Zedd1), ak = as.double(-Zedd2), r = as.double(rho),
       size = as.integer(LLL), singler = as.integer(singler),
       ans = as.double(ans))$ans
  if (any(answer < 0.0))
    warning("some negative values returned")

  answer[is.inf1.neg] <- 0
  answer[is.inf1.pos] <- pnorm(Z2[is.inf1.pos])  # pnorm(Z2[is.inf1.neg])
  answer[is.inf2.neg] <- 0
  answer[is.inf2.pos] <- pnorm(Z1[is.inf2.pos])  # pnorm(Z1[is.inf2.neg])

  answer
}




 pnorm2 <- function(x1, x2,
                    mean1 = 0, mean2 = 0,
                    var1 = 1, var2 = 1,
                    cov12 = 0) {


  warning("decommissioning pnorm2() soon; use ",
          "pbinorm() instead")


  sd1 <- sqrt(var1)
  sd2 <- sqrt(var2)
  rho <- cov12 / (sd1 * sd2)

  if (any(is.na(x1)    | is.na(x2)    |
          is.na(sd1)   | is.na(sd2)   |
          is.na(mean1) | is.na(mean2) | is.na(rho)))
    stop("no NAs allowed in arguments or variables 'x1', 'x2', 'mean1', ",
         "'mean2', 'sd1', 'sd2', 'cov12'")
  if (min(rho) < -1 || max(rho) > +1)
    stop("correlation 'rho' is out of range")


  if (length(mean1) > 1 && length(mean2) == 1 &&
      length(var1) == 1 && length(var2)  == 1 && length(cov12) == 1)
    warning("the call to pnorm2() seems based on the old version ",
            "of the arguments")

  LLL <- max(length(x1), length(x2),
             length(mean1), length(mean2),
             length(sd1), length(sd2),
             length(rho))
  if (length(x1)    != LLL) x1    <- rep(x1,     len = LLL)
  if (length(x2)    != LLL) x2    <- rep(x2,     len = LLL)
  if (length(mean1) != LLL) mean1 <- rep(mean1,  len = LLL)
  if (length(mean2) != LLL) mean2 <- rep(mean2,  len = LLL)
  if (length(sd1)   != LLL) sd1   <- rep(sd1,    len = LLL)
  if (length(sd2)   != LLL) sd2   <- rep(sd2,    len = LLL)
  if (length(rho)   != LLL) rho   <- rep(rho,    len = LLL)

  Z1 <- (x1 - mean1) / sd1
  Z2 <- (x2 - mean2) / sd2

  ans <- Z1
  singler <- ifelse(length(rho) == 1, 1, 0)
  answer <- .C("pnorm2",
       ah = as.double(-Z1), ak = as.double(-Z2), r = as.double(rho),
       size = as.integer(LLL), singler = as.integer(singler),
       ans = as.double(ans))$ans
  if (any(answer < 0.0))
    warning("some negative values returned")
  answer
}










my.dbinom <- function(x,
                      size = stop("no 'size' argument"),
                      prob = stop("no 'prob' argument")) {

  exp(lgamma(size + 1) - lgamma(size - x + 1) - lgamma(x + 1) +
      x * log(prob / (1 - prob)) + size * log1p(-prob))
}



 size.binomial <- function(prob = 0.5, link = "loge") {
  if (any(prob <= 0 | prob >= 1))
    stop("some values of prob out of range")


  link <- as.list(substitute(link))
  earg  <- link2list(link)
  link <- attr(earg, "function.name")



  new("vglmff",
  blurb = c("Binomial with n unknown, prob known (prob = ", prob, ")\n",
            "Links:    ",
            namesof("size", link, tag = TRUE),
            " (treated as real-valued)\n",
            "Variance:  Var(Y) = size * prob * (1-prob);",
            " Var(size) is intractable"),
  initialize = eval(substitute(expression({
    predictors.names <- "size"
    extra$temp2 <- rep( .prob , length = n)

    if (is.null(etastart)) {
      nvec <- (y+0.1)/extra$temp2
      etastart <- theta2eta(nvec, .link )
    }
  }), list( .prob = prob, .link = link ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    nvec <- eta2theta(eta, .link)
    nvec * extra$temp2
  }, list( .link = link ))),
  last = eval(substitute(expression({
    misc$link <- c(size = .link )

    misc$prob <- extra$temp2

  }), list( .link = link ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    nvec <- mu / extra$temp2
    theta2eta(nvec, .link)
  }, list( .link = link ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    nvec <- mu / extra$temp2
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      ll.elts <-
        c(w) * (lgamma(nvec+1) - lgamma(y+1) - lgamma(nvec-y+1) +
                y * log( .prob / (1- .prob )) +
                nvec * log1p(- .prob ))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .prob = prob ))),
  vfamily = c("size.binomial"),
  deriv = eval(substitute(expression({
    nvec <- mu/extra$temp2
    dldnvec <- digamma(nvec+1) - digamma(nvec-y+1) + log1p(-extra$temp2)
    dnvecdeta <- dtheta.deta(nvec, .link)
    c(w) * cbind(dldnvec * dnvecdeta)
  }), list( .link = link ))),
  weight = eval(substitute(expression({
    d2ldnvec2 <- trigamma(nvec+1) - trigamma(nvec-y+1)
    d2ldnvec2[y == 0] <- -sqrt( .Machine$double.eps )
    wz <- -c(w) * dnvecdeta^2 * d2ldnvec2
    wz
  }), list( .link = link ))))
}




 dbetabinom.ab <- function(x, size, shape1, shape2, log = FALSE,
                           Inf.shape = 1e6
                          ) {


  Bigg <- Inf.shape
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  LLL <- max(length(x), length(size), length(shape1), length(shape2))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL)
  if (length(size)   != LLL) size   <- rep(size,   len = LLL)
  if (length(shape1) != LLL) shape1 <- rep(shape1, len = LLL)
  if (length(shape2) != LLL) shape2 <- rep(shape2, len = LLL)

  ans <- x
  ans[TRUE] <- log(0)
  ans[is.na(x)]  <- NA
  ans[is.nan(x)] <- NaN


  ok0 <- !is.na(shape1) & !is.na(shape2) & !is.na(x) & !is.na(size)
  ok <- (round(x) == x) & (x >= 0) & (x <= size) &
        is.finite(shape1) & is.finite(shape2) & ok0
  if (any(ok)) {
    ans[ok] <- lchoose(size[ok], x[ok]) +
               lbeta(shape1[ok]            + x[ok],
                     shape2[ok] + size[ok] - x[ok]) -
               lbeta(shape1[ok], shape2[ok])


    endpt <- (x == size) & ((shape1 < 1/Bigg) | (shape2 < 1/Bigg)) & ok0
    if (any(endpt)) {
      ans[endpt] <- lgamma(size[endpt] + shape1[endpt]) +
                    lgamma(shape1[endpt] + shape2[endpt]) -
                   (lgamma(size[endpt] + shape1[endpt] + shape2[endpt]) +
                    lgamma(shape1[endpt]))
    }




    endpt <- (x == 0) & ((shape1 < 1/Bigg) | (shape2 < 1/Bigg)) & ok0
    if (any(endpt)) {
      ans[endpt] <- lgamma(size[endpt] + shape2[endpt]) +
                    lgamma(shape1[endpt] + shape2[endpt]) -
                   (lgamma(size[endpt] + shape1[endpt] + shape2[endpt]) +
                    lgamma(shape2[endpt]))
    }





    endpt <- ((shape1 > Bigg) | (shape2 > Bigg)) & ok0
    if (any(endpt)) {
      ans[endpt] <- lchoose(size[endpt], x[endpt]) +
                    lgamma(x[endpt] + shape1[endpt]) +
                    lgamma(size[endpt] - x[endpt] + shape2[endpt]) +
                    lgamma(shape1[endpt] + shape2[endpt]) -
                   (lgamma(size[endpt] + shape1[endpt] + shape2[endpt]) +
                    lgamma(shape1[endpt]) +
                    lgamma(shape2[endpt]))
    }
  }  # if (any(ok))



  if (!log.arg) {
    ans <- exp(ans)
  }



  if (FALSE) {
    ok1 <- is.na(shape1)       & is.infinite(shape2)  # rho==0 & prob==0
    ok2 <- is.infinite(shape1) & is.na(shape2)        # rho==0 & prob==1
    ok3 <- is.infinite(shape1) & is.infinite(shape2)  # rho==0 & 0<prob<1
  } else {
    ok1 <-   is.finite(shape1) & is.infinite(shape2)  # rho==0 & prob==0
    ok2 <- is.infinite(shape1) &   is.finite(shape2)  # rho==0 & prob==1
    ok3 <- is.infinite(shape1) & is.infinite(shape2)  # prob undefined

  }

  if (any(ok1))
    ans[ok1] <- dbinom(x = x[ok1], size = size[ok1],
                       prob = shape1[ok1] / (shape1[ok1]+shape2[ok1]),  # 0,
                       log = log.arg)
  if (any(ok2))
    ans[ok2] <- dbinom(x = x[ok2], size = size[ok2],
                       prob = 1,  # Inf / (finite + Inf) == 1
                       log = log.arg)
  if (any(ok3)) {
    ans[ok3] <- dbinom(x = x[ok3], size = size[ok3],
                       prob = shape1[ok3] / (shape1[ok3]+shape2[ok3]),
                       log = log.arg)
  }


  ans[shape1 < 0] <- NaN
  ans[shape2 < 0] <- NaN


  ans
}






 pbetabinom.ab <- function(q, size, shape1, shape2, log.p = FALSE) {

  if (!is.Numeric(q))
    stop("bad input for argument 'q'")
  if (!is.Numeric(size, integer.valued = TRUE))
    stop("bad input for argument 'size'")
  if (!is.Numeric(shape1, positive = TRUE))
    stop("bad input for argument 'shape1'")
  if (!is.Numeric(shape2, positive = TRUE))
    stop("bad input for argument 'shape2'")
  LLL <- max(length(q), length(size), length(shape1), length(shape2))

  if (length(q)       != LLL) q      <- rep(q,      len = LLL)
  if (length(shape1)  != LLL) shape1 <- rep(shape1, len = LLL)
  if (length(shape2)  != LLL) shape2 <- rep(shape2, len = LLL)
  if (length(size)    != LLL) size   <- rep(size,   len = LLL)

  ans <- q   # Retains names(q)
  ans[] <- 0  #  Set all elements to 0

  if (max(abs(size   -   size[1])) < 1.0e-08 &&
      max(abs(shape1 - shape1[1])) < 1.0e-08 &&
      max(abs(shape2 - shape2[1])) < 1.0e-08) {
    if (any(is.infinite(qstar <- floor(q))))
      stop("argument 'q' must be finite")
    temp <- if (max(qstar) >= 0) {
              dbetabinom.ab(0:max(qstar), size = size[1],
                            shape1 = shape1[1],
                            shape2 = shape2[1])
            } else {
              0 * qstar
            }
      unq <- unique(qstar)
    for (ii in unq) {
      index <- (qstar == ii)
      ans[index] <- if (ii >= 0) sum(temp[1:(1+ii)]) else 0
    }
  } else {
    for (ii in 1:LLL) {
      qstar <- floor(q[ii])
      ans[ii] <- if (qstar >= 0) {
                   sum(dbetabinom.ab(x = 0:qstar, size = size[ii],
                                     shape1 = shape1[ii],
                                     shape2 = shape2[ii]))
                 } else 0
    }
  }
  if (log.p) log(ans) else ans
}



 rbetabinom.ab <- function(n, size, shape1, shape2,
                           .dontuse.prob = NULL) {
 #                         checkargs = TRUE

  if (!is.Numeric(size, integer.valued = TRUE))
    stop("bad input for argument 'size'")
  if (any(shape1 < 0, na.rm = TRUE))
    stop("negative values for argument 'shape1' not allowed")
  if (any(shape2 < 0, na.rm = TRUE))
    stop("negative values for argument 'shape2' not allowed")

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n
  if (length(size)   != use.n) size   <- rep(size,   len = use.n)
  if (length(shape1) != use.n) shape1 <- rep(shape1, len = use.n)
  if (length(shape2) != use.n) shape2 <- rep(shape2, len = use.n)

  ans <- rep(NA_real_, len = use.n)
  okay0 <- is.finite(shape1) & is.finite(shape2)
  if (smalln <- sum(okay0))
    ans[okay0] <- rbinom(n = smalln, size = size[okay0],
                         prob = rbeta(n = smalln, shape1 = shape1[okay0],
                                                  shape2 = shape2[okay0]))

  okay1 <- is.na(shape1)       & is.infinite(shape2)  # rho = 0 and prob == 0
  okay2 <- is.infinite(shape1) & is.na(shape2)       # rho = 0 and prob == 1
  okay3 <- is.infinite(shape1) & is.infinite(shape2)  # rho = 0 and 0 < prob < 1

  if (sum.okay1 <- sum(okay1))
    ans[okay1] <- rbinom(n = sum.okay1, size = size[okay1],
                         prob = 0)
  if (sum.okay2 <- sum(okay2))
    ans[okay2] <- rbinom(n = sum.okay2, size = size[okay2],
                         prob = 1)
  if (sum.okay3 <- sum(okay3)) {
    if (length( .dontuse.prob ) != use.n)
      .dontuse.prob   <- rep(.dontuse.prob,   len = use.n)
    ans[okay3] <- rbinom(n = sum.okay3, size = size[okay3],
                         prob = .dontuse.prob[okay3])
  }

  ans
}







 dbetabinom <- function(x, size, prob, rho = 0, log = FALSE) {
  dbetabinom.ab(x = x, size = size, shape1 = prob*(1-rho)/rho,
                shape2 = (1-prob)*(1-rho)/rho, log = log)
}


 pbetabinom <- function(q, size, prob, rho, log.p = FALSE) {
  pbetabinom.ab(q = q, size = size, shape1 = prob*(1-rho)/rho,
                shape2 = (1-prob)*(1-rho)/rho, log.p = log.p)
}


 rbetabinom <- function(n, size, prob, rho = 0) {
  rbetabinom.ab(n = n, size = size, shape1 = prob*(1-rho)/rho,
                shape2 = (1-prob)*(1-rho)/rho,
                .dontuse.prob = prob)
}



 expected.betabin.ab <- function(nvec, shape1, shape2, first) {



  NN <- length(nvec)
  ans <- rep(0.0, len = NN)
  if (first) {
    for (ii in 1:NN) {
      temp639 <- lbeta(shape1[ii], shape2[ii])
      yy <- 0:nvec[ii]
      ans[ii] <- ans[ii] + sum(trigamma(shape1[ii] + yy) *
                exp(lchoose(nvec[ii], yy) +
                    lbeta(shape1[ii]+yy, shape2[ii]+nvec[ii]-yy) -
                    temp639))
    }
  } else {
    for (ii in 1:NN) {
      temp639 <- lbeta(shape1[ii], shape2[ii])
      yy <- 0:nvec[ii]
      ans[ii] <- ans[ii] + sum(trigamma(nvec[ii]+shape2[ii] - yy) *
                exp(lchoose(nvec[ii], yy) +
                    lbeta(shape1[ii]+yy, shape2[ii]+nvec[ii]-yy) -
                    temp639))
    }
  }
  ans
}



betabinomialff.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}




 betabinomialff <-
  function(lshape1 = "loge", lshape2 = "loge",
           ishape1 = 1, ishape2 = NULL, imethod = 1,
           ishrinkage = 0.95, nsimEIM = NULL,
           zero = NULL) {




  lshape1 <- as.list(substitute(lshape1))
  earg1 <- link2list(lshape1)
  lshape1 <- attr(earg1, "function.name")


  lshape2 <- as.list(substitute(lshape2))
  earg2 <- link2list(lshape2)
  lshape2 <- attr(earg2, "function.name")


  if (!is.Numeric(ishape1, positive = TRUE))
    stop("bad input for argument 'ishape1'")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1, 2 or 3")

  if (length(ishape2) && !is.Numeric(ishape2, positive = TRUE))
    stop("bad input for argument 'ishape2'")

  if (!is.null(nsimEIM)) {
    if (!is.Numeric(nsimEIM, length.arg = 1,
                    integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 10)
      warning("'nsimEIM' should be an integer greater than 10, say")
  }


  new("vglmff",
  blurb = c("Beta-binomial model\n",
            "Links:    ",
            namesof("shape1", lshape1, earg = earg1), ", ",
            namesof("shape2", lshape2, earg = earg2), "\n",
            "Mean:     mu = shape1 / (shape1+shape2)", "\n",
            "Variance: mu * (1-mu) * (1+(w-1)*rho) / w, ",
                       "where rho = 1 / (shape1+shape2+1)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("shape1", "shape2"),
         lshape1 = .lshape1 ,
         lshape2 = .lshape2 ,
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w

    if (is.null( .nsimEIM )) {
      save.weights <- control$save.weights <- FALSE
    }

    mustart.orig <- mustart
    eval(binomialff()@initialize)   # Note: n,w,y,mustart is changed 
    if (length(mustart.orig))
      mustart <- mustart.orig  # Retain it if inputted
    predictors.names <-
         c(namesof("shape1", .lshape1 , earg = .earg1 , tag = FALSE),
           namesof("shape2", .lshape2 , earg = .earg2 , tag = FALSE))

    if (!length(etastart)) {

      mustart.use <- if (length(mustart.orig)) mustart.orig else
                     mustart

      shape1 <- rep( .ishape1 , len = n)
      shape2 <- if (length( .ishape2 )) {
                  rep( .ishape2 , len = n)
                } else if (length(mustart.orig)) {
                  shape1 * (1 / mustart.use - 1)
                } else if ( .imethod == 1) {
                  shape1 * (1 / weighted.mean(y, w)  - 1)
                } else if ( .imethod == 2) {
                  temp777 <- .ishrinkage * weighted.mean(y, w) +
                            (1 - .ishrinkage ) * y
                  shape1 * (1 / temp777 - 1)
                } else {
                  shape1 * (1 / weighted.mean(mustart.use, w) - 1)
                }
      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                 y * w # Convert proportions to counts
      if (max(abs(ycounts - round(ycounts))) > 1.0e-6)
         warning("the response (as counts) does not appear to ",
                 "be integer-valued. Am rounding to integer values.")
      ycounts <- round(ycounts)  # Make sure it is an integer
      etastart <- cbind(theta2eta(shape1, .lshape1 , earg = .earg1 ),
                        theta2eta(shape2, .lshape2 , earg = .earg2 ))
      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .earg1 = earg1, .earg2 = earg2,
            .ishape1 = ishape1, .ishape2 = ishape2,
            .nsimEIM = nsimEIM,
            .imethod = imethod, .ishrinkage = ishrinkage ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .earg1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .earg2 )
    shape1 / (shape1 + shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .earg1 = earg1, .earg2 = earg2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("shape1" = .lshape1 , "shape2" = .lshape2 )

    misc$earg <- list("shape1" = .earg1 ,   "shape2" = .earg2   )

    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .earg1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .earg2 )

    misc$rho <- 1 / (shape1 + shape2 + 1)
    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$zero <- .zero
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .earg1 = earg1, .earg2 = earg2,
            .nsimEIM = nsimEIM, .zero = zero ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
               y * w # Convert proportions to counts

    smallno <- 1.0e4 * .Machine$double.eps
    if (max(abs(ycounts - round(ycounts))) > smallno)
      warning("converting 'ycounts' to integer in @loglikelihood")
    ycounts <- round(ycounts)

    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .earg1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .earg2 )
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dbetabinom.ab(x = ycounts, size = nvec, shape1 = shape1,
                       shape2 = shape2, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape1 = lshape1, .lshape2 = lshape2,
           .earg1 = earg1, .earg2 = earg2 ))),
  vfamily = c("betabinomialff"),





  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    w <- pwts
    eta <- predict(object)
    extra <- object@extra
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .earg1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .earg2 )
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
    rbetabinom.ab(nsim * length(shape1), size = nvec,
                  shape1 = shape1,
                  shape2 = shape2)
  }, list( .lshape1 = lshape1, .lshape2 = lshape2, 
           .earg1 = earg1, .earg2 = earg2 ))),





  deriv = eval(substitute(expression({
    nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
    ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
              y * w # Convert proportions to counts
    shape1 <- eta2theta(eta[, 1], .lshape1 , earg = .earg1 )
    shape2 <- eta2theta(eta[, 2], .lshape2 , earg = .earg2 )

    dshape1.deta <- dtheta.deta(shape1, .lshape1 , earg = .earg1 )
    dshape2.deta <- dtheta.deta(shape2, .lshape2 , earg = .earg2 )

    dl.dshape1 <- digamma(shape1+ycounts) -
                  digamma(shape1+shape2+nvec) -
                  digamma(shape1) + digamma(shape1 + shape2)
    dl.dshape2 <- digamma(nvec + shape2 - ycounts) -
                  digamma(shape1 + shape2 + nvec) -
                  digamma(shape2) + digamma(shape1 + shape2)

    (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
    cbind(dl.dshape1 * dshape1.deta,
          dl.dshape2 * dshape2.deta)
  }), list( .lshape1 = lshape1, .lshape2 = lshape2,
            .earg1 = earg1, .earg2 = earg2 ))),
  weight = eval(substitute(expression({
    if (is.null( .nsimEIM)) {
      wz <- matrix(NA_real_, n, dimm(M))  #3=dimm(2)
      wz[, iam(1, 1, M)] <- -(expected.betabin.ab(nvec,shape1,shape2,
                                              TRUE) -
                          trigamma(shape1+shape2+nvec) -
                          trigamma(shape1) + trigamma(shape1+shape2)) *
                          dshape1.deta^2
      wz[, iam(2, 2, M)] <- -(expected.betabin.ab(nvec,shape1,shape2,
                                              FALSE) -
                          trigamma(shape1+shape2+nvec) -
                          trigamma(shape2) + trigamma(shape1+shape2)) *
                          dshape2.deta^2
      wz[, iam(2, 1, M)] <- -(trigamma(shape1+shape2) -
                          trigamma(shape1+shape2+nvec)) *
                          dshape1.deta * dshape2.deta
      wz * (if (is.numeric(extra$orig.w)) extra$orig.w else 1)
    } else {
      run.varcov <- 0
      ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
      dthetas.detas <- cbind(dshape1.deta, dshape2.deta)

      for (ii in 1:( .nsimEIM )) {
        ysim <- rbetabinom.ab(n = n, size = nvec, shape1 = shape1,
                             shape2 = shape2)
                             checkargs = .checkargs
        dl.dshape1 <- digamma(shape1+ysim) -
                      digamma(shape1+shape2+nvec) -
                      digamma(shape1) + digamma(shape1+shape2)
        dl.dshape2 <- digamma(nvec+shape2-ysim) -
                      digamma(shape1+shape2+nvec) -
                      digamma(shape2) + digamma(shape1+shape2)
        rm(ysim)
        temp3 <- cbind(dl.dshape1, dl.dshape2)  # n x M matrix
        run.varcov <- ((ii-1) * run.varcov +
                     temp3[, ind1$row.index]*
                     temp3[, ind1$col.index]) / ii
      }
      wz <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz <- wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
      wz * (if (is.numeric(extra$orig.w)) extra$orig.w else 1)
    }
  }), list( .lshape1 = lshape1, .lshape2 = lshape2, 
            .earg1 = earg1, .earg2 = earg2,
            .nsimEIM = nsimEIM ))))
}



 betageometric <- function(lprob = "logit", lshape = "loge",
                           iprob = NULL, ishape = 0.1,
                           moreSummation = c(2, 100),
                           tolerance = 1.0e-10,
                           zero = NULL) {
  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (!is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (!is.Numeric(moreSummation, positive = TRUE,
                  length.arg = 2, integer.valued = TRUE))
    stop("bad input for argument 'moreSummation'")
  if (!is.Numeric(tolerance, positive = TRUE, length.arg = 1) ||
      1.0 - tolerance >= 1.0)
      stop("bad input for argument 'tolerance'")



  new("vglmff",
  blurb = c("Beta-geometric distribution\n",
            "Links:    ",
            namesof("prob",  lprob,  earg = eprob), ", ",
            namesof("shape", lshape, earg = eshape)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("prob", "shape"),
         lprob  = .lprob ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .lprob = lprob, .lshape = lshape,
           .zero = zero ))),


  initialize = eval(substitute(expression({
    eval(geometric()@initialize)

    predictors.names <-
         c(namesof("prob",  .lprob  , earg = .eprob  , tag = FALSE),
           namesof("shape", .lshape , earg = .eshape , tag = FALSE))

    if (length( .iprob ))
      prob.init <- rep( .iprob , len = n)

    if (!length(etastart) ||
      ncol(cbind(etastart)) != 2) {
      shape.init <- rep( .ishape , len = n)
      etastart <-
        cbind(theta2eta(prob.init,  .lprob ,  earg = .eprob ),
              theta2eta(shape.init, .lshape , earg = .eshape ))
      }
  }), list( .iprob = iprob, .ishape = ishape, .lprob = lprob,
            .eprob = eprob, .eshape = eshape,
            .lshape = lshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob  <- eta2theta(eta[, 1], .lprob ,  earg = .eprob )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    mymu <- (1-prob) / (prob - shape)
    ifelse(mymu >= 0, mymu, NA)
  }, list( .lprob = lprob, .lshape = lshape,
           .eprob = eprob, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$link <-    c("prob" = .lprob , "shape" = .lshape )

    misc$earg <- list("prob" = .eprob , "shape" = .eshape )

    if (intercept.only) {
      misc$shape1 <- shape1[1]  # These quantities computed in @deriv
      misc$shape2 <- shape2[1]
    }
    misc$expected <- TRUE
    misc$tolerance <- .tolerance
    misc$zero <- .zero
    misc$moreSummation = .moreSummation
  }), list( .lprob = lprob, .lshape = lshape,
            .eprob = eprob, .eshape = eshape,
            .tolerance = tolerance,
            .moreSummation = moreSummation, .zero = zero ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    prob  <- eta2theta(eta[, 1], .lprob  , earg = .eprob  )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    ans <- log(prob)
    maxy <- max(y)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
    for (ii in 1:maxy) {
      index <- (ii <= y)
      ans[index] <- ans[index] +
                    log1p(-prob[index] + (ii-1) * shape[index]) -
                    log1p((ii-1) * shape[index])
    }
    ans <- ans - log1p((y+1-1) * shape)




      ll.elts <- w * ans
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lshape = lshape,
           .eprob = eprob, .eshape = eshape ))),
  vfamily = c("betageometric"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    prob  <- eta2theta(eta[, 1], .lprob  , earg = .eprob  )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    rbetageom(nsim * length(shape),
              shape1 = shape, shape2 = shape)
  }, list( .lprob = lprob, .lshape = lshape,
           .eprob = eprob, .eshape = eshape ))),




  deriv = eval(substitute(expression({
    prob  <- eta2theta(eta[, 1], .lprob ,  earg = .eprob )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    shape1 <-      prob  / shape
    shape2 <- (1 - prob) / shape
    dprob.deta  <- dtheta.deta(prob , .lprob  , earg = .eprob  )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dl.dprob <- 1 / prob
    dl.dshape <- 0 * y
    maxy <- max(y)
    for (ii in 1:maxy) {
      index <- (ii <= y)
      dl.dprob[index] <- dl.dprob[index] -
                         1/(1-prob[index]+(ii-1) * shape[index])
      dl.dshape[index] <- dl.dshape[index] +
                         (ii-1)/(1-prob[index]+(ii-1) * shape[index]) -
                         (ii-1)/(1+(ii-1) * shape[index])
    }
    dl.dshape <- dl.dshape - (y+1 -1)/(1+(y+1 -1) * shape)
    c(w) * cbind(dl.dprob * dprob.deta,
                 dl.dshape * dshape.deta)
  }), list( .lprob = lprob, .lshape = lshape,
            .eprob = eprob, .eshape = eshape ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, dimm(M))  #3=dimm(2)
    wz[, iam(1, 1, M)] <- 1 / prob^2
    moresum <- .moreSummation
    maxsummation <- round(maxy * moresum[1] + moresum[2])
    for (ii in 3:maxsummation) {
      temp7 <- 1 - pbetageom(q = ii-1-1, shape1 = shape1,
                                         shape2 = shape2)
      denom1 <- (1-prob+(ii-2)*shape)^2
      denom2 <- (1+(ii-2)*shape)^2
      wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] + temp7 / denom1
      wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)] - (ii-2) * temp7 / denom1
      wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] + (ii-2)^2 * temp7 / denom1 -
                        (ii-1)^2 * temp7 / denom2
      if (max(temp7) < .tolerance ) break
    }
    ii <- 2
    temp7 <- 1 - pbetageom(q=ii-1-1, shape1 = shape1, shape2 = shape2)
    denom1 <- (1-prob+(ii-2)*shape)^2
    denom2 <- (1+(ii-2)*shape)^2

    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] + temp7 / denom1
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] - (ii-1)^2 * temp7 / denom2
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] * dprob.deta^2
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] * dshape.deta^2
    wz[, iam(2, 1, M)] <- wz[, iam(2, 1, M)] * dprob.deta * dshape.deta
    c(w) * wz
  }), list( .lprob = lprob, .lshape = lshape,
            .eprob = eprob, .eshape = eshape,
            .moreSummation = moreSummation,
            .tolerance = tolerance ))))
}




 seq2binomial <- function(lprob1 = "logit", lprob2 = "logit",
                          iprob1 = NULL,    iprob2 = NULL,
                          parallel = FALSE,  # apply.parint = TRUE,
                          zero = NULL) {
  apply.parint <- TRUE

  lprob1 <- as.list(substitute(lprob1))
  eprob1 <- link2list(lprob1)
  lprob1 <- attr(eprob1, "function.name")

  lprob2 <- as.list(substitute(lprob2))
  eprob2 <- link2list(lprob2)
  lprob2 <- attr(eprob2, "function.name")


  if (length(iprob1) &&
     (!is.Numeric(iprob1, positive = TRUE) ||
     max(iprob1) >= 1))
    stop("bad input for argument 'iprob1'")
  if (length(iprob2) &&
     (!is.Numeric(iprob2, positive = TRUE) ||
     max(iprob2) >= 1))
    stop("bad input for argument 'iprob2'")


  new("vglmff",
  blurb = c("Sequential binomial distribution ",
            "(Crowder and Sweeting, 1989)\n",
            "Links:    ",
            namesof("prob1", lprob1, earg = eprob1), ", ",
            namesof("prob2", lprob2, earg = eprob2)),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints,
                           apply.int = .apply.parint )
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .parallel = parallel,
            .apply.parint = apply.parint,
            .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("prob1", "prob2"),
         lprob1 = .lprob1 ,
         lprob2 = .lprob2 ,
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({
    if (!is.vector(w))
      stop("the 'weights' argument must be a vector")

    if (any(abs(w - round(w)) > 0.000001))
      stop("the 'weights' argument does not seem to be integer-valued")


    if (ncol(y <- cbind(y)) != 2)
      stop("the response must be a 2-column matrix")

    if (any(y < 0 | y > 1))
      stop("the response must have values between 0 and 1")


    w <- round(w)
    rvector <- w * y[, 1]
    if (any(abs(rvector - round(rvector)) > 1.0e-8))
      warning("number of successes in column one ",
              "should be integer-valued")
    svector <- rvector * y[, 2]
    if (any(abs(svector - round(svector)) > 1.0e-8))
      warning("number of successes in",
              " column two should be integer-valued")

    predictors.names <-
        c(namesof("prob1", .lprob1,earg =  .eprob1, tag = FALSE),
          namesof("prob2", .lprob2,earg =  .eprob2, tag = FALSE))

    prob1.init <- if (is.Numeric( .iprob1))
                   rep( .iprob1 , len = n) else
                   rep(weighted.mean(y[, 1], w = w), len = n)
    prob2.init <- if (is.Numeric( .iprob2 ))
                   rep( .iprob2 , length = n) else
                   rep(weighted.mean(y[, 2], w = w*y[, 1]),
                       length = n)

    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(prob1.init, .lprob1, earg = .eprob1),
              theta2eta(prob2.init, .lprob2, earg = .eprob2))
    }
  }), list( .iprob1 = iprob1, .iprob2 = iprob2,
            .lprob1 = lprob1, .lprob2 = lprob2,
            .eprob1 = eprob1, .eprob2 = eprob2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob1 <- eta2theta(eta[, 1], .lprob1, earg = .eprob1)
    prob2 <- eta2theta(eta[, 2], .lprob2, earg = .eprob2)
    cbind(prob1, prob2)
  }, list( .lprob1 = lprob1, .lprob2 = lprob2,
           .eprob1 = eprob1, .eprob2 = eprob2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("prob1" = .lprob1 , "prob2" = .lprob2 )

    misc$earg <- list("prob1" = .eprob1 , "prob2" = .eprob2 )

    misc$expected <- TRUE
    misc$zero <- .zero
    misc$parallel <- .parallel
    misc$apply.parint <- .apply.parint
  }), list( .lprob1 = lprob1, .lprob2 = lprob2,
            .eprob1 = eprob1, .eprob2 = eprob2,
            .parallel = parallel,
            .apply.parint = apply.parint,
            .zero = zero ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    prob1 <- eta2theta(eta[, 1], .lprob1, earg = .eprob1)
    prob2 <- eta2theta(eta[, 2], .lprob2, earg = .eprob2)

    smallno <- 100 * .Machine$double.eps
    prob1 <- pmax(prob1, smallno)
    prob1 <- pmin(prob1, 1-smallno)
    prob2 <- pmax(prob2, smallno)
    prob2 <- pmin(prob2, 1-smallno)
    mvector <- w
    rvector <- w * y[, 1]
    svector <- rvector * y[, 2]
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ans1 <-
        dbinom(rvector, size = mvector, prob = prob1, log = TRUE) +
        dbinom(svector, size = rvector, prob = prob2, log = TRUE)

      ll.elts <- ans1
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob1 = lprob1, .lprob2 = lprob2,
           .eprob1 = eprob1, .eprob2 = eprob2 ))),
  vfamily = c("seq2binomial"),
  deriv = eval(substitute(expression({
    prob1 <- eta2theta(eta[, 1], .lprob1, earg = .eprob1)
    prob2 <- eta2theta(eta[, 2], .lprob2, earg = .eprob2)
    smallno <- 100 * .Machine$double.eps
    prob1 <- pmax(prob1, smallno)
    prob1 <- pmin(prob1, 1-smallno)
    prob2 <- pmax(prob2, smallno)
    prob2 <- pmin(prob2, 1-smallno)

    dprob1.deta <- dtheta.deta(prob1, .lprob1, earg = .eprob1)
    dprob2.deta <- dtheta.deta(prob2, .lprob2, earg = .eprob2)

    mvector <- w
    rvector <- w * y[, 1]
    svector <- rvector * y[, 2]

    dl.dprob1 <- rvector / prob1 - (mvector-rvector) / (1-prob1)
    dl.dprob2 <- svector / prob2 - (rvector-svector) / (1-prob2)

    cbind(dl.dprob1 * dprob1.deta, dl.dprob2 * dprob2.deta)
  }), list( .lprob1 = lprob1, .lprob2 = lprob2,
            .eprob1 = eprob1, .eprob2 = eprob2 ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, M)
    wz[, iam(1, 1, M)] <- (dprob1.deta^2) / (prob1 * (1-prob1))
    wz[, iam(2, 2, M)] <- (dprob2.deta^2) * prob1 / (prob2 * (1-prob2))
    c(w) * wz
  }), list( .lprob1 = lprob1, .lprob2 = lprob2,
            .eprob1 = eprob1, .eprob2 = eprob2 ))))
}



 zipebcom   <- function(lmu12 = "cloglog",
                        lphi12 = "logit",
                        loratio = "loge",
                        imu12 = NULL, iphi12 = NULL,
                        ioratio = NULL,
                        zero = c("phi12", "oratio"),
                        tol = 0.001, addRidge = 0.001) {


  lmu12 <- as.list(substitute(lmu12))
  emu12 <- link2list(lmu12)
  lmu12 <- attr(emu12, "function.name")

  lphi12 <- as.list(substitute(lphi12))
  ephi12 <- link2list(lphi12)
  lphi12 <- attr(ephi12, "function.name")

  loratio <- as.list(substitute(loratio))
  eoratio <- link2list(loratio)
  loratio <- attr(eoratio, "function.name")



  if (!is.Numeric(tol, positive = TRUE, length.arg = 1) ||
      tol > 0.1)
      stop("bad input for argument 'tol'") 
  if (!is.Numeric(addRidge, length.arg = 1, positive = TRUE) ||
      addRidge > 0.5)
    stop("bad input for argument 'addRidge'") 

  if (lmu12 != "cloglog")
    warning("argument 'lmu12' should be 'cloglog'")


  new("vglmff",
  blurb = c("Exchangeable bivariate ", lmu12,
            " odds-ratio model based on\n",
            "a zero-inflated Poisson distribution\n\n",
            "Links:    ",
            namesof("mu12",   lmu12,   earg = emu12), ", ",
            namesof("phi12",  lphi12,  earg = ephi12), ", ",
            namesof("oratio", loratio, earg = eoratio)),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
    }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("mu12", "phi12", "oratio"),
         lmu12   = .lmu12 ,
         lphi12  = .lphi12 ,
         loratio = .loratio ,
         zero = .zero )
  }, list( .zero = zero,
           .lmu12 = lmu12, .lphi12 = lphi12, .loratio = loratio
         ))),



  initialize = eval(substitute(expression({
    eval(process.binomial2.data.VGAM)

    predictors.names <- c(
             namesof("mu12",   .lmu12   , earg = .emu12   , tag = FALSE),
             namesof("phi12",  .lphi12  , earg = .ephi12  , tag = FALSE),
             namesof("oratio", .loratio , earg = .eoratio , tag = FALSE))

    propY1.eq.0 <- weighted.mean(y[,'00'], w) + weighted.mean(y[,'01'], w)
    propY2.eq.0 <- weighted.mean(y[,'00'], w) + weighted.mean(y[,'10'], w)
    if (length( .iphi12) && any( .iphi12 > propY1.eq.0))
      warning("iphi12 must be less than the sample proportion of Y1==0")
    if (length( .iphi12) && any( .iphi12 > propY2.eq.0))
      warning("iphi12 must be less than the sample proportion of Y2==0")

    if (!length(etastart)) {
        pstar.init <- ((mu[, 3]+mu[, 4]) + (mu[, 2]+mu[, 4])) / 2
        phi.init <- if (length(.iphi12)) rep(.iphi12, len = n) else
            min(propY1.eq.0 * 0.95, propY2.eq.0 * 0.95, pstar.init/1.5)
        oratio.init <- if (length( .ioratio)) rep( .ioratio, len = n) else
                  mu[, 4]*mu[, 1]/(mu[, 2]*mu[, 3])
        mu12.init <- if (length(.imu12)) rep(.imu12, len = n) else
            pstar.init / (1-phi.init)

        etastart <- cbind(
            theta2eta(mu12.init,   .lmu12 ,   earg = .emu12 ),
            theta2eta(phi.init,    .lphi12,  earg = .ephi12),
            theta2eta(oratio.init, .loratio, earg = .eoratio))
      }
  }), list( .lmu12 = lmu12, .lphi12 = lphi12, .loratio = loratio,
            .emu12 = emu12, .ephi12 = ephi12, .eoratio = eoratio,
            .imu12 = imu12, .iphi12 = iphi12, .ioratio = ioratio ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    A1vec  <- eta2theta(eta[, 1], .lmu12 ,  earg = .emu12 )
    phivec <- eta2theta(eta[, 2], .lphi12, earg = .ephi12)
    pmargin <- matrix((1 - phivec) * A1vec, nrow(eta), 2)
    oratio <- eta2theta(eta[, 3], .loratio, earg = .eoratio)
    a.temp <- 1 + (pmargin[, 1]+pmargin[, 2])*(oratio-1)
    b.temp <- -4 * oratio * (oratio-1) * pmargin[, 1] * pmargin[, 2]
    temp <- sqrt(a.temp^2 + b.temp)
    pj4 <- ifelse(abs(oratio-1) < .tol, pmargin[, 1]*pmargin[, 2],
                 (a.temp-temp)/(2*(oratio-1)))
    pj2 <- pmargin[, 2] - pj4
    pj3 <- pmargin[, 1] - pj4
    cbind("00" = 1-pj4-pj2-pj3, "01" = pj2, "10" = pj3, "11" = pj4)
  }, list( .tol = tol,
           .lmu12 = lmu12, .lphi12 = lphi12, .loratio = loratio,
           .emu12 = emu12, .ephi12 = ephi12, .eoratio = eoratio ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu12"= .lmu12 ,
                      "phi12" = .lphi12,
                      "oratio" = .loratio)

    misc$earg <- list("mu12"= .emu12 ,
                      "phi12" = .ephi12,
                      "oratio" = .eoratio)

    misc$tol <- .tol
    misc$expected <- TRUE
    misc$addRidge <- .addRidge
  }), list( .tol = tol, .addRidge = addRidge,
            .lmu12 = lmu12, .lphi12 = lphi12, .loratio = loratio,
            .emu12 = emu12, .ephi12 = ephi12, .eoratio = eoratio ))),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

        ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                  y * w # Convert proportions to counts
        nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                  round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
          warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                      log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("zipebcom"),
  deriv = eval(substitute(expression({
    A1vec  <- eta2theta(eta[, 1], .lmu12 ,  earg = .emu12 )
    smallno <- .Machine$double.eps^(2/4)
    A1vec[A1vec > 1.0 -smallno] <- 1.0 - smallno

    phivec <- eta2theta(eta[, 2], .lphi12, earg = .ephi12)
    pmargin <- matrix((1 - phivec) * A1vec, n, 2)
    oratio <- eta2theta(eta[, 3], .loratio, earg = .eoratio)

    Vab <- 1 / (1/mu[, 1] + 1/mu[, 2] + 1/mu[, 3] + 1/mu[, 4])
    Vabc <- 1/mu[, 1] + 1/mu[, 2]
    denom3 <- 2 * oratio * mu[, 2] + mu[, 1] + mu[, 4]
    temp1 <- oratio * mu[, 2] + mu[, 4]
    dp11star.dp1unstar <- 2*(1-phivec)*Vab * Vabc
    dp11star.dphi1 <- -2 * A1vec * Vab * Vabc
    dp11star.doratio <- Vab / oratio 
    yandmu <- (y[, 1]/mu[, 1] - y[, 2]/mu[, 2] - y[, 3]/mu[, 3] +
               y[, 4]/mu[, 4])
    dp11.doratio <- Vab / oratio
    check.dl.doratio <- yandmu * dp11.doratio

    cyandmu <- (y[, 2]+y[, 3])/mu[, 2] - 2 * y[, 1]/mu[, 1]
    dl.dmu1 <- dp11star.dp1unstar * yandmu + (1-phivec) * cyandmu
    dl.dphi1 <- dp11star.dphi1 * yandmu - A1vec * cyandmu
    dl.doratio <- check.dl.doratio
    dthetas.detas =
      cbind(dtheta.deta(A1vec,  .lmu12 ,   earg = .emu12 ),
            dtheta.deta(phivec, .lphi12,  earg = .ephi12),
            dtheta.deta(oratio, .loratio, earg = .eoratio))
    c(w) * cbind(dl.dmu1,
                 dl.dphi1,
                 dl.doratio) * dthetas.detas
  }), list( .lmu12 = lmu12, .lphi12 = lphi12, .loratio = loratio,
            .emu12 = emu12, .ephi12 = ephi12, .eoratio = eoratio ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, 4)
    alternwz11 <- 2 * (1-phivec)^2 *
                  (2/mu[, 1] + 1/mu[, 2] - 2*Vab*Vabc^2) *
                  (dthetas.detas[, 1])^2
    wz[, iam(1, 1, M)] <- alternwz11

    alternwz22 <- 2* A1vec^2 *(2/mu[, 1] + 1/mu[, 2] - 2*Vab*Vabc^2) *
                  (dthetas.detas[, 2])^2
    wz[, iam(2, 2, M)] <- alternwz22

    alternwz12 <- -2*A1vec*(1-phivec)*
                  (2/mu[, 1] + 1/mu[, 2] - 2*Vab*Vabc^2) *
                   dthetas.detas[, 1] * dthetas.detas[, 2]
    wz[, iam(1, 2, M)] <- alternwz12

    alternwz33 <- (Vab / oratio^2) * dthetas.detas[, 3]^2
    wz[, iam(3, 3, M)] <- alternwz33

    wz[, 1:2] <- wz[, 1:2] * (1 + .addRidge)
    c(w) * wz
  }), list( .addRidge = addRidge ))))
}






 binom2.Rho <- function(rho = 0, imu1 = NULL, imu2 = NULL, 
                        exchangeable = FALSE, nsimEIM = NULL) {
  lmu12 <- "probit"
  emu12 <- list()

  if (is.Numeric(nsimEIM)) {
    if (!is.Numeric(nsimEIM, length.arg = 1,
                    integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
    if (nsimEIM <= 100)
      warning("'nsimEIM' should be an integer greater than 100")
  }

  new("vglmff",
  blurb = c("Bivariate probit model with rho = ", format(rho), "\n",
            "Links:    ",
            namesof("mu1", lmu12, earg = emu12), ", ",
            namesof("mu2", lmu12, earg = emu12)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(c(1, 1), 2, 1), x = x,
                           bool = .exchangeable ,
                           constraints = constraints,
                           apply.int = TRUE)
  }), list( .exchangeable = exchangeable ))),
  deviance = Deviance.categorical.data.vgam,
  initialize = eval(substitute(expression({
    eval(process.binomial2.data.VGAM)
    predictors.names <- c(
                  namesof("mu1", .lmu12 , earg = .emu12 , short = TRUE),
                  namesof("mu2", .lmu12 , earg = .emu12 , short = TRUE))

    if (is.null( .nsimEIM )) {
         save.weights <- control$save.weights <- FALSE
    }
    if (is.null(etastart)) {
      mu1.init= if (is.Numeric(.imu1))
                rep(.imu1, length = n) else
                mu[, 3] + mu[, 4]
      mu2.init= if (is.Numeric(.imu2))
                rep(.imu2, length = n) else
                mu[, 2] + mu[, 4]
      etastart <-
        cbind(theta2eta(mu1.init, .lmu12 , earg = .emu12 ),
              theta2eta(mu2.init, .lmu12 , earg = .emu12 ))
    }
  }), list( .lmu12 = lmu12, .emu12 = emu12, .nsimEIM = nsimEIM,
            .imu1 = imu1, .imu2 = imu2 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pmargin <- cbind(eta2theta(eta[, 1], .lmu12 , earg = .emu12 ),
                     eta2theta(eta[, 2], .lmu12 , earg = .emu12 ))
    rhovec <- rep( .rho , len = nrow(eta))
    p11 <- pbinorm(eta[, 1], eta[, 2], cov12 = rhovec)
    p01 <- pmin(pmargin[, 2] - p11, pmargin[, 2])
    p10 <- pmin(pmargin[, 1] - p11, pmargin[, 1])
    p00 <- 1 - p01 - p10 - p11
    ansmat <- abs(cbind("00"=p00, "01"=p01, "10"=p10, "11"=p11))
    ansmat / rowSums(ansmat)
  }, list( .lmu12 = lmu12,
           .emu12 = emu12, .rho = rho ))),
  last = eval(substitute(expression({
      misc$link <-    c(mu1 = .lmu12 , mu2 = .lmu12 )

      misc$earg <- list(mu1 = .emu12 , mu2 = .emu12 )

      misc$nsimEIM <- .nsimEIM
      misc$expected <- TRUE
      misc$rho <- .rho
  }), list( .lmu12 = lmu12,
            .emu12 = emu12,
            .rho = rho, .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ycounts <- if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                  y * w # Convert proportions to counts
      nvec <- if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                round(w)

      smallno <- 1.0e4 * .Machine$double.eps
      if (max(abs(ycounts - round(ycounts))) > smallno)
          warning("converting 'ycounts' to integer in @loglikelihood")
      ycounts <- round(ycounts)

      ll.elts <-
        (if (is.numeric(extra$orig.w)) extra$orig.w else 1) *
         dmultinomial(x = ycounts, size = nvec, prob = mu,
                       log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .rho = rho ))),
  vfamily = c("binom2.Rho", "binom2"),
  deriv = eval(substitute(expression({
    pmargin <- cbind(eta2theta(eta[, 1], .lmu12 , earg = .emu12 ),
                     eta2theta(eta[, 2], .lmu12 , earg = .emu12 ))
    rhovec <- rep( .rho , len = nrow(eta))
    p11 <- pbinorm(eta[, 1], eta[, 2], cov12 = rhovec)
    p01 <- pmargin[, 2]-p11
    p10 <- pmargin[, 1]-p11
    p00 <- 1-p01-p10-p11

    ABmat <- (eta[, 1:2] -
              rhovec * eta[, 2:1]) /  sqrt(pmax(1e5 * .Machine$double.eps,
                                                1.0 - rhovec^2))
    PhiA <- pnorm(ABmat[, 1])
    PhiB <- pnorm(ABmat[, 2])
    onemPhiA <- pnorm(ABmat[, 1], lower.tail = FALSE)
    onemPhiB <- pnorm(ABmat[, 2], lower.tail = FALSE)

    smallno <- 1000 * .Machine$double.eps
    p00[p00 < smallno] <- smallno
    p01[p01 < smallno] <- smallno
    p10[p10 < smallno] <- smallno
    p11[p11 < smallno] <- smallno

    dprob00 <- dibinorm(eta[, 1], eta[, 2], cov12 = rhovec)
    dl.dprob1 <- PhiB*(y[, 4]/p11-y[, 2]/p01) +
                onemPhiB*(y[, 3]/p10-y[, 1]/p00)
    dl.dprob2 <- PhiA*(y[, 4]/p11-y[, 3]/p10) +
                onemPhiA*(y[, 2]/p01-y[, 1]/p00)
    dprob1.deta <- dtheta.deta(pmargin[, 1], .lmu12 , earg = .emu12 )
    dprob2.deta <- dtheta.deta(pmargin[, 2], .lmu12 , earg = .emu12 )
    dthetas.detas <- cbind(dprob1.deta, dprob2.deta)

    c(w) * cbind(dl.dprob1, dl.dprob2) * dthetas.detas
  }), list( .lmu12 = lmu12, .emu12 = emu12, .rho = rho ))),
  weight = eval(substitute(expression({
    if (is.null( .nsimEIM)) {
      ned2l.dprob1prob1 <- PhiB^2 *(1/p11+1/p01) + onemPhiB^2 *(1/p10+1/p00)
      ned2l.dprob2prob2 <- PhiA^2 *(1/p11+1/p10) + onemPhiA^2 *(1/p01+1/p00)
      ned2l.dprob1prob2 <- PhiA * (PhiB/p11 - onemPhiB/p10) +
                       onemPhiA * (onemPhiB/p00 - PhiB/p01)
      wz <- matrix(0, n, dimm(M))  # 6=dimm(M)
      wz[, iam(1, 1, M)] <- ned2l.dprob1prob1 * dprob1.deta^2
      wz[, iam(2, 2, M)] <- ned2l.dprob2prob2 * dprob2.deta^2
      wz[, iam(1, 2, M)] <- ned2l.dprob1prob2 * dprob1.deta * dprob2.deta
    } else {
      run.varcov <- 0
      ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
      for (ii in 1:( .nsimEIM )) {
        ysim <- rbinom2.rho(n = n, mu1 = pmargin[, 1],
                                   mu2 = pmargin[, 2],
                         twoCols = FALSE, rho = rhovec)
        dl.dprob1 <- PhiB * (ysim[, 4]/p11-ysim[, 2]/p01) +
                             onemPhiB * (ysim[, 3]/p10-ysim[, 1]/p00)
        dl.dprob2 <- PhiA * (ysim[, 4]/p11-ysim[, 3]/p10) +
                             onemPhiA * (ysim[, 2]/p01-ysim[, 1]/p00)
  
        rm(ysim)
        temp3 <- cbind(dl.dprob1, dl.dprob2)
        run.varcov <- ((ii-1) * run.varcov +
                      temp3[, ind1$row.index] *
                      temp3[, ind1$col.index]) / ii
      }
      wz <- if (intercept.only)
        matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz <- wz * dthetas.detas[, ind1$row] *
                 dthetas.detas[, ind1$col]
    }
    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}




 binom2.rho.ss <-
               function(lrho = "rhobit",
                        lmu = "probit",  # added 20120817
                        imu1 = NULL, imu2 = NULL, irho = NULL,
                        imethod = 1,
                        zero = 3,
                        exchangeable = FALSE,
                        grho = seq(-0.95, 0.95, by = 0.05)) {



  lrho <- as.list(substitute(lrho))
  e.rho <- link2list(lrho)
  l.rho <- attr(e.rho, "function.name")

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  if (lmu != "probit")
    warning("argument 'lmu' should be 'probit'. Changing it.")

    lmu12 <- "probit"  # But emu may contain some arguments.
    emu12 <- emu  # list()


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Bivariate probit model with sample selection\n",
            "Links:    ",
            namesof("mu1", lmu12, earg = emu12), ", ",
            namesof("mu2", lmu12, earg = emu12), ", ",
            namesof("rho", l.rho, earg = e.rho)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(c(1, 1, 0, 0, 0, 1), 3, 2), x = x,
                           bool = .exchangeable ,
                           constraints = constraints,
                           apply.int = TRUE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .exchangeable = exchangeable, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         multipleResponses = FALSE,
         parameters.names = c("mu1", "mu2", "rho"),
         zero = .zero )
  }, list( .zero = zero ))),


  initialize = eval(substitute(expression({

    if (!is.matrix(y))
      stop("response must be a 2- or 3-column matrix")
    ncoly <- ncol(y)

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.min = 1,
              ncol.w.max = 1,
              ncol.y.min = 2,
              ncol.y.max = 3,
              Is.integer.y = TRUE,
              Is.nonnegative.y = TRUE,
              out.wy = TRUE,
              colsyperw = ncoly,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    if (!all(c(y) == 0 | c(y) == 1))
      stop("response matrix must have values 0 and 1 only")


    if (ncoly == 2) {
      extra$ymat2col <- y
      y <- cbind("0"  = 1 - y[, 1],
                 "10" = y[, 1] * (1 - y[, 2]),
                 "11" = y[, 1] *      y[, 2])
    } else {
      if (!all(rowSums(y) == 1))
        stop("response matrix must have two 0s and one 1 in each row")
      y1vec <- 1 - y[, 1]  # Not a 0 means a 1.
      y2vec <- ifelse(y1vec == 1, y[, 3], 0)
      extra$ymat2col <- cbind(y1vec, y2vec)
    }


    predictors.names <- c(
        namesof("mu1", .lmu12 , earg = .emu12 , short = TRUE),
        namesof("mu2", .lmu12 , earg = .emu12 , short = TRUE),
        namesof("rho", .l.rho , earg = .e.rho,  short = TRUE))


    ycounts <- y
    nvec <- 1



    if (!length(etastart)) {
      if (length(mustart)) {
        mu1.init <- mustart[, 1]
        mu2.init <- mustart[, 2]
      } else if ( .imethod == 1) {
        mu1.init <- weighted.mean(extra$ymat2col[, 1], c(w))
        index1 <- (extra$ymat2col[, 1] == 1)
        mu2.init <- weighted.mean(extra$ymat2col[index1, 2], w[index1, 1])
        mu1.init <- rep(mu1.init, len = n)
        mu2.init <- rep(mu2.init, len = n)

      } else if ( .imethod == 2) {
 warning("not working yet2")
          glm1.fit <- glm(ycounts ~ x - 1,
                          weights = c(w),
                          fam = binomial("probit"))
          glm2.fit <- glm(ycounts[, 2:1] ~ x - 1,
                          weights = c(w),
                          fam = binomial("probit"))
          mu1.init <- fitted(glm1.fit)
          mu2.init <- fitted(glm2.fit)
      } else {
        stop("bad value for argument 'imethod'")
      }

      if (length( .imu1 ))
        mu1.init <- rep( .imu1 , length = n)
      if (length( .imu2 ))
        mu2.init <- rep( .imu2 , length = n)



      binom2.rho.ss.Loglikfun <-
          function(rhoval, y, x, w, extraargs) {
          init.mu1 <-    extraargs$initmu1
          init.mu2 <-    extraargs$initmu2
          ymat2col <-    extraargs$ymat2col
          nvec     <-    extraargs$nvec
          eta1 <- qnorm(init.mu1)
          eta2 <- qnorm(init.mu2)

          smallno <- 1000 * .Machine$double.eps
          p11 <- pmax(smallno, pbinorm(eta1, eta2, cov12 = rhoval))
          p10 <- pmax(smallno, pnorm( eta1) - p11)
          p0  <- pmax(smallno, pnorm(-eta1))

          mumat <- abs(cbind("0"  = p0,
                             "10" = p10,
                             "11" = p11))  # rows sum to unity

          smallpos <- 1.0e-100
          mumat[mumat < smallpos] <- smallpos
          ycounts <- y  # n x 3
          use.mu <- mumat  # cbind(p0, p10, p11)

          retval <-
          sum(c(w) *
              dmultinomial(x = ycounts, size = nvec, prob = use.mu,  # mumat,
                           log = TRUE, dochecking = FALSE))
          retval
        }
        rho.grid <- .grho  # seq(-0.95, 0.95, len = 31)
        try.this <- grid.search(rho.grid, objfun = binom2.rho.ss.Loglikfun,
                                y = y, x = x, w = w, extraargs = list(
                                ymat2col = extra$ymat2col,
                                initmu1  = mu1.init,
                                initmu2  = mu2.init,
                                nvec = nvec ))

      rho.init <- if (is.Numeric( .irho ))
                   rep( .irho , len = n) else {
          try.this
      }



      etastart <- cbind(theta2eta(mu1.init, .lmu12 , earg = .emu12 ),
                        theta2eta(mu2.init, .lmu12 , earg = .emu12 ),
                        theta2eta(rho.init, .l.rho , earg = .e.rho ))
    }
    mustart <- NULL  # Since etastart has been computed and/or no @linkfun.
  }), list( .lmu12 = lmu12, .l.rho = l.rho,
            .emu12 = emu12, .e.rho = e.rho, 
                            .grho = grho,
                            .irho = irho,
            .imethod = imethod,
            .imu1 = imu1, .imu2 = imu2 ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    rhovec <- eta2theta(eta[, 3], .l.rho , earg = .e.rho )

    smallno <- 1000 * .Machine$double.eps
    p11 <- pmax(smallno, pbinorm(eta[, 1], eta[, 2], cov12 = rhovec))
    p10 <- pmax(smallno, pnorm( eta[, 1]) - p11)
    p0  <- pmax(smallno, pnorm(-eta[, 1]))
    sumprob <- p11 + p10 + p0
    p11 <- p11 / sumprob
    p10 <- p10 / sumprob
    p0  <- p0  / sumprob

    ansmat <- abs(cbind("0"  = p0,  # p0 == P(Y_1 = 0)
                        "10" = p10,
                        "11" = p11))
    ansmat
  }, list( .lmu12 = lmu12, .l.rho = l.rho,
           .emu12 = emu12, .e.rho = e.rho ))),
  last = eval(substitute(expression({
    misc$link <-    c(mu1 = .lmu12 , mu2 = .lmu12 , rho = .l.rho )

    misc$earg <- list(mu1 = .emu12 , mu2 = .emu12 , rho = .e.rho )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lmu12 = lmu12, .l.rho = l.rho,
            .emu12 = emu12, .e.rho = e.rho ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      ycounts <- y  # n x 3
      nvec <- 1

      smallno <- 1000 * .Machine$double.eps
      rhovec <- eta2theta(eta[, 3], .l.rho , earg = .e.rho )
      p11 <- pmax(smallno, pbinorm(eta[, 1], eta[, 2], cov12 = rhovec))
      p10 <- pmax(smallno, pnorm( eta[, 1]) - p11)
      p0  <- pmax(smallno, pnorm(-eta[, 1]))
      sumprob <- p11 + p10 + p0
      p11 <- p11 / sumprob
      p10 <- p10 / sumprob
      p0  <- p0  / sumprob


      ll.elts <-
        c(w) * dmultinomial(x = ycounts, size = nvec, prob = mu,  # use.mu,
                            log = TRUE, dochecking = FALSE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .l.rho = l.rho, .e.rho = e.rho ))),
  vfamily = c("binom2.rho.ss", "binom2"),

  deriv = eval(substitute(expression({
    nvec <- 1
    ycounts <- extra$ymat2col

    pmargin <- cbind(eta2theta(eta[, 1], .lmu12 , earg = .emu12 ),
                     eta2theta(eta[, 2], .lmu12 , earg = .emu12 ))
    rhovec <-        eta2theta(eta[, 3], .l.rho , earg = .e.rho )


    smallno <- 1000 * .Machine$double.eps
    p11 <- pmax(smallno, pbinorm(eta[, 1], eta[, 2], cov12 = rhovec))
    p10 <- pmax(smallno, pnorm( eta[, 1]) - p11)
    p0  <- pmax(smallno, pnorm(-eta[, 1]))
    sumprob <- p11 + p10 + p0
    p11 <- p11 / sumprob
    p10 <- p10 / sumprob
    p0  <- p0  / sumprob


    BAmat <- (eta[, 1:2] -
              rhovec * eta[, 2:1]) /  sqrt(pmax(1e5 * .Machine$double.eps,
                                                1.0 - rhovec^2))


    PhiA     <- pnorm(BAmat[, 2])
    PhiB     <- pnorm(BAmat[, 1])
    onemPhiA <- pnorm(BAmat[, 2], lower.tail = FALSE)
    onemPhiB <- pnorm(BAmat[, 1], lower.tail = FALSE)


  mycode <- FALSE  # zz
  mycode <- TRUE   # zz

 if (mycode) {
    dprob00 <- dibinorm(eta[, 1], eta[, 2], cov12 = rhovec)
    dl.dprob1 <-     PhiA *      ycounts[, 1] *      ycounts[, 2]  / p11 +
                 onemPhiA *      ycounts[, 1] * (1 - ycounts[, 2]) / p10 -
                            (1 - ycounts[, 1]) / p0
    dl.dprob2 <-     PhiB * (    ycounts[, 1] *      ycounts[, 2]  / p11 -
                                 ycounts[, 1] * (1 - ycounts[, 2]) / p10)
    dl.drho   <-  dprob00 * (    ycounts[, 1] *      ycounts[, 2]  / p11 -
                                 ycounts[, 1] * (1 - ycounts[, 2]) / p10)

    dprob1.deta <- dtheta.deta(pmargin[, 1], .lmu12 , earg = .emu12 )
    dprob2.deta <- dtheta.deta(pmargin[, 2], .lmu12 , earg = .emu12 )
    drho...deta <- dtheta.deta(rhovec,       .l.rho , earg = .e.rho )

    ans.deriv <- c(w) * cbind(dl.dprob1 * dprob1.deta,
                              dl.dprob2 * dprob2.deta,
                              dl.drho   * drho...deta)
 }  # else {
    eta1 <- eta[, 1]  # dat1 %*% params[1:X1.d2]
    eta2 <- eta[, 2]  # dat2 %*% params[(X1.d2 + 1):(X1.d2 + X2.d2)]
    corr.st <- eta[, 3]  # params[(X1.d2 + X2.d2 + 1)]
    corr <- rhovec # tanh(corr.st)

    dat <- ycounts

    y1.y2  <-      dat[, 1] *      dat[, 2]
    y1.cy2 <-      dat[, 1] * (1 - dat[, 2])
    cy1    <- (1 - dat[, 1])

    d.r <- 1/sqrt(pmax(10000 * .Machine$double.eps, 1 - corr^2))
    A <- pnorm((eta2 - corr * eta1) * d.r)
    A.c <- 1 - A
    B <- pnorm((eta1 - corr * eta2) * d.r)
    p11 <- pmax(pbinorm(eta1, eta2, cov12 = corr), 1000 * .Machine$double.eps)
    p10 <- pmax(pnorm( eta1) - p11, 1000 * .Machine$double.eps)
    p0  <- pmax(pnorm(-eta1), 1000 * .Machine$double.eps)
    d.n1 <- dnorm(eta1)
    d.n2 <- dnorm(eta2)
    d.n1n2 <- dibinorm(eta1, eta2, cov12 = corr)
    drh.drh.st <- 4 * exp(2 * corr.st)/(exp(2 * corr.st) + 1)^2

    dl.dbe1 <- d.n1 * (y1.y2/p11 * A + y1.cy2/p10 * A.c - cy1/p0)
    dl.dbe2 <- d.n2 * B * (y1.y2/p11 - y1.cy2/p10)
    dl.drho <- d.n1n2 * (y1.y2/p11 - y1.cy2/p10) * drh.drh.st

    ans.deriv2 <- c(w) * cbind(dl.dbe1, dl.dbe2, dl.drho)
 # }










    ans.deriv
  }), list( .lmu12 = lmu12, .l.rho = l.rho,
            .emu12 = emu12, .e.rho = e.rho ))),

  weight = eval(substitute(expression({



 if (mycode) {
    ned2l.dprob1prob1 <-      PhiA^2 / p11 +
                          onemPhiA^2 / p10 +
                                   1 / p0
    ned2l.dprob2prob2 <-   (1/p11 + 1/p10) * PhiB^2
    ned2l.drho2       <-   (1/p11 + 1/p10) * dprob00^2

    ned2l.dprob1prob2 <-    PhiA * PhiB / p11 - onemPhiA * PhiB / p10
    ned2l.dprob1rho   <-   (PhiA/p11 -  onemPhiA/p10) * dprob00
    ned2l.dprob2rho   <-   (1/p11 + 1/p10) * PhiB * dprob00

    wz <- matrix(0, n, dimm(M))  # 6=dimm(M)
    wz[, iam(1, 1, M)] <- ned2l.dprob1prob1 * dprob1.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dprob2prob2 * dprob2.deta^2
    wz[, iam(3, 3, M)] <- ned2l.drho2       * drho...deta^2
    wz[, iam(1, 2, M)] <- ned2l.dprob1prob2 * dprob1.deta * dprob2.deta
    wz[, iam(1, 3, M)] <- ned2l.dprob1rho   * dprob1.deta * drho...deta
    wz[, iam(2, 3, M)] <- ned2l.dprob2rho   * dprob2.deta * drho...deta
  }  # else {

    ned2l.be1.be1 <- (A^2/p11 + A.c^2/p10 + 1/p0)      * d.n1^2
    ned2l.be2.be2 <- (  1/p11 +     1/p10) * B^2       * d.n2^2
    ned2l.rho.rho <- (  1/p11 +     1/p10) * d.n1n2^2  * drh.drh.st^2

    ned2l.be1.be2 <- (A *  B/p11  - A.c *  B/p10)  * d.n1   * d.n2
    ned2l.be1.rho <- (A * (1/p11) - A.c * (1/p10)) * d.n1n2 * d.n1 * drh.drh.st
    ned2l.be2.rho <-  B * (1/p11  +        1/p10)  * d.n1n2 * d.n2 * drh.drh.st




    WZ <- matrix(0, n, dimm(M))  # 6=dimm(M)
    WZ[, iam(1, 1, M)] <- ned2l.be1.be1
    WZ[, iam(2, 2, M)] <- ned2l.be2.be2
    WZ[, iam(3, 3, M)] <- ned2l.rho.rho
    WZ[, iam(1, 2, M)] <- ned2l.be1.be2
    WZ[, iam(1, 3, M)] <- ned2l.be1.rho
    WZ[, iam(2, 3, M)] <- ned2l.be2.rho

    c(w) * wz
  }), list( .zero = zero ))))
}




