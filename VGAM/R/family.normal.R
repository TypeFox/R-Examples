# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.











VGAM.weights.function <- function(w, M, n) {


  ncolw <- ncol(as.matrix(w))
  if (ncolw == 1) {
    wz <- matrix(w, nrow = n, ncol = M)  # w_i * diag(M)
  } else if (ncolw == M) {
    wz <- as.matrix(w)
  } else if (ncolw < M && M > 1) {
    stop("ambiguous input for 'weights'")
  } else if (ncolw > M*(M+1)/2) {
    stop("too many columns")
  } else {
    wz <- as.matrix(w)
  }
  wz
}












 gaussianff <- function(dispersion = 0, parallel = FALSE, zero = NULL) {

  if (!is.Numeric(dispersion, length.arg = 1) ||
      dispersion < 0)
    stop("bad input for argument 'dispersion'")
  estimated.dispersion <- dispersion == 0


  new("vglmff",
  blurb = c("Vector linear/additive model\n",
            "Links:    identitylink for Y1,...,YM"),

  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel , 
                           constraints = constraints)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel = parallel, .zero = zero ))),

  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M <- if (is.matrix(y)) ncol(y) else 1
    n <- if (is.matrix(y)) nrow(y) else length(y)
    wz <- VGAM.weights.function(w = w, M = M, n = n)
    if (residuals) {
      if (M > 1) {
        U <- vchol(wz, M = M, n = n) 
        temp <- mux22(U, y-mu, M = M, upper = TRUE, as.matrix = TRUE)
        dimnames(temp) <- dimnames(y)
        temp
      } else (y-mu) * sqrt(wz)
    } else {
      ResSS.vgam(y-mu, wz = wz, M = M)
    }
  },

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    if (is.R())
      assign("CQO.FastAlgorithm", TRUE, envir = VGAM::VGAMenv) else
      CQO.FastAlgorithm <<- TRUE
    if (any(function.name == c("cqo", "cao")) &&
       (length( .zero ) ||
       (is.logical( .parallel ) && .parallel )))
        stop("cannot handle non-default arguments for cqo() and cao()")


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    M <- if (is.matrix(y)) ncol(y) else 1
    dy <- dimnames(y)

    predictors.names <- if (!is.null(dy[[2]])) dy[[2]] else
                       paste("Y", 1:M, sep = "")

    if (!length(etastart)) 
      etastart <- 0 * y
  }), list( .parallel = parallel, .zero = zero ))),
  linkinv = function(eta, extra = NULL) eta, 
  last = eval(substitute(expression({
    dy <- dimnames(y)
    if (!is.null(dy[[2]]))
        dimnames(fit$fitted.values) <- dy
    dpar <- .dispersion
    if (!dpar) {
      wz <- VGAM.weights.function(w = w, M = M, n = n)
      temp5 <- ResSS.vgam(y-mu, wz = wz, M = M)
        dpar <- temp5 / (length(y) -
        (if (is.numeric(ncol(X.vlm.save))) ncol(X.vlm.save) else 0))
    }
    misc$dispersion <- dpar
    misc$default.dispersion <- 0
    misc$estimated.dispersion <- .estimated.dispersion

    misc$link <- rep("identitylink", length = M)
    names(misc$link) <- predictors.names

    misc$earg <- vector("list", M)
    for (ilocal in 1:M)
      misc$earg[[ilocal]] <- list()
    names(misc$link) <- predictors.names


    if (is.R()) {
      if (exists("CQO.FastAlgorithm", envir = VGAM::VGAMenv))
        rm("CQO.FastAlgorithm", envir = VGAM::VGAMenv)
    } else {
      while (exists("CQO.FastAlgorithm"))
        remove("CQO.FastAlgorithm")
    }

    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .dispersion = dispersion,
            .estimated.dispersion = estimated.dispersion ))),
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M <- if (is.matrix(y)) ncol(y) else 1
    n <- if (is.matrix(y)) nrow(y) else length(y)
    wz <- VGAM.weights.function(w = w, M = M, n = n)


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      temp1 <- ResSS.vgam(y-mu, wz = wz, M = M)




      ll.elts <-
      if (M == 1 || ncol(wz) == M) {

        -0.5 * temp1 + 0.5 *    (log(wz)) - n * (M / 2) * log(2*pi)
      } else {
        if (all(wz[1, ] == apply(wz, 2, min)) &&
            all(wz[1, ] == apply(wz, 2, max))) {
          onewz <- m2a(wz[1, , drop = FALSE], M = M)
          onewz <- onewz[,, 1]  # M x M


          logdet <- determinant(onewz)$modulus
          logretval <- -0.5 * temp1 + 0.5 * n * logdet -
                       n * (M / 2) * log(2*pi)

        distval <- stop("variable 'distval' not computed yet")
        logretval <- -(ncol(onewz) * log(2 * pi) + logdet + distval)/2
        logretval
      } else {
        logretval <- -0.5 * temp1 - n * (M / 2) * log(2*pi)
        for (ii in 1:n) {
          onewz <- m2a(wz[ii, , drop = FALSE], M = M)
          onewz <- onewz[,, 1]  # M x M
          logdet <- determinant(onewz)$modulus
            logretval <- logretval + 0.5 * logdet
          }
          logretval
        }
      }

      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  linkfun = function(mu, extra = NULL) mu,
  vfamily = "gaussianff",
  deriv = expression({
    wz <- VGAM.weights.function(w = w, M = M, n = n)
    mux22(cc = t(wz), xmat = y-mu, M = M, as.matrix = TRUE)
  }),
  weight = expression({
    wz
  }))
}










dposnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(mean), length(sd))
  if (length(x)    != L) x    <- rep(x,    len = L)
  if (length(mean) != L) mean <- rep(mean, len = L)
  if (length(sd)   != L) sd   <- rep(sd,   len = L)

  if (log.arg) {
    ifelse(x < 0, log(0), dnorm(x, mean = mean, sd = sd, log = TRUE) -
           pnorm(mean / sd, log.p = TRUE))
  } else {
    ifelse(x < 0, 0, dnorm(x = x, mean = mean, sd = sd) / pnorm(mean / sd))
  }
}



pposnorm <- function(q, mean = 0, sd = 1,
                     lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  ans <- (pnorm(q, mean = mean, sd = sd) -
          pnorm(0, mean = mean, sd = sd)) / pnorm(mean / sd)
  ans[q <= 0] <- 0

  if (lower.tail) {
    if (log.p) log(ans) else ans
  } else {
    if (log.p) log1p(-ans) else 1-ans
  }
}



qposnorm <- function(p, mean = 0, sd = 1,
                     lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)   # 20150102 KaiH
  
  if (lower.tail) {
    if (log.arg) p <- exp(p)
  } else {
    p <- if (log.arg) -expm1(p) else 1 - p
  }

  qnorm(p = p + (1 - p) * pnorm(0, mean = mean, sd = sd),
        mean = mean, sd = sd)
}



rposnorm <- function(n, mean = 0, sd = 1) {
  qnorm(p = runif(n, min = pnorm(0, mean = mean, sd = sd)),
        mean = mean, sd = sd)
}




if (FALSE)
 posnormal.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}




 posnormal <- function(lmean = "identitylink", lsd = "loge",
                       eq.mean = FALSE, eq.sd = FALSE,
                       gmean = exp((-5:5)/2), gsd = exp((-1:5)/2),
                       imean = NULL, isd = NULL, probs.y = 0.10,
                       imethod = 1,
                       nsimEIM = NULL, zero = "sd") {





  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")

  if (!is.logical(eq.mean) || length(eq.mean) != 1)
    stop("bad input for argument 'eq.mean'")
  if (!is.logical(eq.sd  ) || length(eq.sd  ) != 1)
    stop("bad input for argument 'eq.sd'")

  if (length(isd) &&
      !is.Numeric(isd, positive = TRUE))
    stop("bad input for argument 'isd'")


  if (length(nsimEIM))
    if (!is.Numeric(nsimEIM, length.arg = 1,
                    integer.valued = TRUE) ||
        nsimEIM <= 10)
      stop("argument 'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("Positive (univariate) normal distribution\n\n",
          "Links:    ",
          namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
          namesof("sd",   lsd,   earg = esd,   tag = TRUE)),



  constraints = eval(substitute(expression({


    constraints.orig <- constraints
    M1 <- 2
    NOS <- M / M1

    cm1.m <-
    cmk.m <- kronecker(diag(NOS), rbind(1, 0))
    con.m <- cm.VGAM(kronecker(matrix(1, NOS, 1), rbind(1, 0)),
                     x = x,
                     bool = .eq.mean ,  #
                     constraints = constraints.orig,
                     apply.int = TRUE,
                     cm.default           = cmk.m,
                     cm.intercept.default = cm1.m)


    cm1.s <-
    cmk.s <- kronecker(diag(NOS), rbind(0, 1))
    con.s <- cm.VGAM(kronecker(matrix(1, NOS, 1), rbind(0, 1)),
                     x = x,
                     bool = .eq.sd ,  #
                     constraints = constraints.orig,
                     apply.int = TRUE,
                     cm.default           = cmk.s,
                     cm.intercept.default = cm1.s)


    con.use <- con.m
    for (klocal in 1:length(con.m)) {


      con.use[[klocal]] <- interleave.cmat(con.m[[klocal]], con.s[[klocal]])

    }
    constraints <- con.use

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)


  }), list( .zero    = zero,
            .eq.sd   = eq.sd,
            .eq.mean = eq.mean ))),



  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         eq.mean = .eq.mean ,
         eq.sd   = .eq.sd   ,
         multipleResponses = TRUE,
         parameters.names = c("mean", "sd"),
         zero = .zero )
  }, list( .zero = zero,
           .eq.mean = eq.mean,
           .eq.sd   = eq.sd
         ))),


  initialize = eval(substitute(expression({
    M1 <- 2
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
    NOS <- ncol(y)
    M <- NOS * M1

    mean.names  <- param.names("mean",     NOS)
    sdev.names  <- param.names("sd",       NOS)

    predictors.names <-
      c(namesof(mean.names , .lmean     , earg = .emean     , tag = FALSE),
        namesof(sdev.names , .lsd       , earg = .esd       , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]

    if (!length(etastart)) {
      init.me <- matrix( if (length( .imean )) .imean else NA_real_,
                        n, NOS, byrow = TRUE)
      init.sd <-  matrix( if (length( .isd  )) .isd   else NA_real_,
                        n, NOS, byrow = TRUE)

      mean.grid.orig <- .gmean
      sdev.grid.orig <- .gsd


      for (jay in 1:NOS) {
        yvec <- y[, jay]
        wvec <- w[, jay]
        if (any(is.na(init.me[, jay]))) {
          init.me[, jay] <- if ( .imethod == 1) {
            weighted.mean(yvec, wvec)
          } else if ( .imethod == 2) {
            quantile(yvec, probs = .probs.y )
          } else if ( .imethod == 3) {
            median(yvec)
          }
        }
        if (any(is.na(init.sd[, jay])))
          init.sd[, jay] <- sd(yvec)


        ll.posnormal <- function(sdev.val, y, x, w, extraargs) {
          ans <-
          sum(c(w) * dposnorm(x = y, mean = extraargs$Mean,
                              sd = sdev.val, log = TRUE))
          ans
        }


        sdev.grid <- sdev.grid.orig * init.sd[1, jay]
        mean.grid <- mean.grid.orig * init.me[1, jay]
        mean.grid <- sort(c(-mean.grid,
                             mean.grid))
        allmat1 <- expand.grid(Mean = mean.grid)
        allmat2 <- matrix(NA_real_, nrow(allmat1), 2)

         for (iloc in 1:nrow(allmat1)) {
            allmat2[iloc, ] <-
              grid.search(sdev.grid, objfun = ll.posnormal,
                           y = yvec, x = x, w = wvec,
                           ret.objfun = TRUE,  # 2nd value is the loglik
                           extraargs = list(Mean = allmat1[iloc, "Mean"]))
         }
        ind5 <- which.max(allmat2[, 2])  # 2nd value is the loglik

        if (!length( .imean ))
          init.me[, jay] <- allmat1[ind5, "Mean"]
        if (!length( .isd   ))
          init.sd[, jay] <- allmat2[ind5, 1]
      }  # jay




      etastart <- cbind(theta2eta(init.me, .lmean , earg = .emean ),
                        theta2eta(init.sd, .lsd ,   earg = .esd   ))
      etastart <- etastart[, interleave.VGAM(M, M1 = M1)]

    }
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .gmean = gmean, .gsd = gsd,
            .imean = imean, .isd = isd,
            .imethod = imethod, .probs.y = probs.y
           ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mymu <- eta2theta(eta[, c(TRUE, FALSE)], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, c(FALSE, TRUE)], .lsd   , earg = .esd  )
    mymu + mysd * dnorm(-mymu/mysd) / pnorm(mymu/mysd)
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd
         ))),
  last = eval(substitute(expression({
    misc$link <- c(rep( .lmean , length = NOS),
                   rep( .lsd   , length = NOS))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mean.names, sdev.names)
    names(misc$link) <- temp.names[interleave.VGAM(M, M1 = M1)]

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .emean
      misc$earg[[M1*ii  ]] <- .esd
    }

    misc$expected          <- TRUE
    misc$multipleResponses <- TRUE

    misc$nsimEIM <- .nsimEIM
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mymu <- eta2theta(eta[, c(TRUE, FALSE)], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, c(FALSE, TRUE)], .lsd   , earg = .esd  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dposnorm(x = y, m = mymu, sd = mysd, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd ))),
  vfamily = c("posnormal"),



  simslot = eval(substitute(
  function(object, nsim) {
    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    mymu <- eta2theta(eta[, c(TRUE, FALSE)], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, c(FALSE, TRUE)], .lsd   , earg = .esd   )
    rposnorm(nsim * length(mymu), mean = mymu, sd = mysd)
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd ))),







  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta[, c(TRUE, FALSE)], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, c(FALSE, TRUE)], .lsd   , earg = .esd   )


    zedd <- (y-mymu) / mysd
    temp0 <-   mymu  / mysd
    imratio <- dnorm(temp0) / pnorm(temp0)

    dl.dmu <- (zedd - imratio) / mysd
    dl.dsd <- (temp0 * imratio + zedd^2 - 1) / mysd

    dmu.deta <- dtheta.deta(mymu, .lmean , earg = .emean )
    dsd.deta <- dtheta.deta(mysd, .lsd   , earg = .esd   )
    dthetas.detas <- cbind(dmu.deta, dsd.deta)
    myderiv <- c(w) * dthetas.detas * cbind(dl.dmu, dl.dsd)
    myderiv <- myderiv[, interleave.VGAM(M, M1 = M1)]
    myderiv
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    if (length( .nsimEIM )) {



      NOS <- M / M1
      dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M1 = M1)]

      wz <- matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal' 

      ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)

      for (spp. in 1:NOS) {
        run.varcov <- 0
        Mymu <- mymu[, spp.]
        Mysd <- mysd[, spp.]

      for (ii in 1:( .nsimEIM )) {
        ysim <- rposnorm(n, m = Mymu, sd = Mysd)


        zedd <- (ysim-Mymu) / Mysd
        dl.dmu <- (zedd - imratio) / Mysd
        dl.dsd <- (temp0 * imratio + zedd^2 - 1) / Mysd

        
        temp7 <- cbind(dl.dmu, dl.dsd)
        run.varcov <- run.varcov +
                      temp7[, ind1$row.index] *
                      temp7[, ind1$col.index]
      }
      run.varcov <- cbind(run.varcov / .nsimEIM )



      wz1 <- if (intercept.only)
          matrix(colMeans(run.varcov),
                 nrow = n, ncol = ncol(run.varcov), byrow = TRUE) else
          run.varcov

      wz1 <- wz1 * dThetas.detas[, M1 * (spp. - 1) + ind1$row] *
                   dThetas.detas[, M1 * (spp. - 1) + ind1$col]


      for (jay in 1:M1)
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = M1)]
        }
      }  # End of for (spp.) loop



      wz <- w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = M / M1)


    } else {

      ned2l.dmu2 <- (1 - imratio * (temp0 + imratio)) / mysd^2
      ned2l.dmusd <- imratio * (1 + temp0 * (temp0 + imratio)) / mysd^2
      ned2l.dsd2 <- (2 - imratio * (temp0 * (1 + temp0 *
                    (temp0 + imratio)))) / mysd^2
  
      wz <- array(c(c(w) * ned2l.dmu2  * dmu.deta^2,
                    c(w) * ned2l.dsd2  * dsd.deta^2,
                    c(w) * ned2l.dmusd * dmu.deta * dsd.deta),
                  dim = c(n, M/M1, M1*(M1+1)/2))
      wz <- arwz2wz(wz, M = M, M1 = M1)
    }
    wz
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .nsimEIM = nsimEIM ))))
}





dbetanorm <- function(x, shape1, shape2, mean = 0, sd = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  logden <-
    dnorm(x = x, mean = mean, sd = sd, log = TRUE) +
    (shape1-1) * pnorm(q = x, mean = mean, sd = sd, log.p = TRUE) +
    (shape2-1) * pnorm(q = x, mean = mean, sd = sd, log.p = TRUE,
                       lower.tail = FALSE) -
    lbeta(shape1, shape2)

  logden[is.infinite(x)] <- log(0)  # 20141210 KaiH
  if (log.arg) logden else exp(logden)
}




pbetanorm <- function(q, shape1, shape2, mean = 0, sd = 1,
                      lower.tail = TRUE, log.p = FALSE) {
  pbeta(q = pnorm(q = q, mean = mean, sd = sd),
        shape1 = shape1, shape2 = shape2,
        lower.tail = lower.tail, log.p = log.p)
}



qbetanorm <- function(p, shape1, shape2, mean = 0, sd = 1,
                      lower.tail = TRUE, log.p = FALSE) {
  qnorm(p = qbeta(p = p, shape1 = shape1, shape2 = shape2,
                  lower.tail = lower.tail, log.p = log.p),
        mean = mean, sd = sd)
}



rbetanorm <- function(n, shape1, shape2, mean = 0, sd = 1) {
  qnorm(p = qbeta(p = runif(n), shape1 = shape1, shape2 = shape2),
        mean = mean, sd = sd)
}




dtikuv <- function(x, d, mean = 0, sigma = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(d, length.arg = 1) ||
      max(d) >= 2)
    stop("bad input for argument 'd'")

  L <- max(length(x), length(mean), length(sigma))
  if (length(x)     != L) x     <- rep(x,     len = L)
  if (length(mean)  != L) mean  <- rep(mean,  len = L)
  if (length(sigma) != L) sigma <- rep(sigma, len = L)


  hh <- 2 - d
  KK <- 1 / (1 + 1/hh + 0.75/hh^2)
  logden <- dnorm(x = x, mean = mean, sd = sigma, log = TRUE) + log(KK) +
    2 * log1p(((x-mean)/sigma)^2 / (2*hh))
  logden[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) logden else exp(logden)
}



ptikuv <- function(q, d, mean = 0, sigma = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(d, length.arg = 1) ||
      max(d) >= 2)
    stop("bad input for argument 'd'")

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)  # 20141231 KaiH

  L <- max(length(q), length(mean), length(sigma))
  if (length(q)     != L) q     <- rep(q,     len = L)
  if (length(mean)  != L) mean  <- rep(mean,  len = L)
  if (length(sigma) != L) sigma <- rep(sigma, len = L)

  zedd1 <- 0.5 * ((q - mean) / sigma)^2
  ans <- q*0 + 0.5
  hh <- 2 - d
  KK <- 1 / (1 + 1/hh + 0.75/hh^2)
  if (any(lhs <- q < mean)) {
    ans[lhs] <- ( KK/(2*sqrt(pi))) * (
    gamma(0.5) * (1 - pgamma(zedd1[lhs], 0.5)) +
    2 * gamma(1.5) * (1 - pgamma(zedd1[lhs], 1.5)) / hh +
    gamma(2.5) * (1 - pgamma(zedd1[lhs], 2.5)) / hh^2)
  }
  if (any(rhs <- q > mean)) {
    ans[rhs] <- 1.0 - Recall(q = (2*mean[rhs] - q[rhs]), d = d,
                             mean = mean[rhs], sigma = sigma[rhs])
  }

  if (lower.tail) {
    if (log.arg) log(ans) else ans
  } else {
    if (log.arg) log1p(-ans) else 1 - ans
  }
}




qtikuv <- function(p, d, mean = 0, sigma = 1, 
                   lower.tail = TRUE, log.p = FALSE, ...) {
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (!is.Numeric(d, length.arg = 1) || max(d) >= 2)
    stop("bad input for argument 'd'")

  orig.p <- p
  if (lower.tail) {
    if (log.p) p <- exp(p)
  } else {
    p <- if (log.p) -expm1(p) else 1 - p
  }
  
  L <- max(length(p), length(mean), length(sigma))
  if (length(p)     != L) p     <- rep(p,     len = L)
  if (length(mean)  != L) mean  <- rep(mean,  len = L)
  if (length(sigma) != L) sigma <- rep(sigma, len = L)
  ans <- rep(0.0, len = L)


  myfun <- function(x, d, mean = 0, sigma = 1, p)
    ptikuv(q = x, d = d, mean = mean, sigma = sigma) - p

  for (ii in 1:L) {
    Lower <- ifelse(p[ii] <= 0.5, mean[ii] - 3 * sigma[ii], mean[ii])
    while (ptikuv(q = Lower, d = d, mean = mean[ii],
                  sigma = sigma[ii]) > p[ii])
      Lower <- Lower - sigma[ii]
    Upper <- ifelse(p[ii] >= 0.5, mean[ii] + 3 * sigma[ii], mean[ii])
    while (ptikuv(q = Upper, d = d, mean = mean[ii],
                  sigma = sigma[ii]) < p[ii])
      Upper <- Upper + sigma[ii]
    ans[ii] <- uniroot(f = myfun, lower = Lower, upper = Upper,
                       d = d, p = p[ii],
                       mean = mean[ii], sigma = sigma[ii], ...)$root
  }


  if (log.p) {
    ans[orig.p > 0] <- NaN
  } else {
    ans[orig.p < 0] <- NaN
    ans[orig.p > 1] <- NaN
  }

  ans
}


rtikuv <- function(n, d, mean = 0, sigma = 1, Smallno = 1.0e-6) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n
  if (!is.Numeric(d, length.arg = 1) || max(d) >= 2)
    stop("bad input for argument 'd'")
  if (!is.Numeric(mean, length.arg = 1))
    stop("bad input for argument 'mean'")
  if (!is.Numeric(sigma, length.arg = 1))
    stop("bad input for argument 'sigma'")
  if (!is.Numeric(Smallno, positive = TRUE, length.arg = 1) ||
      Smallno > 0.01 ||
      Smallno < 2 * .Machine$double.eps)
      stop("bad input for argument 'Smallno'")
  ans <- rep(0.0, len = use.n)

  ptr1 <- 1; ptr2 <- 0
  hh <- 2 - d
  KK <- 1 / (1 + 1/hh + 0.75/hh^2)
  ymax <- ifelse(hh < 2,
                 dtikuv(x = mean + sigma*sqrt(4 - 2*hh),
                        d = d, mean = mean, sigma = sigma),
                 KK / (sqrt(2 * pi) * sigma))
  while (ptr2 < use.n) {
    Lower <- mean - 5 * sigma
    while (ptikuv(q = Lower, d = d, mean = mean, sigma = sigma) > Smallno)
      Lower <- Lower - sigma
    Upper <- mean + 5 * sigma
    while (ptikuv(q = Upper, d = d, mean = mean, sigma = sigma) < 1-Smallno)
      Upper <- Upper + sigma
    x <- runif(2*use.n, min = Lower, max = Upper)
    index <- runif(2*use.n, max = ymax) <
             dtikuv(x, d = d, mean = mean, sigma = sigma)
    sindex <- sum(index)
    if (sindex) {
      ptr2 <- min(use.n, ptr1 + sindex - 1)
      ans[ptr1:ptr2] <- (x[index])[1:(1+ptr2-ptr1)]
      ptr1 <- ptr2 + 1
    }
  }
  ans
}




 tikuv <- function(d, lmean = "identitylink", lsigma = "loge",
                   isigma = NULL, zero = "sigma") {


  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")



  if (!is.Numeric(d, length.arg = 1) || max(d) >= 2)
      stop("bad input for argument 'd'")



  new("vglmff",
  blurb = c("Short-tailed symmetric [Tiku and Vaughan (1999)] ",
            "distribution\n",
          "Link:     ",
          namesof("mean",  lmean,  earg = emean), ", ",
          namesof("sigma", lsigma, earg = esigma),
          "\n", "\n",
          "Mean:     mean"),
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
         parameters.names = c("mean", "sigma"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y)


    predictors.names <- 
      c(namesof("mean",  .lmean  , earg = .emean  , tag = FALSE),
        namesof("sigma", .lsigma , earg = .esigma , tag = FALSE))


    if (!length(etastart)) {
      sigma.init <- if (length( .isigma )) rep( .isigma , length = n) else {
        hh <- 2 - .d
        KK <- 1 / (1 + 1/hh + 0.75/hh^2)
        K2 <- 1 + 3/hh + 15/(4*hh^2)
        rep(sqrt(var(y) / (KK*K2)), len = n)
      }
      mean.init <- rep(weighted.mean(y, w), len = n) 
      etastart <-
        cbind(theta2eta(mean.init,  .lmean  , earg = .emean  ),
              theta2eta(sigma.init, .lsigma , earg = .esigma ))
    }
  }),list( .lmean = lmean, .lsigma = lsigma,
                           .isigma = isigma, .d = d,
           .emean = emean, .esigma = esigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmean , earg = .emean )
  }, list( .lmean = lmean,
           .emean = emean, .esigma = esigma ))),
  last = eval(substitute(expression({
      misc$link <-    c("mean" = .lmean , "sigma"= .lsigma )

      misc$earg <- list("mean" = .emean , "sigma"= .esigma )

      misc$expected <- TRUE
      misc$d <- .d 
  }), list( .lmean = lmean, .lsigma = lsigma, .d = d,
            .emean = emean, .esigma = esigma ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mymu  <- eta2theta(eta[, 1], .lmean  , earg = .emean  )
    sigma <- eta2theta(eta[, 2], .lsigma , earg = .esigma )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtikuv(x = y, d = .d , mean = mymu,
                               sigma = sigma, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmean = lmean, .lsigma = lsigma, .d = d,
           .emean = emean, .esigma = esigma ))),
  vfamily = c("tikuv"),















  deriv = eval(substitute(expression({
    mymu  <- eta2theta(eta[, 1], .lmean ,  earg = .emean )
    sigma <- eta2theta(eta[, 2], .lsigma, earg = .esigma)

    dmu.deta <- dtheta.deta(mymu, .lmean , earg = .emean )
    dsigma.deta <- dtheta.deta(sigma, .lsigma, earg = .esigma)

    zedd <- (y - mymu) / sigma
    hh <- 2 - .d 
    gzedd <- zedd / (1 + 0.5*zedd^2 / hh)

    dl.dmu <- zedd / sigma - 2 * gzedd / (hh*sigma)
    dl.dsigma <- (zedd^2 - 1 - 2 * zedd * gzedd / hh) / sigma

    c(w) * cbind(dl.dmu    * dmu.deta,
                 dl.dsigma * dsigma.deta)
  }), list( .lmean = lmean, .lsigma = lsigma, .d = d,
            .emean = emean, .esigma = esigma ))),
  weight = eval(substitute(expression({
    ayy <- 1 / (2*hh)
    Dnos <- 1 - (2/hh) * (1 - ayy) / (1 + 2*ayy + 3*ayy^2)
    Dstar <- -1 + 3 * (1 + 2*ayy + 11*ayy^2) / (1 + 2*ayy + 3*ayy^2)

    ned2l.dmymu2 <- Dnos / sigma^2
    ned2l.dnu2   <- Dstar / sigma^2

    wz <- matrix(NA_real_, n, M)  # diagonal matrix
    wz[, iam(1, 1, M)] <- ned2l.dmymu2 * dmu.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dnu2 * dsigma.deta^2
    c(w) * wz
  }), list( .lmean = lmean, .lsigma = lsigma,
            .emean = emean, .esigma = esigma ))))
}




dfoldnorm <- function(x, mean = 0, sd = 1, a1 = 1, a2 = 1,
                      log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  ans <- dnorm(x = x/(a1*sd) - mean/sd) / (a1*sd) +
         dnorm(x = x/(a2*sd) + mean/sd) / (a2*sd)
  ans[x < 0] <- 0

  ans[a1 <= 0 | a2 <= 0] <- NA
  ans[sd <= 0] <- NA

  if (log.arg) log(ans) else ans
}



pfoldnorm <- function(q, mean = 0, sd = 1, a1 = 1, a2 = 1,
                      lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  if (lower.tail) {
    if (log.p) {
      ans <- log(pnorm(q =  q/(a1*sd) - mean/sd) - 
                 pnorm(q = -q/(a2*sd) - mean/sd))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- pnorm(q =  q/(a1*sd) - mean/sd) - 
             pnorm(q = -q/(a2*sd) - mean/sd)
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log(pnorm(q =  q/(a1*sd) - mean/sd, lower.tail = FALSE) + 
                 pnorm(q = -q/(a2*sd) - mean/sd))
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- pnorm(q =  q/(a1*sd) - mean/sd, lower.tail = FALSE) + 
             pnorm(q = -q/(a2*sd) - mean/sd)
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  } 
  ans[a1 <= 0 | a2 <= 0] <- NaN
  ans[sd <= 0] <- NaN
  ans
}



qfoldnorm <- function(p, mean = 0, sd = 1, a1 = 1, a2 = 1,
                      lower.tail = TRUE, log.p = FALSE, ...) {

  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)
  
  if (lower.tail) {
    if (log.arg) p <- exp(p)
  } else {
    p <- if (log.arg) -expm1(p) else 1 - p
  }
  
  L <- max(length(p), length(mean), length(sd), length(a1), length(a2))
  if (length(p)    != L) p    <- rep(p,    len = L)
  if (length(mean) != L) mean <- rep(mean, len = L)
  if (length(sd)   != L) sd   <- rep(sd,   len = L)
  if (length(a1)   != L) a1   <- rep(a1,   len = L)
  if (length(a2)   != L) a2   <- rep(a2,   len = L)
  ans  <- rep(0.0 , len = L)

  myfun <- function(x, mean = 0, sd = 1, a1 = 1, a2 = 2, p)
    pfoldnorm(q = x, mean = mean, sd = sd, a1 = a1, a2 = a2) - p

  for (ii in 1:L) {
    mytheta <- mean[ii] / sd[ii]
    EY <- sd[ii] * ((a1[ii] + a2[ii]) *
          (mytheta * pnorm(mytheta) + dnorm(mytheta)) -
          a2[ii] * mytheta)
    Upper <- 2 * EY
    while (pfoldnorm(q = Upper, mean = mean[ii], sd = sd[ii],
                     a1 = a1[ii], a2 = a2[ii]) < p[ii])
      Upper <- Upper + sd[ii]
    ans[ii] <- uniroot(f = myfun, lower = 0, upper = Upper,
                       mean = mean[ii], sd = sd[ii],
                       a1 = a1[ii], a2 = a2[ii],
                       p = p[ii], ...)$root
  }

  ans[a1 <= 0 | a2 <= 0] <- NaN
  ans[sd <= 0] <- NaN

  ans
}



rfoldnorm <- function(n, mean = 0, sd = 1, a1 = 1, a2 = 1) {
  X <- rnorm(n, mean = mean, sd = sd)
  ans <- pmax(a1 * X, -a2*X)
  ans[a1 <= 0 | a2 <= 0] <- NA
  ans[sd <= 0] <- NA
  ans
}




 foldnormal <- function(lmean = "identitylink", lsd = "loge",
                        imean = NULL,       isd = NULL,
                        a1 = 1, a2 = 1,
                        nsimEIM = 500, imethod = 1, zero = NULL) {
  if (!is.Numeric(a1, positive = TRUE, length.arg = 1) ||
      !is.Numeric(a2, positive = TRUE, length.arg = 1))
    stop("bad input for arguments 'a1' and 'a2'")
  if (any(a1 <= 0 | a2 <= 0))
    stop("arguments 'a1' and 'a2' must each be a positive value")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")



  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")




  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("argument 'nsimEIM' should be an integer greater than 10")
  if (length(imean) && !is.Numeric(imean))
    stop("bad input for 'imean'")

  if (length(isd) && !is.Numeric(isd, positive = TRUE))
    stop("bad input for 'isd'")


  new("vglmff",
  blurb = c("(Generalized) folded univariate normal distribution\n\n",
            "Link:     ",
            namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
            namesof("sd",   lsd,   earg = esd,   tag = TRUE)),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         a1 = .a1 ,
         a2 = .a2 ,
         multiple.responses = FALSE,
         parameters.names = c("mean", "sd"),
         zero = .zero ,
         nsimEIM = .nsimEIM )
  }, list( .zero = zero,
           .a1 = a1, .a2 = a2,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
        c(namesof("mean", .lmean , earg = .emean, tag = FALSE),
          namesof("sd",   .lsd ,   earg = .esd,   tag = FALSE))

    if (!length(etastart)) {
      junk <- lm.wfit(x = x, y = c(y), w = c(w))


 if (FALSE) {
      if ((ncol(cbind(w)) != 1) || any(w != round(w)))
        stop("'weights' must be a vector or a one-column matrix ",
             "with integer values")
      m1d <- meany <- weighted.mean(y, w)
      m2d <- weighted.mean(y^2, w)
      stddev <- sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
      Ahat <- m1d^2 / m2d
      thetahat <- sqrt(max(1/Ahat -1, 0.1))
      mean.init <- rep(if (length( .imean)) .imean else
                thetahat * sqrt((stddev^2 + meany^2) * Ahat), len = n)
      sd.init <- rep(if (length( .isd)) .isd else
                sqrt((stddev^2 + meany^2) * Ahat), len = n)
}


      stddev <- sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
      meany <- weighted.mean(y, w)
      mean.init <- rep(if (length( .imean )) .imean else
          {if ( .imethod == 1) median(y) else meany}, len = n)
      sd.init <- rep(if (length( .isd )) .isd else
          {if ( .imethod == 1)  stddev else 1.2*sd(c(y))}, len = n)
      etastart <- cbind(theta2eta(mean.init, .lmean , earg = .emean ),
                        theta2eta(sd.init,   .lsd ,   earg = .esd ))
    }
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .imean = imean, .isd = isd,
            .a1 = a1, .a2 = a2, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mymu <- eta2theta(eta[, 1], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, 2], .lsd   , earg = .esd   )
    mytheta <- mymu / mysd
    mysd * (( .a1 + .a2 ) * (mytheta * pnorm(mytheta) +
        dnorm(mytheta)) - .a2 * mytheta)
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd,
           .a1 = a1, .a2 = a2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmean , "sd" = .lsd )

    misc$earg <- list("mu" = .emean , "sd" = .esd )

    misc$multipleResponses <- FALSE
    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$simEIM <- TRUE
    misc$imethod <- .imethod
    misc$a1 <- .a1
    misc$a2 <- .a2
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd,
            .imethod = imethod, .nsimEIM = nsimEIM,
            .a1 = a1, .a2 = a2 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mymu <- eta2theta(eta[, 1], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, 2], .lsd   , earg = .esd   )
    a1vec <- .a1
    a2vec <- .a2
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {

      ll.elts <-
        c(w) * dfoldnorm(y, mean = mymu, sd = mysd,
                         a1 = a1vec, a2 = a2vec, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmean = lmean, .lsd = lsd,
           .emean = emean, .esd = esd,
           .a1 = a1, .a2 = a2 ))),
  vfamily = c("foldnormal"),
  deriv = eval(substitute(expression({
    M1 <- 2
    mymu <- eta2theta(eta[, 1], .lmean , earg = .emean )
    mysd <- eta2theta(eta[, 2], .lsd ,   earg = .esd )

    dmu.deta <- dtheta.deta(mymu, .lmean , earg = .emean )
    dsd.deta <- dtheta.deta(mysd, .lsd ,   earg = .esd )

    a1vec <- .a1
    a2vec <- .a2
    d3 <- deriv3(~ log((exp(-0.5*(y/(a1vec*mysd) - mymu/mysd)^2)/a1vec +
                        exp(-0.5*(y/(a2vec*mysd) +
                           mymu/mysd)^2)/a2vec)/(mysd*sqrt(2*pi))),
                name = c("mymu", "mysd"), hessian = FALSE)
    eval.d3 <- eval(d3)
    dl.dthetas <-  attr(eval.d3, "gradient")  # == cbind(dl.dmu, dl.dsd)
    DTHETA.detas <- cbind(dmu.deta, dsd.deta)
    c(w) * DTHETA.detas * dl.dthetas
  }), list( .lmean = lmean, .lsd = lsd,
            .emean = emean, .esd = esd, .a1 = a1, .a2 = a2 ))),
  weight = eval(substitute(expression({
    de3 <- deriv3(~ log((exp(-0.5*(ysim/(a1vec*mysd) -
                             mymu/mysd)^2)/a1vec +
                        exp(-0.5*(ysim/(a2vec*mysd) +
                             mymu/mysd)^2)/a2vec)/(mysd*sqrt(2*pi))),
                  name = c("mymu", "mysd"), hessian = TRUE)
    run.mean <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- abs(rnorm(n, m = mymu, sd = mysd))
      ysim <- rfoldnorm(n = n, mean = mymu, sd = mysd,
                     a1 = a1vec, a2 = a2vec)
      eval.de3 <- eval(de3)
      d2l.dthetas2 <- attr(eval.de3, "hessian")
      rm(ysim)

      temp3 <- matrix(0, n, dimm(M))
      for (ss in 1:M)
        for (tt in ss:M)
          temp3[, iam(ss,tt, M)] <-  -d2l.dthetas2[, ss,tt]

      run.mean <- ((ii-1) * run.mean + temp3) / ii
    }

    wz <- if (intercept.only)
        matrix(colMeans(run.mean), n, dimm(M), byrow = TRUE) else
        run.mean

    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz <- wz * DTHETA.detas[, index0$row] * DTHETA.detas[, index0$col]

  }), list( .nsimEIM = nsimEIM, .a1 = a1, .a2 = a2 ))))
}





lqnorm.control <- function(trace = TRUE, ...) {
    list(trace = trace)
}





lqnorm <- function(qpower = 2,
                   link = "identitylink",
                   imethod = 1, imu = NULL, ishrinkage = 0.95) {


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  if (!is.Numeric(qpower, length.arg = 1) || qpower <= 1)
    stop("bad input for argument 'qpower'")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
      ishrinkage < 0 ||
      ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")



    new("vglmff",
    blurb = c("Minimizing the q-norm of residuals\n",
              "Links:    ",
              namesof("Y1", link, earg = earg, tag = TRUE)),
    initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    M <- if (is.matrix(y)) ncol(y) else 1
    dy <- dimnames(y)


    predictors.names <- if (!is.null(dy[[2]])) dy[[2]] else
                        paste("mu", 1:M, sep = "")
    predictors.names <- namesof(predictors.names, link = .link,
                                earg = .earg, short = TRUE)


    if (!length(etastart))  {
      meany <- weighted.mean(y, w)
      mean.init <- rep(if (length( .i.mu )) .i.mu else {
        if ( .imethod == 2) median(y) else 
        if ( .imethod == 1) meany else
          .ishrinkage * meany + (1 - .ishrinkage ) * y
      }, len = n)
      etastart <- theta2eta(mean.init, link = .link, earg = .earg)
    }
  }), list( .imethod = imethod, .i.mu = imu,
            .ishrinkage = ishrinkage,
            .link = link, .earg = earg ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu <- eta2theta(eta, link = .link , earg = .earg )
    mu
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    dy <- dimnames(y)
    if (!is.null(dy[[2]]))
        dimnames(fit$fitted.values) = dy
    misc$link <- rep( .link, length = M)
    names(misc$link) <- predictors.names

    misc$earg <- list(mu = .earg)

    misc$qpower <- .qpower
    misc$imethod <- .imethod
    misc$objectiveFunction <- sum( c(w) * (abs(y - mu))^(.qpower) )
  }), list( .qpower = qpower,
            .link = link, .earg = earg,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link = .link, earg = .earg)
  }, list( .link = link, .earg = earg ))),
  vfamily = "lqnorm",
  deriv = eval(substitute(expression({
    dmu.deta <- dtheta.deta(theta=mu, link = .link, earg = .earg )
    myresid <- y - mu
    signresid <- sign(myresid)
    temp2 <- (abs(myresid))^(.qpower-1)
    .qpower * c(w) * temp2 * signresid * dmu.deta
  }), list( .qpower = qpower, .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    temp3 <- (abs(myresid))^(.qpower-2)
    wz <- .qpower * (.qpower - 1) * c(w) * temp3 * dmu.deta^2
    wz
  }), list( .qpower = qpower, .link = link, .earg = earg ))))
}







dtobit <- function(x, mean = 0, sd = 1,
                   Lower = 0, Upper = Inf, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(mean), length(sd),
           length(Lower), length(Upper))
  if (length(x)     != L) x     <- rep(x,     len = L)
  if (length(mean)  != L) mean  <- rep(mean,  len = L)
  if (length(sd)    != L) sd    <- rep(sd,    len = L)
  if (length(Lower) != L) Lower <- rep(Lower, len = L)
  if (length(Upper) != L) Upper <- rep(Upper, len = L)

  if (!all(Lower < Upper, na.rm = TRUE))
    stop("all(Lower < Upper) is not TRUE")

  ans <- dnorm(x = x, mean = mean, sd = sd, log = log.arg)
  ans[x <  Lower] <- if (log.arg) log(0.0) else 0.0
  ans[x >  Upper] <- if (log.arg) log(0.0) else 0.0


  ind3 <- x == Lower
  ans[ind3] <- pnorm(q = Lower[ind3], mean = mean[ind3], sd = sd[ind3],
                     log.p = log.arg)

  ind4 <- x == Upper
  ans[ind4] <- pnorm(q = Upper[ind4], mean = mean[ind4], sd = sd[ind4],
                     lower.tail = FALSE, log.p = log.arg)

  ans
}



ptobit <- function(q, mean = 0, sd = 1, Lower = 0, Upper = Inf,
                   lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail) != 1)
    stop("argument 'lower.tail' must be a single logical")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("argument 'log.p' must be a single logical")


  if (!all(Lower < Upper, na.rm = TRUE))
    stop("all(Lower < Upper) is not TRUE")


  ans <- pnorm(q = q, mean = mean, sd = sd, 
               lower.tail = lower.tail, log.p = log.p)
  ind1 <- (q <  Lower)
  ans[ind1] <- if (lower.tail) ifelse(log.p, log(0.0), 0.0) else
                               ifelse(log.p, log(1.0), 1.0)
  ind2 <- (Upper <= q)
  ans[ind2] <- if (lower.tail) ifelse(log.p, log(1.0), 1.0) else
                               ifelse(log.p, log(0.0), 0.0)
  ans
}




qtobit <- function(p, mean = 0, sd = 1,
                   Lower = 0, Upper = Inf,
                   lower.tail = TRUE, log.p = FALSE) {


  if (!all(Lower < Upper, na.rm = TRUE))
    stop("all(Lower < Upper) is not TRUE")

  # 20150127 KaiH; add lower.tail = lower.tail, log.p = log.p
  ans <- qnorm(p, mean = mean, sd = sd, 
               lower.tail = lower.tail, log.p = log.p)
  pnorm.Lower <- ptobit(q = Lower, mean = mean, sd = sd, 
                        lower.tail = lower.tail, log.p = log.p)
  pnorm.Upper <- ptobit(q = Upper, mean = mean, sd = sd, 
                        lower.tail = lower.tail, log.p = log.p)

if (FALSE) {
  if (lower.tail) {
    ind1 <- (p <= pnorm.Lower)
    ans[ind1] <- Lower[ind1]
    ind2 <- (pnorm.Upper <= p)
    ans[ind2] <- Upper[ind2] 
  } else {
    ind1 <- (p >= pnorm.Lower)
    ans[ind1] <- Lower[ind1]
    ind2 <- (pnorm.Upper >= p)
    ans[ind2] <- Upper[ind2] 
  }
} else {
  ans <- qnorm(p = p, mean = mean, sd = sd,
               lower.tail = lower.tail, log.p = log.p)
  ans <- pmax(ans, Lower)
  ans <- pmin(ans, Upper)
}

  ans
}






rtobit <- function(n, mean = 0, sd = 1, Lower = 0, Upper = Inf) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n
  L <- use.n
  if (length(mean)  != L) mean  <- rep(mean,  len = L)
  if (length(sd)    != L) sd    <- rep(sd,    len = L)
  if (length(Lower) != L) Lower <- rep(Lower, len = L)
  if (length(Upper) != L) Upper <- rep(Upper, len = L)

  if (!all(Lower < Upper, na.rm = TRUE))
    stop("all(Lower < Upper) is not TRUE")

  ans <- rnorm(n = use.n, mean = mean, sd = sd)
  cenL <- (ans < Lower)
  cenU <- (ans > Upper)
  if (FALSE) {
    ans[cenL] <- Lower[cenL]
    ans[cenU] <- Upper[cenU]
  } else {
    ans <- pmax(ans, Lower)
    ans <- pmin(ans, Upper)
  }
  
  attr(ans, "Lower") <- Lower
  attr(ans, "Upper") <- Upper
  attr(ans, "cenL") <- cenL
  attr(ans, "cenU") <- cenU
  ans
}




 tobit <- function(Lower = 0, Upper = Inf,  # See the trick described below.
                   lmu = "identitylink",  lsd = "loge",
                   imu = NULL,        isd = NULL,
                   type.fitted = c("uncensored", "censored", "mean.obs"),
                   byrow.arg = FALSE,
                   imethod = 1, zero = "sd") {









  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if ( # length(Lower) != 1 || length(Upper) != 1 ||
      !is.numeric(Lower) ||
      !is.numeric(Upper) ||
      any(Lower >= Upper))
    stop("arguments 'Lower' and 'Upper' must be numeric and ",
         "satisfy Lower < Upper")


  if (mode(type.fitted) != "character" && mode(type.fitted) != "name")
    type.fitted <- as.character(substitute(type.fitted))
  type.fitted <- match.arg(type.fitted,
                           c("uncensored", "censored", "mean.obs"))[1]


  stdTobit <- all(Lower == 0.0) &&
              all(is.infinite(Upper)) &&
              all(lmu == "identitylink")


  new("vglmff",
  blurb = c("Tobit model (censored normal)\n\n",
            "Links:    ",
            namesof("mu", lmu, earg = emu, tag = TRUE), "; ",
            namesof("sd", lsd, earg = esd, tag = TRUE), "\n",
            "Mean:                 mu", "\n",
            "Conditional variance: sd^2"),
  constraints = eval(substitute(expression({

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)

  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         type.fitted = .type.fitted ,
         zero = .zero ,
         multiple.responses = TRUE,
         parameters.names = c("mu", "sd"),
         byrow.arg = .byrow.arg ,
         stdTobit = .stdTobit ,
         expected = TRUE )
  }, list( .zero = zero,
           .byrow.arg = byrow.arg,
           .stdTobit = stdTobit,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    M1 <- 2


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)
    M <- M1 * ncoly

    Lowmat <- matrix( .Lower , nrow = n, ncol = ncoly, byrow = .byrow.arg )
    Uppmat <- matrix( .Upper , nrow = n, ncol = ncoly, byrow = .byrow.arg )

    extra$type.fitted <- .type.fitted
    extra$censoredL <- (y <= Lowmat)
    extra$censoredU <- (y >= Uppmat)
    if (any(matTF <- (y < Lowmat))) {
      warning("replacing response values less than 'Lower' by 'Lower'")
      y[matTF] <- Lowmat[matTF]
    }
    if (any(matTF <- (y > Uppmat))) {
      warning("replacing response values greater than 'Upper' by 'Upper'")
      y[matTF] <- Uppmat[matTF]
    }

    temp1.names <- param.names("mu", ncoly)
    temp2.names <- param.names("sd", ncoly)
    predictors.names <-
      c(namesof(temp1.names, .lmu , earg = .emu , tag = FALSE),
        namesof(temp2.names, .lsd , earg = .esd , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {
      anyc <- cbind(extra$censoredL | extra$censoredU)
      i11 <- if ( .imethod == 1) anyc else
             matrix(FALSE, n, 1)  # can be all data

      mu.init <-
      sd.init <- matrix(0.0, n, ncoly)
      for (jay in 1:ncol(y)) {
        if ( .imethod >  2) {
          mu.init[, jay] <- (y[, jay] + weighted.mean(y[, jay], w[, jay]))/2
          sd.init[, jay] <- pmax(weighted.mean((y[, jay] - mu.init[, jay])^2,
                                                w[, jay])^0.5,
                                 0.001)
        } else {  # .imethod <= 2

          use.i11 <- i11[, jay]

          if (sum(!use.i11) < ncol(x)) {
            use.i11 <- rep(FALSE, length = n)
          }
          mylm <- lm.wfit(x = x[!use.i11,     , drop = FALSE],
                          y = y[!use.i11, jay],
                          w = w[!use.i11, jay])

                     

          sd.init[, jay] <- sqrt( sum(w[!use.i11, jay] * mylm$resid^2)
                                / mylm$df.residual ) * 1.5
          mu.init[!use.i11, jay] <- mylm$fitted.values
          if (any(anyc[, jay]))
            mu.init[anyc[, jay], jay] <- x[anyc[, jay],, drop = FALSE] %*%
                                         mylm$coeff
        }  # .imethod <= 2
      }  # for (jay in 1:ncol(y))

      if (length( .Imu ))
        mu.init <- matrix( .Imu , n, ncoly, byrow = .byrow.arg )
      if (length( .isd ))
        sd.init <- matrix( .isd , n, ncoly, byrow = .byrow.arg )

      etastart <- cbind(theta2eta(mu.init, .lmu , earg = .emu ),
                        theta2eta(sd.init, .lsd , earg = .esd ))

      etastart <- etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
    }   # if (!length(etastart))
 }), list( .Lower = Lower, .Upper = Upper,
           .lmu = lmu, .lsd = lsd,
           .emu = emu, .esd = esd,
           .Imu = imu, .isd = isd,
           .type.fitted = type.fitted,
           .stdTobit = stdTobit,
           .byrow.arg = byrow.arg,
           .imethod = imethod ))),
  linkinv = eval(substitute( function(eta, extra = NULL) {
    M1 <- 2
    ncoly <- ncol(eta) / M1
    mum <- eta2theta(eta[, M1*(1:ncoly)-1, drop = FALSE],
                     .lmu , earg = .emu )

    type.fitted <-
      if (length(extra$type.fitted)) {
        extra$type.fitted
      } else {
        warning("cannot find 'type.fitted'. Returning 'uncensored'.")
        "uncensored"
      }

    type.fitted <- match.arg(type.fitted,
                             c("uncensored", "censored", "mean.obs"))[1]

    if ( type.fitted == "uncensored")
      return(mum)

    Lowmat <- matrix( .Lower , nrow = nrow(eta), ncol = ncoly,
                      byrow = .byrow.arg )
    Uppmat <- matrix( .Upper , nrow = nrow(eta), ncol = ncoly,
                      byrow = .byrow.arg )
    if ( type.fitted == "censored") {
      mum[mum < Lowmat] <- Lowmat[mum < Lowmat]
      mum[mum > Uppmat] <- Uppmat[mum > Uppmat]
      mum
    } else {


      sdm <- eta2theta(eta[, M1*(1:ncoly)-0, drop = FALSE],
                       .lsd , earg = .esd )
      zeddL <- (Lowmat - mum) / sdm
      zeddU <- (Uppmat - mum) / sdm
      Phi.L <- pnorm(zeddL)
      phi.L <- dnorm(zeddL)
      Phi.U <- pnorm(zeddU)
      phi.U <- dnorm(zeddU)

      mum * (Phi.U - Phi.L) +
      sdm * (phi.L - phi.U) +
      ifelse(is.infinite(Lowmat), 0, Lowmat *      Phi.U ) +
      ifelse(is.infinite(Uppmat), 0, Uppmat * (1 - Phi.U))
    }
  }, list( .lmu = lmu, .lsd = lsd,
           .emu = emu, .esd = esd,
           .byrow.arg = byrow.arg,
           .Lower = Lower, .Upper = Upper ))),
  last = eval(substitute(expression({

    temp0303 <- c(rep( .lmu , length = ncoly),
                  rep( .lsd , length = ncoly))
    names(temp0303) <- c(param.names("mu", ncoly),
                         param.names("sd", ncoly))
    temp0303 <- temp0303[interleave.VGAM(M, M1 = M1)]
    misc$link <- temp0303  # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .emu
      misc$earg[[M1*ii  ]] <- .esd
    }

    misc$multipleResponses <- TRUE
    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$M1 <- M1
    misc$stdTobit <- .stdTobit
    misc$Lower <- Lowmat
    misc$Upper <- Uppmat


  }), list( .lmu = lmu, .lsd = lsd,
            .emu = emu, .esd = esd,
            .imethod = imethod,
            .stdTobit = stdTobit,
            .Lower = Lower,
            .Upper = Upper ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 2
    y <- cbind(y)
    ncoly <- ncol(y)

    cenL <- extra$censoredL
    cenU <- extra$censoredU
    cen0 <- !cenL & !cenU  # uncensored obsns
    Lowmat <- matrix( .Lower , nrow = nrow(eta), ncol = ncoly,
                      byrow = .byrow.arg )
    Uppmat <- matrix( .Upper , nrow = nrow(eta), ncol = ncoly,
                      byrow = .byrow.arg )


    mum <- eta2theta(eta[, M1*(1:ncoly)-1, drop = FALSE],
                     .lmu , earg = .emu )
    sdm <- eta2theta(eta[, M1*(1:ncoly)-0, drop = FALSE],
                     .lsd , earg = .esd )

    ell0 <- dnorm(     y[cen0], mean = mum[cen0], sd = sdm[cen0],
                  log = TRUE)
    ellL <- pnorm(Lowmat[cenL], mean = mum[cenL], sd = sdm[cenL],
                  log.p = TRUE, lower.tail = TRUE)
    ellU <- pnorm(Uppmat[cenU], mean = mum[cenU], sd = sdm[cenU],
                  log.p = TRUE, lower.tail = FALSE)

    wmat <- matrix(w, nrow = nrow(eta), ncol = ncoly)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- y  # Right dimension only
      ll.elts[cen0] <- wmat[cen0] * ell0
      ll.elts[cenL] <- wmat[cenL] * ellL
      ll.elts[cenU] <- wmat[cenU] * ellU
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmu = lmu, .lsd = lsd,
           .emu = emu, .esd = esd,
           .byrow.arg = byrow.arg,
           .Lower = Lower, .Upper = Upper ))),
  vfamily = c("tobit"),









  deriv = eval(substitute(expression({
    M1 <- 2
    y <- cbind(y)
    ncoly <- ncol(y)

    mills.ratio1 <- function(x) {
      ans <- exp(dnorm(x, log = TRUE) - pnorm(x, log = TRUE))
      ans[x < -1e2] <- -x / (1 - 1/x^2 + 3 / x^4)
      ans
    }


  mills.ratio2 <- function(x) {
    ans <- exp(2 * dnorm(x, log = TRUE) - pnorm(x, log = TRUE))
    ans[x < -40] <- 0
    ans
  }



moment.k.dnorm <- function(z, k = 0) {
  if (any(k < 0))
    stop("this function works only for non-negative 'k'")
  ans <- dnorm(z) * z^k
  ans[is.infinite(z)] <- 0
  ans
}



moment.millsratio2 <- function(zedd) {
  ans <- exp(2 * (log(abs(zedd)) + dnorm(zedd, log = TRUE)) -
             pnorm(zedd, log = TRUE))
  ans[is.infinite(zedd)] <- 0  # Needed for zedd == Inf and -Inf
  ans
}



    Lowmat <- matrix( .Lower , nrow = nrow(eta), ncol = ncoly,
                      byrow = .byrow.arg )
    Uppmat <- matrix( .Upper , nrow = nrow(eta), ncol = ncoly,
                      byrow = .byrow.arg )

    cenL <- extra$censoredL
    cenU <- extra$censoredU
    cen0 <- !cenL & !cenU  # uncensored obsns

    mum <- eta2theta(eta[, M1*(1:ncoly)-1, drop = FALSE],
                     .lmu , earg = .emu )
    sdm <- eta2theta(eta[, M1*(1:ncoly)-0, drop = FALSE],
                     .lsd , earg = .esd )

    zedd <- (y - mum) / sdm
    dl.dmu <- zedd / sdm
    dl.dsd <- (zedd^2 - 1) / sdm

    dmu.deta <- dtheta.deta(mum, .lmu , earg = .emu )
    dsd.deta <- dtheta.deta(sdm, .lsd , earg = .esd )

    if (any(cenL)) {
      mumL <- Lowmat - mum
      temp21L <- mumL[cenL] / sdm[cenL]
      fred21 <- mills.ratio1(temp21L)
      dl.dmu[cenL] <- -fred21 / sdm[cenL]
      dl.dsd[cenL] <-  fred21 * (-temp21L / sdm[cenL])
    }
    if (any(cenU)) {
      mumU <- Uppmat - mum
      temp21U <- mumU[cenU] / sdm[cenU]
      fred21 <- -mills.ratio1(-temp21U)
      dl.dmu[cenU] <- -fred21 / sdm[cenU]  # Negated
      dl.dsd[cenU] <-  fred21 * (-temp21U / sdm[cenU])
    }

    dthetas.detas <- cbind(dmu.deta, dsd.deta)
    dThetas.detas <- dthetas.detas[, interleave.VGAM(M, M1 = M1)]

    myderiv <- cbind(c(w) * dl.dmu,
                     c(w) * dl.dsd) * dthetas.detas
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lmu = lmu, .lsd = lsd,
            .emu = emu, .esd = esd,
            .byrow.arg = byrow.arg,
            .Lower = Lower, .Upper = Upper ))),
  weight = eval(substitute(expression({



    v.large <-  3.5
    v.small <- -5.0  # pnorm(-5) == 3e-07

    v.large <-  5.5
    v.small <- -6.5  # pnorm(-5) == 3e-07

    if ( .stdTobit ) {
      wz  <- matrix(0.0, n, M + M - 1)  # wz is 'tridiagonal'
      wz1 <- matrix(0.0, n, dimm(M1))
      ind1 <- iam(NA, NA, M = M1, both = TRUE, diag = TRUE)

      for (spp. in 1:ncoly) {
        zedd0 <- (            mum[, spp.]) / sdm[, spp.]
        phivec  <- dnorm(zedd0)
        Phivec  <- pnorm(zedd0)
        phicPhi <- mills.ratio1(-zedd0)

        wz1[, iam(1, 2, M = M1)] <- phivec * (1 + zedd0 *
                                    (zedd0 - phicPhi))


        wz1[, iam(1, 1, M = M1)] <- Phivec +
                                    mills.ratio2(-zedd0) +
                                    moment.k.dnorm(-zedd0, k = 1)
        wz1[, iam(2, 2, M = M1)] <- 2 * Phivec +
                                    moment.k.dnorm(-zedd0, k = 2) *
                                    mills.ratio1(-zedd0) +
                                    moment.k.dnorm(-zedd0, k = 1) +
                                    moment.k.dnorm(-zedd0, k = 3)



        if (FALSE && any(index1 <- (zedd0 < v.small))) {
          wz1[index1, iam(1, 1, M = M1)] <- 1e-7
          wz1[index1, iam(1, 2, M = M1)] <- 0
          wz1[index1, iam(2, 2, M = M1)] <- 1e-7
        }
        if (FALSE && any(index1 <- (zedd0 > v.large))) {
          wz1[index1, iam(1, 1, M = M1)] <- 1
          wz1[index1, iam(1, 2, M = M1)] <- 0
          wz1[index1, iam(2, 2, M = M1)] <- 2
        }


      wz1 <- wz1 * dThetas.detas[, M1 * (spp. - 1) + ind1$row] *
                   dThetas.detas[, M1 * (spp. - 1) + ind1$col]

      for (jay in 1:M1)
        for (kay in jay:M1) {
          cptr <- iam((spp. - 1) * M1 + jay,
                      (spp. - 1) * M1 + kay,
                      M = M)
          wz[, cptr] <- wz1[, iam(jay, kay, M = M1)]
        }
      }  # End of for (spp.) loop

    } else {  # Not a standard Tobit model ,,,,,,,,,,,,,,,,,,,,,,,,,,,,



      A.i <- (Lowmat - mum) / sdm
      B.i <- (Uppmat - mum) / sdm
      phivec.A  <- dnorm(A.i)
      phivec.B  <- dnorm(B.i)
      Phivec.A  <- pnorm(A.i)
      Phivec.B  <- pnorm(B.i)
      Phivec.BB <- pnorm(-B.i)
      phiPhi.A  <- mills.ratio1( A.i)
      phicPhi.B <- mills.ratio1(-B.i)


                         


      ned2l.dmumu <- Phivec.B - Phivec.A +
                     moment.k.dnorm( A.i, k = 1) + mills.ratio2( A.i) +
                     moment.k.dnorm(-B.i, k = 1) + mills.ratio2(-B.i)
      ned2l.dsdsd <- 2 * (Phivec.B - Phivec.A) +
                     3 * (moment.k.dnorm( A.i, k = 1) +
                          moment.k.dnorm(-B.i, k = 1)) -

                     2 * moment.k.dnorm(-B.i, k = 1) +
                     moment.k.dnorm(-B.i, k = 3) +
                     moment.millsratio2(-B.i) -
                         
                     2 * moment.k.dnorm( A.i, k = 1) +
                     moment.k.dnorm( A.i, k = 3) +
                     moment.millsratio2( A.i)
      ned2l.dmusd <- phivec.A - phivec.B +
                     moment.k.dnorm( A.i, k = 2) +
                     moment.k.dnorm( A.i, k = 1) * mills.ratio1( A.i) +
                     moment.k.dnorm( B.i, k = 2) +
                     moment.k.dnorm(-B.i, k = 1) * mills.ratio1(-B.i)



      if (TRUE && any(index1 <- (A.i < v.small))) {
        ned2l.dmusd[index1] <- 0
      }
      if (TRUE && any(index1 <- (B.i > v.large))) {
        ned2l.dmusd[index1] <- 0
      }


      wz <- array(c(ned2l.dmumu * dmu.deta^2,
                    ned2l.dsdsd * dsd.deta^2,
                    ned2l.dmusd * dmu.deta * dsd.deta),
                    dim = c(n, M / M1, 3))
      wz <- arwz2wz(wz, M = M, M1 = M1)


    }  # Not a standard Tobit model

    w.wz.merge(w = w / sdm^2, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .lmu = lmu, .Lower = Lower, .Upper = Upper,
            .lsd = lsd,
            .stdTobit = stdTobit ))))
}  # End of tobit()







 normal1 <-
 uninormal <- function(lmean = "identitylink", lsd = "loge", lvar = "loge",
                       var.arg = FALSE,
                       imethod = 1,
                       isd = NULL,
                       parallel = FALSE,
                       smallno = 1.0e-5,
                       zero = "sd") {





  apply.parint <- FALSE


  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsdev <- as.list(substitute(lsd))
  esdev <- link2list(lsdev)
  lsdev <- attr(esdev, "function.name")

  lvare <- as.list(substitute(lvar))
  evare <- link2list(lvare)
  lvare <- attr(evare, "function.name")







  if (!is.Numeric(smallno, length.arg = 1,
                  positive = TRUE))
      stop("argument 'smallno' must be positive and close to 0")
  if (smallno > 0.1) {
    warning("replacing argument 'smallno' with 0.1")
    smallno <- 0.1
  }

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
      stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.logical(var.arg) ||
      length(var.arg) != 1)
    stop("argument 'var.arg' must be a single logical")
  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")


  if (is.logical(parallel) && parallel && length(zero))
    stop("set 'zero = NULL' if 'parallel = TRUE'")


  new("vglmff",
  blurb = c("Univariate normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, earg = emean, tag = TRUE), "; ",
            if (var.arg)
            namesof("var",  lvare, earg = evare, tag = TRUE) else
            namesof("sd" ,  lsdev, earg = esdev, tag = TRUE),
            "\n",
            if (var.arg) "Variance: var" else "Variance: sd^2"),



  constraints = eval(substitute(expression({

    constraints <-
      cm.VGAM(matrix(1, M, 1), x = x,
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
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("mean", if ( .var.arg ) "var" else "sd"),
         var.arg = .var.arg ,
         parallel = .parallel ,
         zero = .zero )
  }, list( .zero = zero ,
           .parallel = parallel ,
           .var.arg = var.arg ))),

  initialize = eval(substitute(expression({
    orig.y <- y








    if (length(attr(orig.y, "Prior.Weights"))) {
      if (any(c(w) != 1))
        warning("replacing the 'weights' argument by the 'Prior.Weights'",
                "attribute of the response (probably due to Qvar()")


      w <- attr(orig.y, "Prior.Weights")


      extra$attributes.y <- attributes(orig.y)

    } else {
    }






    temp5 <-
    w.y.check(w = w, y = y,
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



    mynames1 <- param.names("mean", ncoly)
    mynames2 <- param.names(if ( .var.arg ) "var" else "sd", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lmean , earg = .emean , tag = FALSE),
          if ( .var.arg ) 
          namesof(mynames2, .lvare , earg = .evare , tag = FALSE) else
          namesof(mynames2, .lsdev , earg = .esdev , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]
    extra$predictors.names <- predictors.names


    if (!length(etastart)) {
      sdev.init <- mean.init <- matrix(0, n, ncoly)
      for (jay in 1:ncoly) {
        jfit <- lm.wfit(x = x,  y = y[, jay], w = w[, jay])
        mean.init[, jay] <- if ( .lmean == "loge")
                            pmax(1/1024, y[, jay]) else
          if ( .imethod == 1) median(y[, jay]) else
          if ( .imethod == 2) weighted.mean(y[, jay], w = w[, jay]) else
          if ( .imethod == 3) weighted.mean(y[, jay], w = w[, jay]) *
                              0.5 + y[, jay] * 0.5 else
                              mean(jfit$fitted)

        sdev.init[, jay] <-
          if ( .imethod == 1) {
            sqrt( sum(w[, jay] *
                (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else if ( .imethod == 2) {
            if (jfit$df.resid > 0)
              sqrt( sum(w[, jay] * jfit$resid^2) / jfit$df.resid ) else
              sqrt( sum(w[, jay] * jfit$resid^2) / sum(w[, jay]) )
          } else if ( .imethod == 3) {
            sqrt( sum(w[, jay] * 
                  (y[, jay] - mean.init[, jay])^2) / sum(w[, jay]) )
          } else {
            sqrt( sum(w[, jay] * abs(y[, jay] -
                                     mean.init[, jay])) / sum(w[, jay]) )
          }

        if (any(sdev.init[, jay] <= sqrt( .Machine$double.eps ) ))
          sdev.init[, jay] <- 1.01

      }


      if (length( .isdev )) {
        sdev.init <- matrix( .isdev , n, ncoly, byrow = TRUE)
      }


      etastart <-
        cbind(theta2eta(mean.init,   .lmean , earg = .emean ),
              if ( .var.arg )
              theta2eta(sdev.init^2, .lvare , earg = .evare ) else
              theta2eta(sdev.init  , .lsdev , earg = .esdev ))
      etastart <-
        etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]

      colnames(etastart) <- predictors.names
    }
  }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
            .emean = emean, .esdev = esdev, .evare = evare,
                            .isdev = isd,
            .var.arg = var.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    M1 <- extra$M1
    ncoly <- extra$ncoly



    if ( .lmean == "explink") {
      if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
        warning("turning some columns of 'eta' positive in @linkinv")
        for (ii in 1:ncoly)
          eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
      }
    }


    eta2theta(eta[, M1*(1:ncoly) - 1], .lmean , earg = .emean )
  }, list( .lmean = lmean,
           .emean = emean, .esdev = esdev , .evare = evare,
           .smallno = smallno ))),

  last = eval(substitute(expression({
    M1 <- extra$M1

    temp.names <- c(mynames1, mynames2)
    temp.names <- temp.names[interleave.VGAM(M1 * ncoly, M1 = M1)]
    misc$link <- rep( .lmean , length = M1 * ncoly)
    misc$earg <- vector("list", M1 * ncoly)
    names(misc$link) <- names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$link[ M1*ii-1 ] <- .lmean
      misc$link[ M1*ii   ] <- if ( .var.arg ) .lvare else .lsdev
      misc$earg[[M1*ii-1]] <- .emean
      misc$earg[[M1*ii  ]] <- if ( .var.arg ) .evare else .esdev
    }

    misc$var.arg <- .var.arg
    misc$M1 <- M1
    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$multipleResponses <- TRUE
    misc$parallel <- .parallel
    misc$apply.parint <- .apply.parint
    misc$smallno <- .smallno
  }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
            .emean = emean, .esdev = esdev, .evare = evare,
            .parallel = parallel, .apply.parint = apply.parint,
            .smallno = smallno,
            .var.arg = var.arg, .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    ncoly <- extra$ncoly
    M1 <- extra$M1

    if ( .lmean == "explink") {
      if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
        warning("turning some columns of 'eta' positive in @loglikelihood")
        for (ii in 1:ncoly)
          eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
      }
    }

    if ( .var.arg ) {
      Varm <- eta2theta(eta[, M1*(1:ncoly)], .lvare , earg = .evare )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, M1*(1:ncoly)], .lsdev , earg = .esdev )
    }
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnorm(y, m = mu, sd = sdev, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lsdev = lsdev, .lvare = lvare,
           .esdev = esdev, .evare = evare,
           .lmean = lmean,
           .smallno = smallno,
           .var.arg = var.arg ))),
  vfamily = c("uninormal"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    mymu <- fitted(object)
    eta <- predict(object)
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, c(FALSE, TRUE)], .lvare , earg = .evare )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, c(FALSE, TRUE)], .lsdev , earg = .esdev )
    }
    rnorm(nsim * length(mymu), mean = mymu, sd = sdev)
  }, list( .lsdev = lsdev, .lvare = lvare,
           .esdev = esdev, .evare = evare,
           .lmean = lmean,
           .smallno = smallno,
           .var.arg = var.arg ))),




  deriv = eval(substitute(expression({
    ncoly <- extra$ncoly
    M1 <- extra$M1


    if ( .lmean == "explink") {
      if (any(eta[, M1*(1:ncoly) - 1] <= 0)) {
        warning("turning some columns of 'eta' positive in @deriv")
        for (ii in 1:ncoly)
          eta[, M1*ii - 1] <- pmax( .smallno , eta[, M1*ii - 1])
      }
    }


    mymu <- eta2theta(  eta[, M1*(1:ncoly) - 1], .lmean , earg = .emean )
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, M1*(1:ncoly)    ], .lvare , earg = .evare )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, M1*(1:ncoly)    ], .lsdev , earg = .esdev )
    }

    dl.dmu <- (y - mymu) / sdev^2
    if ( .var.arg ) {
      dl.dva <- -0.5 / Varm + 0.5 * (y - mymu)^2 / sdev^4
    } else {
      dl.dsd <- -1.0 / sdev +       (y - mymu)^2 / sdev^3
    }

    dmu.deta <- dtheta.deta(mymu,   .lmean , earg = .emean )
    if ( .var.arg ) {
      dva.deta <- dtheta.deta(Varm, .lvare , earg = .evare )
    } else {
      dsd.deta <- dtheta.deta(sdev, .lsdev , earg = .esdev )
    }

    ans <- c(w) *
           cbind(dl.dmu * dmu.deta,
                 if ( .var.arg ) dl.dva * dva.deta else
                                 dl.dsd * dsd.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]






    ans
  }), list( .lmean = lmean, .lsdev = lsdev, .lvare = lvare,
            .emean = emean, .esdev = esdev, .evare = evare,
            .smallno = smallno,
            .var.arg = var.arg ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, M)  # Diagonal matrix


    ned2l.dmu2 <- 1 / sdev^2

    if ( .var.arg ) {
      ned2l.dva2 <- 0.5 / Varm^2
    } else {
      ned2l.dsd2 <- 2 / sdev^2
    }

    wz[, M1*(1:ncoly) - 1] <- ned2l.dmu2 * dmu.deta^2
    wz[, M1*(1:ncoly)    ] <- if ( .var.arg ) {
      ned2l.dva2 * dva.deta^2
    } else {
      ned2l.dsd2 * dsd.deta^2
    }

    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = ncoly)
  }), list( .var.arg = var.arg ))))
}  #  End of uninormal()









 normal.vcm <-
  function(link.list = list("(Default)" = "identitylink"),
           earg.list = list("(Default)" = list()),
           lsd = "loge", lvar = "loge",
           esd = list(), evar = list(),
           var.arg = FALSE,
           imethod = 1,
           icoefficients = NULL,
           isd = NULL,
           zero = "sd") {






  orig.esd  <- esd
  orig.evar <- evar

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")

  lvar <- as.list(substitute(lvar))
  evar <- link2list(lvar)
  lvar <- attr(evar, "function.name")






  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
      stop("argument 'imethod' must be 1 or 2 or 3 or 4")


  if (!is.logical(var.arg) || length(var.arg) != 1)
    stop("argument 'var.arg' must be a single logical")


  new("vglmff",
  blurb = c("Univariate normal distribution with ",
            "varying coefficients\n\n",
            "Links:    ",
            "G1: g1(coeff:v1), ",
            "G2: g2(coeff:v2)",
            ", ..., ",
            if (var.arg)
            namesof("var",  lvar, earg = evar, tag = TRUE) else
            namesof("sd" ,  lsd,  earg = esd,  tag = TRUE), "; ",
            "\n",
            if (var.arg) "Variance: var" else "Variance: sd^2"),

  constraints = eval(substitute(expression({



    M1 <- NA
  if (FALSE) {
    dotzero <- .zero
    if (is.character(dotzero) && dotzero == "M")
      dotzero <- M

    M1 <- NA
    eval(negzero.expression.VGAM)
  } else {
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = M)  # 20151222; Okay for one response?
  }
  }), list( .zero = zero 
          ))),

  infos = eval(substitute(function(...) {
    list(M1 = NA,
         Q1 = 1,
         multipleResponses = FALSE,  # zz unsure
         parameters.names = as.character(NA),  # zz unsure
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    asgn2 <- attr(Xm2, "assign")
    nasgn2 <- names(asgn2)




  link.list.lengths <- unlist(lapply(asgn2, length))


  link.list <- .link.list
  earg.list <- .earg.list
  if (FALSE) {
    if (length(link.list) > 0)
      if (length(nasgn2) != length(names(link.list)) ||
          !all(sort(nasgn2) == sort(names(link.list))))
        stop("names of 'link.list' do not match argument 'form2'")
    if (length(earg.list) > 0)
      if (length(nasgn2) != length(names(earg.list)) ||
          !all(sort(nasgn2) == sort(names(earg.list))))
        stop("names of 'earg.list' do not match argument 'form2'")
  }



  link.list.ordered <- vector("list", ncol(Xm2))
  earg.list.ordered <- vector("list", ncol(Xm2))




  if (sum(names(link.list) == "(Default)") > 1)
    stop("only one default allowed in argument 'link.list'!")
  if (sum(names(earg.list) == "(Default)") > 1)
    stop("only one default allowed in argument 'earg.list'!")
  default.link <- if (any(names(link.list) == "(Default)"))
    link.list[["(Default)"]] else "identitylink"
  default.earg <- if (any(names(earg.list) == "(Default)"))
    earg.list[["(Default)"]] else list()


  names(link.list.ordered) <-
  names(earg.list.ordered) <- colnames(Xm2)
  i.ptr <- 1
  for (jlocal in 1:length(nasgn2)) {
    for (klocal in 1:link.list.lengths[jlocal]) {
      link.list.ordered[[i.ptr]] <-
        if (any(names(link.list) == nasgn2[jlocal]))
          link.list[[(nasgn2[jlocal])]] else
          default.link
      earg.list.ordered[[i.ptr]] <-
        if (any(names(earg.list) == nasgn2[jlocal]))
          earg.list[[(nasgn2[jlocal])]] else
          default.earg
      i.ptr <- i.ptr + 1
    }
  }
  link.list <- link.list.ordered
  earg.list <- earg.list.ordered
  extra$link.list <- link.list
  extra$earg.list <- earg.list


    

  if (any(is.multilogit <- (unlist(link.list.ordered) == "multilogit"))) {
    if (sum(is.multilogit) < 2)
      stop("at least two 'multilogit' links need to be specified, ",
           "else none")
    col.index.is.multilogit <- (1:length(is.multilogit))[is.multilogit]
    extra$col.index.is.multilogit <- col.index.is.multilogit
    extra$is.multilogit <- is.multilogit
  }

    


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,  # M-1 ?
              out.wy = TRUE,
              colsyperw = 1,  # Use M-1, not 1, for plotvgam(y=TRUE)
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$ncoly <- ncoly <- ncol(y)
    extra$M <- M <- ncol(Xm2) + 1 -
                    (length(extra$is.multilogit) > 0)
    M1 <- NA  # Since this cannot be determined apriori.

    extra$M1 <- M1
    extra$Xm2 <- Xm2  # Needed for @linkinv
    extra$depvar <- y





  mynames1 <- paste("coeff:", colnames(Xm2), sep = "")


  for (jlocal in 1:length(mynames1)) {
    mynames1[jlocal] <- namesof(mynames1[jlocal],
                                link = link.list[[jlocal]],
                                earg = earg.list[[jlocal]], short = TRUE)
  }
  extra$all.mynames1 <- all.mynames1 <- mynames1

  if (LLL <- length(extra$is.multilogit)) {
    mynames1 <- mynames1[-max(extra$col.index.is.multilogit)]
  }

    mynames2 <- param.names(if ( .var.arg ) "var" else "sd", ncoly)
    predictors.names <-
        c(mynames1,
          if ( .var.arg ) 
          namesof(mynames2, .lvar  , earg = .evar  , tag = FALSE) else
          namesof(mynames2, .lsd   , earg = .esd   , tag = FALSE))
    extra$predictors.names <- predictors.names


    if (!length(etastart)) {

      jfit <- lm.wfit(x = Xm2,  y = c(y), w = c(w))
      jfit.coeff <- jfit$coeff




      if (icoefficients.given <- is.numeric( .icoefficients ))
        jfit.coeff <- rep( .icoefficients , length = length(jfit.coeff))



      if (!icoefficients.given)
      for (jlocal in 1:length(nasgn2)) {
        if (link.list[[jlocal]] %in%
            c("cauchit", "probit", "cloglog", "logit",
              "logc", "golf", "polf", "nbolf") &&
            abs(jfit.coeff[jlocal] - 0.5) >= 0.5)
          jfit.coeff[jlocal] <- 0.5 +
            sign(jfit.coeff[jlocal] - 0.5) * 0.25

        if (link.list[[jlocal]] %in% c("rhobit", "fisherz") &&
            abs(jfit.coeff[jlocal]) >= 1)
          jfit.coeff[jlocal] <- sign(jfit.coeff[jlocal]) * 0.5

        if (link.list[[jlocal]] == "loglog" &&
            abs(jfit.coeff[jlocal]) <= 1)
          jfit.coeff[jlocal] <- 1 + 1/8

        if (link.list[[jlocal]] == "logoff" &&
            is.numeric(LLL <- (earg.list[[jlocal]])$offset) &&
            jfit.coeff[jlocal] <= -LLL)
          jfit.coeff[jlocal] <- max((-LLL) * 1.05,
                                    (-LLL) * 0.95, -LLL + 1)

        if (link.list[[jlocal]] == "loge" &&
            jfit.coeff[jlocal] <= 0.001)
          jfit.coeff[jlocal] <- 1/8
      }

      if (!icoefficients.given)
      if (LLL <- length(extra$is.multilogit)) {
        raw.coeffs <- jfit.coeff[extra$col.index.is.multilogit]
        possum1 <- (0.01 + abs(raw.coeffs)) / sum(0.01 + abs(raw.coeffs))
        jfit.coeff[extra$is.multilogit] <- possum1
      }


      thetamat.init <- matrix(jfit.coeff, n, length(jfit.coeff),
                              byrow = TRUE)
      etamat.init <- 1 * thetamat.init  # May delete a coln later
      for (jlocal in 1:ncol(etamat.init)) {
        earg.use <- if (!length(extra$earg.list)) {
          list(theta = NULL)
        } else {
          extra$earg.list[[jlocal]]
        }

        if (length(extra$is.multilogit) && !extra$is.multilogit[jlocal])
          etamat.init[, jlocal] <-
            theta2eta(thetamat.init[, jlocal],
                      link = extra$link.list[[jlocal]],
                      earg = earg.use)
      }

      if (LLL <- length(extra$col.index.is.multilogit)) {
        etamat.init[, extra$col.index.is.multilogit[-LLL]] <-
          multilogit(thetamat.init[, extra$col.index.is.multilogit])
        etamat.init <- etamat.init[, -max(extra$col.index.is.multilogit)]
      }
      

      mean.init <- jfit$fitted
      sdev.init <-
          if ( .imethod == 1) {
            sqrt( sum(w * (y - mean.init)^2) / sum(w) )
          } else if ( .imethod == 2) {
            if (jfit$df.resid > 0)
              sqrt( sum(w * jfit$resid^2) / jfit$df.resid ) else
              sqrt( sum(w * jfit$resid^2) / sum(w) )
          } else if ( .imethod == 3) {
            sqrt( sum(w * (y - mean.init)^1.5) / sum(w) )
          } else {
            sqrt( sum(w * abs(y - mean.init)) / sum(w) )
          }

      inflation.factor <- 1.5
      sdev.init <- sdev.init * inflation.factor
      sdev.init[sdev.init <= sqrt( .Machine$double.eps )] <- 0.01

      if (length( .isdev )) {
        sdev.init <- matrix( .isdev , n, ncoly, byrow = TRUE)
      }

      etastart <-
        cbind(etamat.init,  # eta.equi.probs,
              if ( .var.arg )
              theta2eta(sdev.init^2, .lvar , earg = .evar ) else
              theta2eta(sdev.init  , .lsd  , earg = .esd  ))

      colnames(etastart) <- predictors.names
    }
  }), list( .link.list = link.list,
            .earg.list = earg.list,
            .lsd = lsd, .lvar = lvar,
            .esd = esd, .evar = evar,
            .orig.esd = orig.esd, .orig.evar = orig.evar,
            .var.arg = var.arg,
            .isdev = isd,
            .icoefficients = icoefficients,
            .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    M <- ncol(eta)

  coffs <- eta[, -M, drop = FALSE]

  if (LLL <- length(extra$col.index.is.multilogit)) {
    last.one <- extra$col.index.is.multilogit[LLL]
    coffs <- cbind(coffs[, 1:(last.one-1)],
                   probs.last.multilogit = 0,
                   if (last.one == M) NULL else
                   coffs[, last.one:ncol(coffs)])
    colnames(coffs) <- extra$all.mynames1 
  }


  for (jlocal in 1:ncol(coffs)) {
    earg.use <- if (!length(extra$earg.list[[jlocal]])) {
      list(theta = NULL)
    } else {
      extra$earg.list[[jlocal]]
    }

    if (length(extra$is.multilogit) && !extra$is.multilogit[jlocal]) {
      iskip <- (jlocal > max(extra$col.index.is.multilogit))
      coffs[, jlocal] <- eta2theta(eta[, jlocal - iskip],
                                   link = extra$link.list[[jlocal]],
                                   earg = earg.use)
    }
  }


    if (LLL <- length(extra$col.index.is.multilogit)) {
      coffs[, extra$col.index.is.multilogit] <-
        multilogit(eta[, extra$col.index.is.multilogit[-LLL], drop = FALSE],
               inverse = TRUE)
    }

    rowSums(extra$Xm2 * coffs)
  }, list( .link.list = link.list,
           .earg.list = earg.list,
           .esd = esd , .evar = evar ))),

  last = eval(substitute(expression({
    M1 <- extra$M1


    misc$link <- c(link.list.ordered,
                   "sd" = if ( .var.arg ) .lvar else .lsd )


    temp.earg.list <- c(earg.list.ordered,
                        "sd" = if ( .var.arg ) list( .orig.evar ) else
                                               list( .orig.esd  ))
    misc$earg <- temp.earg.list


    misc$var.arg <- .var.arg
    misc$M1 <- M1
    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$multipleResponses <- FALSE
    misc$icoefficients <- .icoefficients
  }), list( .link.list = link.list,
            .earg.list = earg.list,
            .lsd = lsd, .lvar = lvar,
            .esd = esd, .evar = evar,
            .orig.esd = orig.esd, .orig.evar = orig.evar,
            .icoefficients = icoefficients,
            .var.arg = var.arg, .imethod = imethod ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    if ( .var.arg ) {
      Varm <- eta2theta(eta[, ncol(eta)], .lvar , earg = .evar )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, ncol(eta)], .lsd  , earg = .esd  )
    }
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dnorm(y, m = mu, sd = sdev, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lsd = lsd, .lvar = lvar,
           .esd = esd, .evar = evar,
           .var.arg = var.arg ))),
  vfamily = c("normal.vcm"),



  deriv = eval(substitute(expression({

    if ( .var.arg ) {
      Varm <- eta2theta(eta[, M], .lvar , earg = .evar )
      sdev <- sqrt(Varm)
    } else {
      sdev <- eta2theta(eta[, M], .lsd  , earg = .esd  )
    }

    zedd <- (y - mu) / sdev
    dl.dmu <- c(zedd / sdev)  #   dl.dmu <- (y - mymu) / sdev^2

    dmu.dcoffs <- Xm2

    mymu <- mu


    coffs <- eta[, -M, drop = FALSE]  # Exclude log(sdev) or log(var)

    if (LLL <- length(extra$is.multilogit)) {
      last.one <- max(extra$col.index.is.multilogit)
      coffs <- cbind(coffs[, 1:(last.one-1)],
                     probsLastmultilogit = 0,
                     if (last.one == M) NULL else
                     coffs[, last.one:ncol(coffs)])
      colnames(coffs) <- extra$all.mynames1
    }

    dcoffs.deta <- coffs  # Includes any last "multilogit"

    for (jlocal in 1:ncol(coffs)) {
      earg.use <- if (!length(extra$earg.list[[jlocal]])) {
        list(theta = NULL)
      } else {
        extra$earg.list[[jlocal]]
      }

      if (!length(extra$is.multilogit) ||
          !extra$is.multilogit[jlocal]) {
        iskip <- length(extra$is.multilogit) &&
                 (jlocal  > max(extra$col.index.is.multilogit))
        coffs[, jlocal] <- eta2theta(eta[, jlocal - iskip],
                                     link = extra$link.list[[jlocal]],
                                     earg = earg.use)
      }
    }

    if (LLL <- length(extra$col.index.is.multilogit)) {
      coffs[, extra$col.index.is.multilogit] <-
        multilogit(eta[, extra$col.index.is.multilogit[-LLL], drop = FALSE],
               inverse = TRUE)
    }


  for (jlocal in 1:ncol(coffs)) {
    if (!length(extra$is.multilogit) ||
        !extra$is.multilogit[jlocal]) {
      earg.use <- if (!length(extra$earg.list[[jlocal]])) {
        list(theta = NULL)
      } else {
        extra$earg.list[[jlocal]]
      }
      dcoffs.deta[, jlocal] <-
        dtheta.deta(coffs[, jlocal], 
                    link = extra$link.list[[jlocal]],
                    earg = earg.use)
    }
  }



    if ( .var.arg ) {
      dl.dva <- -0.5 / Varm + 0.5 * (y - mymu)^2 / sdev^4
    } else {
      dl.dsd <- -1.0 / sdev +       (y - mymu)^2 / sdev^3
    }

    if ( .var.arg ) {
      dva.deta <- dtheta.deta(Varm, .lvar , earg = .evar )
    } else {
      dsd.deta <- dtheta.deta(sdev, .lsd  , earg = .esd )
    }

    
    dMu.deta <- dmu.dcoffs * dcoffs.deta  # n x pLM, but may change below
    if (LLL <- length(extra$col.index.is.multilogit)) {
      dMu.deta[, extra$col.index.is.multilogit[-LLL]] <-
         coffs[, extra$col.index.is.multilogit[-LLL]] *
        (dmu.dcoffs[, extra$col.index.is.multilogit[-LLL]] -
         rowSums(dmu.dcoffs[, extra$col.index.is.multilogit]  *
                      coffs[, extra$col.index.is.multilogit]))
      dMu.deta <- dMu.deta[, -extra$col.index.is.multilogit[LLL]]
    }
    

    dl.deta <- if ( .var.arg )
               c(w) * cbind(dl.dmu * dMu.deta,
                            "var" = c(dl.dva * dva.deta)) else
               c(w) * cbind(dl.dmu * dMu.deta,
                            "sd"  = c(dl.dsd * dsd.deta))
 
    dl.deta
  }), list( .link.list = link.list, .lsd = lsd, .lvar = lvar,
            .earg.list = earg.list, .esd = esd, .evar = evar,
            .var.arg = var.arg ))),

      



  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, dimm(M))  # Treated as a general full matrix


    wz[, iam(M, M, M = M)] <- if ( .var.arg ) {
      ned2l.dva2 <- 0.5 / Varm^2
      ned2l.dva2 * dva.deta^2
    } else {
      ned2l.dsd2 <- 2 / sdev^2
      ned2l.dsd2 * dsd.deta^2
    }



    if (length(extra$col.index.is.multilogit)) {
      LLL <- max(extra$col.index.is.multilogit)
      dmu.dcoffs <- dmu.dcoffs[, -LLL]
      dcoffs.deta <- dcoffs.deta[, -LLL]
    }


    index <- iam(NA, NA, M  , both = TRUE, diag = TRUE)
    indtw <- iam(NA, NA, M-1, both = TRUE, diag = TRUE)
    ned2l.dmu2 <- 1 / sdev^2

 



 
    if ((LLL <- length(extra$col.index.is.multilogit))) {
       dmu.dcoffs[, extra$col.index.is.multilogit[-LLL]] <-
         dMu.deta[, extra$col.index.is.multilogit[-LLL]]
      dcoffs.deta[, extra$col.index.is.multilogit[-LLL]] <- 1
     }
  
    twz  <- crossprod(dmu.dcoffs * sqrt(c(w))) / sum(w)

    twz <- matrix(twz[cbind(indtw$row.index,
                            indtw$col.index)],
                  n, dimm(M-1), byrow = TRUE)
    if (length(indtw$row.index) != dimm(M-1))
      stop("dim of twz incorrect")

    twz <- twz *
           dcoffs.deta[, indtw$row.index, drop = FALSE] *
           dcoffs.deta[, indtw$col.index, drop = FALSE] *
           ned2l.dmu2

    for (ilocal in 1:length(indtw$row.index))
      wz[, iam(indtw$row.index[ilocal],
               indtw$col.index[ilocal], M = M)] <-
     twz[, iam(indtw$row.index[ilocal],
               indtw$col.index[ilocal], M = M-1)]


    c(w) * wz
  }), list( .var.arg = var.arg ))))
}  # End of normal.vcm()








 lognormal <- function(lmeanlog = "identitylink", lsdlog = "loge",
                       zero = "sdlog") {




  lmulog <- as.list(substitute(lmeanlog))
  emulog <- link2list(lmulog)
  lmulog <- attr(emulog, "function.name")

  lsdlog <- as.list(substitute(lsdlog))
  esdlog <- link2list(lsdlog)
  lsdlog <- attr(esdlog, "function.name")






  new("vglmff",
  blurb = c("Two-parameter (univariate) lognormal distribution\n\n",
          "Links:    ",
          namesof("meanlog", lmulog, earg = emulog, tag = TRUE), ", ",
          namesof("sdlog",   lsdlog, earg = esdlog, tag = TRUE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         lmeanlog = .lmeanlog ,
         lsdlog   = .lsdlog ,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("meanlog", "sdlog"),
         zero = .zero )
  }, list( .zero = zero,
           .lmeanlog = lmeanlog,
           .lsdlog   = lsdlog
         ))),


  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE)



    predictors.names <-
        c(namesof("meanlog", .lmulog , earg = .emulog , tag = FALSE),
          namesof("sdlog",   .lsdlog , earg = .esdlog , tag = FALSE))

    if (!length(etastart)) {
      mylm <- lm.wfit(x = x, y = c(log(y)), w = c(w))
      sdlog.y.est <- sqrt( sum(c(w) * mylm$resid^2) / mylm$df.residual )
      etastart <- cbind(
        meanlog = rep(theta2eta(log(median(y)), .lmulog ,
                                earg = .emulog ), length = n),
        sdlog   = rep(theta2eta(sdlog.y.est, .lsdlog ,
                                earg = .esdlog ), length = n))
    }
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
    sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
    exp(mulog + 0.5 * sdlog^2)
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),
  last = eval(substitute(expression({
    misc$link <-    c("meanlog" = .lmulog , "sdlog" = .lsdlog )

    misc$earg <- list("meanlog" = .emulog , "sdlog" = .esdlog )

    misc$expected <- TRUE
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
    sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlnorm(y, meanlog = mulog, sdlog = sdlog,
                               log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),
  vfamily = c("lognormal"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    mulog <- eta2theta(eta[, c(TRUE, FALSE)], .lmulog , earg = .emulog )
    sdlog <- eta2theta(eta[, c(FALSE, TRUE)], .lsdlog , earg = .esdlog )
    rlnorm(nsim * length(mulog),
           meanlog = mulog, sdlog = sdlog)
  }, list( .lmulog = lmulog, .lsdlog = lsdlog,
           .emulog = emulog, .esdlog = esdlog ))),



  deriv = eval(substitute(expression({
    mulog <- eta2theta(eta[, 1], .lmulog , earg = .emulog )
    sdlog <- eta2theta(eta[, 2], .lsdlog , earg = .esdlog )

    dmulog.deta <- dtheta.deta(mulog, .lmulog , earg = .emulog )
    dsdlog.deta <- dtheta.deta(sdlog, .lsdlog ,   earg = .esdlog )

    dl.dmulog <- (log(y) - mulog) / sdlog^2
    dl.dsdlog <- -1 / sdlog + (log(y) - mulog)^2 / sdlog^3

    c(w) * cbind(dl.dmulog * dmulog.deta, 
                 dl.dsdlog * dsdlog.deta)
  }), list( .lmulog = lmulog, .lsdlog = lsdlog,
            .emulog = emulog, .esdlog = esdlog ))),
  weight = expression({
    wz <- matrix(NA_real_, n, 2)  # Diagonal!
    ned2l.dmulog2 <- 1 / sdlog^2
    ned2l.dsdlog2 <- 2 * ned2l.dmulog2

    wz[, iam(1, 1, M)] <- ned2l.dmulog2 * dmulog.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dsdlog2 * dsdlog.deta^2

    wz = c(w) * wz
    wz
  }))
}






dskewnorm <- function(x, location = 0, scale = 1, shape = 0, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  zedd <- (x - location) / scale
  loglik <- log(2) + dnorm(zedd, log = TRUE) +
            pnorm(shape * zedd, log.p = TRUE) - log(scale)

  loglik[is.infinite(x)] <- log(0)  # 20141209 KaiH

  if (log.arg) {
    loglik
  } else {
    exp(loglik)
  }
}



rskewnorm <- function(n, location = 0, scale = 1, shape = 0) {


  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n

  rho <- shape / sqrt(1 + shape^2)
  u0 <- rnorm(use.n)
  v  <- rnorm(use.n)
  u1 <- rho * u0 + sqrt(1 - rho^2) * v




  ans <- location + scale * sign(u0) * u1

  ans[scale <= 0] <- NA
  ans
}






 skewnormal <- function(lshape = "identitylink",
                        ishape = NULL,
                        nsimEIM = NULL) {


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (length(nsimEIM) &&
     (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10))
    stop("argument 'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("1-parameter skew-normal distribution\n\n",
          "Link:     ",
          namesof("shape", lshape , earg = eshape), "\n",
          "Mean:     shape * sqrt(2 / (pi * (1 + shape^2 )))\n",
          "Variance: 1-mu^2"),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("shape"),
         nsimEIM = .nsimEIM)
  }, list( .nsimEIM = nsimEIM ))),
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
      namesof("shape", .lshape , earg = .eshape , tag = FALSE)

    if (!length(etastart)) {
      init.shape <- if (length( .ishape ))
        rep( .ishape , len = n) else {
        temp <- y
        index <- abs(y) < sqrt(2/pi)-0.01
        temp[!index] <- y[!index]
        temp[index] <- sign(y[index]) / sqrt(2/(pi*y[index]*y[index])-1)
        temp
      }
      etastart <- matrix(init.shape, n, ncol(y))
    }
  }), list( .lshape = lshape, .eshape = eshape,
            .ishape = ishape ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    alpha <- eta2theta(eta, .lshape , earg = .eshape )
    alpha * sqrt(2/(pi * (1+alpha^2 )))
  }, list( .eshape = eshape, .lshape = lshape ))),
  last = eval(substitute(expression({



    misc$link <-    c(shape = .lshape) 

    misc$earg <- list(shape = .eshape )

    misc$nsimEIM <- .nsimEIM
      misc$expected <- (length( .nsimEIM ) > 0)
  }), list( .eshape = eshape, .lshape = lshape,
            .nsimEIM = nsimEIM ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    alpha <- mu / sqrt(2/pi - mu^2)
    theta2eta(alpha, .lshape , earg = .eshape )
  }, list( .eshape = eshape, .lshape = lshape ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
      alpha <- eta2theta(eta, .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dskewnorm(x = y, location = 0, scale = 1,
                                  shape = alpha, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .eshape = eshape, .lshape = lshape ))), 
  vfamily = c("skewnormal"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    alpha <- eta2theta(eta, .lshape , earg = .eshape )
    rskewnorm(nsim * length(alpha), location = 0, scale = 1,
              shape = alpha)
  }, list( .eshape = eshape, .lshape = lshape ))), 






  deriv = eval(substitute(expression({
    alpha <- eta2theta(eta, .lshape , earg = .eshape )

    zedd <- y*alpha
    tmp76 <- pnorm(zedd)
    tmp86 <- dnorm(zedd)
    dl.dshape <- tmp86 * y / tmp76

    dshape.deta <- dtheta.deta(alpha, .lshape , earg = .eshape )

    c(w) * dl.dshape * dshape.deta
  }), list( .eshape = eshape,
            .lshape = lshape ))),
  weight = eval(substitute(expression({
    if ( length( .nsimEIM )) {
      run.mean <- 0
      for (ii in 1:( .nsimEIM)) {
        ysim <- rsnorm(n, location = 0, scale = 1, shape = alpha)
        zedd <- ysim*alpha
        tmp76 <- pnorm(zedd)
        tmp86 <- dnorm(zedd)
        d2l.dshape2 <- -ysim*ysim*tmp86*(tmp76*zedd+tmp86)/tmp76^2
        rm(ysim)
        run.mean <- ((ii-1) * run.mean + d2l.dshape2) / ii
      }
      if (intercept.only)
        run.mean <- mean(run.mean)
      wz <-  -c(w) * (dshape.deta^2) * run.mean
    } else {
      d2shape.deta2 <- d2theta.deta2(alpha, .lshape , earg = .eshape )
      d2l.dshape2 <- -y*y * tmp86 * (tmp76 * zedd + tmp86) / tmp76^2
      wz <- -(dshape.deta^2) * d2l.dshape2 - d2shape.deta2 * dl.dshape
      wz <- c(w) * wz
    }
    wz
  }), list( .eshape = eshape,
            .lshape = lshape, .nsimEIM = nsimEIM ))))
}






