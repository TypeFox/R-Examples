# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












mix2normal.control <- function(trace = TRUE, ...) {
  list(trace = trace)
}


 mix2normal <-
    function(lphi = "logit",
             lmu = "identitylink",
             lsd = "loge",
             iphi = 0.5,
             imu1 = NULL, imu2 = NULL,
             isd1 = NULL, isd2 = NULL,
             qmu = c(0.2, 0.8),
             eq.sd = TRUE,
             nsimEIM = 100,
             zero = "phi") {
  lphi <- as.list(substitute(lphi))
  ephi <- link2list(lphi)
  lphi <- attr(ephi, "function.name")

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")


  emu1 <- emu2 <- emu
  esd1 <- esd2 <- esd


  if (!is.Numeric(qmu, length.arg = 2,
                  positive = TRUE) ||
      any(qmu >= 1))
    stop("bad input for argument 'qmu'")


  if (length(iphi) &&
     (!is.Numeric(iphi, length.arg = 1,
                  positive = TRUE) ||
      iphi>= 1))
      stop("bad input for argument 'iphi'")
  if (length(imu1) && !is.Numeric(imu1))
    stop("bad input for argument 'imu1'")
  if (length(imu2) && !is.Numeric(imu2))
    stop("bad input for argument 'imu2'")
  if (length(isd1) && !is.Numeric(isd1, positive = TRUE))
    stop("bad input for argument 'isd1'")
  if (length(isd2) && !is.Numeric(isd2, positive = TRUE))
    stop("bad input for argument 'isd2'")


  if (!is.logical(eq.sd) || length(eq.sd) != 1)
    stop("bad input for argument 'eq.sd'")
  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("Mixture of two univariate normals\n\n",
            "Links:    ",
            namesof("phi", lphi, earg = ephi, tag = FALSE), ", ", 
            namesof("mu1",  lmu, earg = emu1, tag = FALSE), ", ",
            namesof("sd1",  lsd, earg = esd1, tag = FALSE), ", ",
            namesof("mu2",  lmu, earg = emu2, tag = FALSE), ", ",
            namesof("sd2",  lsd, earg = esd2, tag = FALSE), "\n",
            "Mean:     phi*mu1 + (1 - phi)*mu2\n",
            "Variance: phi*sd1^2 + (1 - phi)*sd2^2 + ",
                      "phi*(1 - phi)*(mu1-mu2)^2"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(rbind(diag(4), c(0, 0, 1, 0)), x = x,
                           bool = .eq.sd ,
                           constraints = constraints,
                           apply.int = TRUE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 5)
  }), list( .zero = zero, .eq.sd = eq.sd ))),


  infos = eval(substitute(function(...) {
    list(M1 = 5,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("phi", "mu1", "sd1", "mu2", "sd2"),
         nsimEIM = .nsimEIM ,
         lphi      = .lphi   ,
         lmu1      = .lmu    ,
         lsd1      = .lsd    ,
         lmu2      = .lmu    ,
         lsd2      = .lsd    ,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM,
           .lphi = lphi,
           .lmu  = lmu , .lsd = lsd
         ))),


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
        namesof("phi", .lphi , earg = .ephi , tag = FALSE),
        namesof("mu1", .lmu  , earg = .emu1 , tag = FALSE),
        namesof("sd1", .lsd  , earg = .esd1 , tag = FALSE),
        namesof("mu2", .lmu  , earg = .emu2 , tag = FALSE),
        namesof("sd2", .lsd  , earg = .esd2 , tag = FALSE))



    if (!length(etastart)) {
      qy <- quantile(y, prob = .qmu )
      init.phi <- rep(if (length( .iphi )) .iphi else   0.5, length = n)
      init.mu1 <- rep(if (length( .imu1 )) .imu1 else qy[1], length = n)
      init.mu2 <- rep(if (length( .imu2 )) .imu2 else qy[2], length = n)
      ind.1 <- if (init.mu1[1] < init.mu2[1])
                1:round(n* init.phi[1]) else
                round(n* init.phi[1]):n
      ind.2 <- if (init.mu1[1] < init.mu2[1])
                round(n* init.phi[1]):n else
                1:round(n* init.phi[1])
      sorty <- sort(y)
      init.sd1 <- rep(if (length( .isd1 )) .isd1 else sd(sorty[ind.1]),
                      len = n)
      init.sd2 <- rep(if (length( .isd2 )) .isd2 else sd(sorty[ind.2]),
                      len = n)
      if ( .eq.sd ) {
        init.sd1 <-
        init.sd2 <- (init.sd1 + init.sd2) / 2
        if (!all.equal( .esd1, .esd2 ))
          stop("'esd1' and 'esd2' must be equal if 'eq.sd = TRUE'")
      }
      etastart <- cbind(
                  theta2eta(init.phi, .lphi , earg = .ephi ),
                  theta2eta(init.mu1,  .lmu , earg = .emu1 ),
                  theta2eta(init.sd1,  .lsd , earg = .esd1 ),
                  theta2eta(init.mu2,  .lmu , earg = .emu2 ),
                  theta2eta(init.sd2,  .lsd , earg = .esd2 ))
    }
  }), list(.lphi = lphi, .lmu = lmu,
           .iphi = iphi, .imu1 = imu1, .imu2 = imu2,
           .ephi = ephi, .emu1 = emu1, .emu2 = emu2,
           .esd1 = esd1, .esd2 = esd2, .eq.sd = eq.sd,
           .lsd = lsd, .isd1 = isd1, .isd2 = isd2, .qmu = qmu))),
  linkinv = eval(substitute(function(eta, extra = NULL){
      phi <- eta2theta(eta[, 1], link = .lphi , earg = .ephi )
      mu1 <- eta2theta(eta[, 2], link =  .lmu , earg = .emu1 )
      mu2 <- eta2theta(eta[, 4], link =  .lmu , earg = .emu2 )
      phi * mu1 + (1 - phi) * mu2
  }, list( .lphi = lphi, .lmu = lmu,
           .ephi = ephi, .emu1 = emu1, .emu2 = emu2,
           .esd1 = esd1, .esd2 = esd2 ))),
  last = eval(substitute(expression({
    misc$link <-    c("phi" = .lphi , "mu1" = .lmu ,
                      "sd1" = .lsd  , "mu2" = .lmu , "sd2" = .lsd )

    misc$earg <- list("phi" = .ephi , "mu1" = .emu1 ,
                      "sd1" = .esd1 , "mu2" = .emu2 , "sd2" = .esd2 )

    misc$expected <- TRUE
    misc$eq.sd <- .eq.sd
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list(.lphi = lphi, .lmu = lmu, .lsd = lsd, .eq.sd = eq.sd,
           .ephi = ephi, .emu1 = emu1, .emu2 = emu2,
           .esd1 = esd1, .esd2 = esd2,
           .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    phi <- eta2theta(eta[, 1], link = .lphi , earg = .ephi )
    mu1 <- eta2theta(eta[, 2], link = .lmu  , earg = .emu1 )
    sd1 <- eta2theta(eta[, 3], link = .lsd  , earg = .esd1 )
    mu2 <- eta2theta(eta[, 4], link = .lmu  , earg = .emu2 )
    sd2 <- eta2theta(eta[, 5], link = .lsd  , earg = .esd2 )
    f1 <- dnorm(y, mean = mu1, sd = sd1)
    f2 <- dnorm(y, mean = mu2, sd = sd2)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * log(phi*f1 + (1 - phi)*f2)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(.lphi = lphi, .lmu = lmu,
          .ephi = ephi, .emu1 = emu1, .emu2 = emu2,
          .esd1 = esd1, .esd2 = esd2,
          .lsd = lsd ))),
  vfamily = c("mix2normal"),
  deriv = eval(substitute(expression({
    phi <- eta2theta(eta[, 1], link = .lphi , earg = .ephi )
    mu1 <- eta2theta(eta[, 2], link = .lmu  , earg = .emu1 )
    sd1 <- eta2theta(eta[, 3], link = .lsd  , earg = .esd1 )
    mu2 <- eta2theta(eta[, 4], link = .lmu  , earg = .emu2 )
    sd2 <- eta2theta(eta[, 5], link = .lsd  , earg = .esd2 )
    dphi.deta <- dtheta.deta(phi, link = .lphi , earg = .ephi )
    dmu1.deta <- dtheta.deta(mu1, link = .lmu  , earg = .emu1 )
    dmu2.deta <- dtheta.deta(mu2, link = .lmu  , earg = .emu2 )
    dsd1.deta <- dtheta.deta(sd1, link = .lsd  , earg = .esd1 )
    dsd2.deta <- dtheta.deta(sd2, link = .lsd  , earg = .esd2 )
    f1 <- dnorm(y, mean = mu1, sd = sd1)
    f2 <- dnorm(y, mean = mu2, sd = sd2)
    pdf <- phi*f1 + (1 - phi)*f2
    z1 <- (y-mu1) / sd1
    z2 <- (y-mu2) / sd2
    df1.dmu1 <- z1 * f1 / sd1
    df2.dmu2 <- z2 * f2 / sd2
    df1.dsd1 <- (z1^2 - 1) * f1 / sd1
    df2.dsd2 <- (z2^2 - 1) * f2 / sd2
    dl.dphi <- (f1-f2) / pdf
    dl.dmu1 <- phi * df1.dmu1 / pdf
    dl.dmu2 <- (1 - phi) * df2.dmu2 / pdf
    dl.dsd1 <- phi * df1.dsd1 / pdf
    dl.dsd2 <- (1 - phi) * df2.dsd2 / pdf
    c(w) * cbind(dl.dphi * dphi.deta,
                 dl.dmu1 * dmu1.deta,
                 dl.dsd1 * dsd1.deta,
                 dl.dmu2 * dmu2.deta,
                 dl.dsd2 * dsd2.deta)
  }), list(.lphi = lphi, .lmu = lmu, .lsd = lsd,
           .ephi = ephi, .emu1 = emu1, .emu2 = emu2,
           .esd1 = esd1, .esd2 = esd2,
           .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({

    d3 <- deriv3(~ log(     phi  * dnorm((ysim-mu1)/sd1) / sd1 +
                       (1 - phi) * dnorm((ysim-mu2)/sd2) / sd2),
        c("phi","mu1","sd1","mu2","sd2"), hessian = TRUE)
    run.mean <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- ifelse(runif(n) < phi, rnorm(n, mu1, sd1),
                                     rnorm(n, mu2, sd2))

        eval.d3 <- eval(d3)
      d2l.dthetas2 <-  attr(eval.d3, "hessian")
      rm(ysim)

      temp3 <- matrix(0, n, dimm(M))
      for (ss in 1:M)
        for (tt in ss:M)
          temp3[,iam(ss,tt, M)] <-  -d2l.dthetas2[, ss, tt]

      run.mean <- ((ii-1) * run.mean + temp3) / ii
    }
    wz <- if (intercept.only)
      matrix(colMeans(run.mean), n, dimm(M), byrow = TRUE) else
      run.mean

    dtheta.detas <- cbind(dphi.deta,
                          dmu1.deta,
                          dsd1.deta,
                          dmu2.deta,
                          dsd2.deta)
    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list(.lphi = lphi, .lmu = lmu, .nsimEIM = nsimEIM ))))
}




mix2poisson.control <- function(trace = TRUE, ...) {
    list(trace = trace)
}


 mix2poisson <- function(lphi = "logit", llambda = "loge",
                         iphi = 0.5, il1 = NULL, il2 = NULL,
                         qmu = c(0.2, 0.8), nsimEIM = 100,
                         zero = "phi") {

  lphi <- as.list(substitute(lphi))
  ephi <- link2list(lphi)
  lphi <- attr(ephi, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  el1 <- el2 <- elambda



  if (!is.Numeric(qmu, length.arg = 2, positive = TRUE) ||
      any(qmu >= 1))
    stop("bad input for argument 'qmu'")
  if (length(iphi) &&
     (!is.Numeric(iphi, length.arg = 1, positive = TRUE) ||
     iphi >= 1))
    stop("bad input for argument 'iphi'")
  if (length(il1) && !is.Numeric(il1))
    stop("bad input for argument 'il1'")
  if (length(il2) && !is.Numeric(il2))
    stop("bad input for argument 'il2'")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("Mixture of two Poisson distributions\n\n",
            "Links:    ",
            namesof("phi",lphi, earg = ephi), ", ", 
            namesof("lambda1", llambda, earg = el1, tag = FALSE), ", ",
            namesof("lambda2", llambda, earg = el2, tag = FALSE), "\n",
            "Mean:     phi*lambda1 + (1 - phi)*lambda2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("phi", "lambda1", "lambda2"),
         nsimEIM = .nsimEIM ,
         lphi      = .lphi   ,
         llambda1      = .llambda    ,
         llambda2      = .llambda    ,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM,
           .lphi = lphi,
           .llambda = llambda
         ))),


  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              Is.integer.y = TRUE,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    predictors.names <-
      c(namesof("phi",     .lphi    , earg = .ephi , tag = FALSE),
        namesof("lambda1", .llambda , earg = .el1  , tag = FALSE),
        namesof("lambda2", .llambda , earg = .el2  , tag = FALSE))

    if (!length(etastart)) {
      qy <- quantile(y, prob =  .qmu)
      init.phi <-     rep(if (length( .iphi )) .iphi else 0.5,   length = n)
      init.lambda1 <- rep(if (length( .il1  )) .il1  else qy[1], length = n)
      init.lambda2 <- rep(if (length( .il2  )) .il2  else qy[2], length = n)

      if (!length(etastart))  
        etastart <- cbind(theta2eta(init.phi, .lphi , earg = .ephi ),
                          theta2eta(init.lambda1, .llambda , earg = .el1 ),
                          theta2eta(init.lambda2, .llambda , earg = .el2 ))
    }
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2,
           .iphi = iphi, .il1 = il1, .il2 = il2,
           .qmu = qmu))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    phi     <- eta2theta(eta[, 1], link = .lphi ,    earg = .ephi )
    lambda1 <- eta2theta(eta[, 2], link = .llambda , earg = .el1  )
    lambda2 <- eta2theta(eta[, 3], link = .llambda , earg = .el2  )
    phi * lambda1 + (1 - phi) * lambda2
  }, list(.lphi = lphi, .llambda = llambda,
          .ephi = ephi, .el1 = el1, .el2 = el2 ))),
  last = eval(substitute(expression({
    misc$link <-
         c("phi" = .lphi , "lambda1" = .llambda , "lambda2" = .llambda )

    misc$earg <-
      list("phi" = .ephi , "lambda1" = .el1 ,     "lambda2" = .el2 )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2,
           .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    phi     <- eta2theta(eta[, 1], link = .lphi    , earg = .ephi )
    lambda1 <- eta2theta(eta[, 2], link = .llambda , earg = .el1  )
    lambda2 <- eta2theta(eta[, 3], link = .llambda , earg = .el2  )
    f1 <- dpois(y, lam = lambda1)
    f2 <- dpois(y, lam = lambda2)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * log(phi*f1 + (1 - phi)*f2)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2 ))),
  vfamily = c("mix2poisson"),
  deriv = eval(substitute(expression({
    phi     <- eta2theta(eta[, 1], link = .lphi    , earg = .ephi )
    lambda1 <- eta2theta(eta[, 2], link = .llambda , earg = .el1  )
    lambda2 <- eta2theta(eta[, 3], link = .llambda , earg = .el2  )

    dphi.deta     <- dtheta.deta(phi,     link = .lphi    , earg = .ephi )
    dlambda1.deta <- dtheta.deta(lambda1, link = .llambda , earg = .el1  )
    dlambda2.deta <- dtheta.deta(lambda2, link = .llambda , earg = .el2  )

    f1 <- dpois(x = y, lam = lambda1)
    f2 <- dpois(x = y, lam = lambda2)
    pdf <- phi*f1 + (1 - phi)*f2
    df1.dlambda1 <- dpois(y-1, lam = lambda1) - f1
    df2.dlambda2 <- dpois(y-1, lam = lambda2) - f2
    dl.dphi <- (f1-f2) / pdf
    dl.dlambda1 <- phi * df1.dlambda1 / pdf
    dl.dlambda2 <- (1 - phi) * df2.dlambda2 / pdf

    c(w) * cbind(dl.dphi * dphi.deta,
                 dl.dlambda1 * dlambda1.deta,
                 dl.dlambda2 * dlambda2.deta)
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2,
           .nsimEIM = nsimEIM ))),
  weight = eval(substitute(expression({
    run.mean <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- ifelse(runif(n) < phi, rpois(n, lambda1),
                                     rpois(n, lambda2))
      f1 <- dpois(x = ysim, lam = lambda1)
      f2 <- dpois(x = ysim, lam = lambda2)
      pdf <- phi*f1 + (1 - phi)*f2

      df1.dlambda1 <- dpois(ysim-1, lam = lambda1) - f1
      df2.dlambda2 <- dpois(ysim-1, lam = lambda2) - f2

      dl.dphi <- (f1 - f2) / pdf
      dl.dlambda1 <- phi * df1.dlambda1 / pdf
      dl.dlambda2 <- (1 - phi) * df2.dlambda2 / pdf

      d2f1.dlambda12 <- dpois(ysim-2, lambda1) -
                     2*dpois(ysim-1, lambda1) +
                       dpois(ysim, lambda1)
      d2f2.dlambda22 <- dpois(ysim-2, lambda2) -
                     2*dpois(ysim-1, lambda2) +
                       dpois(ysim, lambda2)
      d2l.dphi2 <-  dl.dphi^2
      d2l.dlambda12 <- phi * (phi * df1.dlambda1^2 / pdf -
                       d2f1.dlambda12) / pdf
      d2l.dlambda22 <- (1 - phi) * ((1 - phi) * df2.dlambda2^2 / pdf -
                       d2f2.dlambda22) / pdf
      d2l.dlambda1lambda2 <-  phi * (1 - phi) *
                              df1.dlambda1 * df2.dlambda2 / pdf^2
      d2l.dphilambda1 <- df1.dlambda1 * (phi*(f1-f2)/pdf - 1) / pdf
      d2l.dphilambda2 <- df2.dlambda2 * ((1 - phi)*(f1-f2)/pdf - 1) / pdf

      rm(ysim)
      temp3 <- matrix(0, n, dimm(M))
      temp3[, iam(1, 1, M = 3)] <- d2l.dphi2
      temp3[, iam(2, 2, M = 3)] <- d2l.dlambda12
      temp3[, iam(3, 3, M = 3)] <- d2l.dlambda22
      temp3[, iam(1, 2, M = 3)] <- d2l.dphilambda1
      temp3[, iam(1, 3, M = 3)] <- d2l.dphilambda2
      temp3[, iam(2, 3, M = 3)] <- d2l.dlambda1lambda2
      run.mean <- ((ii-1) * run.mean + temp3) / ii
    }

    wz <- if (intercept.only)
          matrix(colMeans(run.mean), n, dimm(M), byrow = TRUE) else
          run.mean

    dtheta.detas <- cbind(dphi.deta, dlambda1.deta, dlambda2.deta)
    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]

    c(w) * wz
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2,
           .nsimEIM = nsimEIM ))))
}





mix2exp.control <- function(trace = TRUE, ...) {
  list(trace = trace)
}



 mix2exp <- function(lphi = "logit", llambda = "loge",
                     iphi = 0.5, il1 = NULL, il2 = NULL,
                     qmu = c(0.8, 0.2), nsimEIM = 100,
                     zero = "phi") {
  lphi <- as.list(substitute(lphi))
  ephi <- link2list(lphi)
  lphi <- attr(ephi, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  el1 <- el2 <- elambda


  if (!is.Numeric(qmu, length.arg = 2, positive = TRUE) ||
      any(qmu >= 1))
    stop("bad input for argument 'qmu'")
  if (length(iphi) &&
     (!is.Numeric(iphi, length.arg = 1, positive = TRUE) ||
      iphi >= 1))
    stop("bad input for argument 'iphi'")
  if (length(il1) && !is.Numeric(il1))
    stop("bad input for argument 'il1'")
  if (length(il2) && !is.Numeric(il2))
    stop("bad input for argument 'il2'")




  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("'nsimEIM' should be an integer greater than 10")


  new("vglmff",
  blurb = c("Mixture of two univariate exponentials\n\n",
            "Links:    ",
            namesof("phi",     lphi,    earg = ephi, tag = FALSE), ", ", 
            namesof("lambda1", llambda, earg = el1 , tag = FALSE), ", ",
            namesof("lambda2", llambda, earg = el2 , tag = FALSE), "\n",
            "Mean:     phi / lambda1 + (1 - phi) / lambda2\n"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),



  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("phi", "lambda1", "lambda2"),
         nsimEIM = .nsimEIM ,
         lphi      = .lphi   ,
         llambda1      = .llambda    ,
         llambda2      = .llambda    ,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM,
           .lphi = lphi,
           .llambda = llambda
         ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    predictors.names <-
      c(namesof("phi",     .lphi    , earg = .ephi , tag = FALSE),
        namesof("lambda1", .llambda , earg = .el1  , tag = FALSE),
        namesof("lambda2", .llambda , earg = .el2  , tag = FALSE))

    if (!length(etastart)) {
      qy <- quantile(y, prob =  .qmu)
      init.phi <-     rep(if (length(.iphi)) .iphi else 0.5, length = n)
      init.lambda1 <- rep(if (length(.il1)) .il1 else 1/qy[1], length = n)
      init.lambda2 <- rep(if (length(.il2)) .il2 else 1/qy[2], length = n)
      if (!length(etastart))  
        etastart <- cbind(theta2eta(init.phi,     .lphi    , earg = .ephi ),
                          theta2eta(init.lambda1, .llambda , earg = .el1  ),
                          theta2eta(init.lambda2, .llambda , earg = .el2  ))
      }
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2,
           .iphi = iphi, .il1 = il1, .il2 = il2,
           .qmu = qmu))),
  linkinv = eval(substitute(function(eta, extra = NULL){
    phi     <- eta2theta(eta[, 1], link = .lphi    , earg = .ephi )
    lambda1 <- eta2theta(eta[, 2], link = .llambda , earg = .el1  )
    lambda2 <- eta2theta(eta[, 3], link = .llambda , earg = .el2  )
    phi / lambda1 + (1 - phi) / lambda2
  }, list(.lphi = lphi, .llambda = llambda,
          .ephi = ephi, .el1 = el1, .el2 = el2 ))),
  last = eval(substitute(expression({
    misc$link <-
         c("phi" = .lphi , "lambda1" = .llambda , "lambda2" = .llambda )

    misc$earg <-
      list("phi" = .ephi , "lambda1" = .el1 ,     "lambda2" = .el2 )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list(.lphi = lphi, .llambda = llambda, .nsimEIM = nsimEIM,
           .ephi = ephi, .el1 = el1, .el2 = el2 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    phi     <- eta2theta(eta[, 1], link = .lphi    , earg = .ephi )
    lambda1 <- eta2theta(eta[, 2], link = .llambda , earg = .el1  )
    lambda2 <- eta2theta(eta[, 3], link = .llambda , earg = .el2  )

    f1 <- dexp(y, rate=lambda1)
    f2 <- dexp(y, rate=lambda2)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * log(phi*f1 + (1 - phi)*f2)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list(.lphi = lphi, .llambda = llambda,
          .ephi = ephi, .el1 = el1, .el2 = el2 ))),
  vfamily = c("mix2exp"),
  deriv = eval(substitute(expression({
    phi     <- eta2theta(eta[, 1], link = .lphi    , earg = .ephi )
    lambda1 <- eta2theta(eta[, 2], link = .llambda , earg = .el1  )
    lambda2 <- eta2theta(eta[, 3], link = .llambda , earg = .el2  )

    dphi.deta     <- dtheta.deta(phi,     link = .lphi    , earg = .ephi )
    dlambda1.deta <- dtheta.deta(lambda1, link = .llambda , earg = .el1  )
    dlambda2.deta <- dtheta.deta(lambda2, link = .llambda , earg = .el2  )

    f1 <- dexp(x = y, rate = lambda1)
    f2 <- dexp(x = y, rate = lambda2)
    pdf <- phi*f1 + (1 - phi)*f2
    df1.dlambda1 <- exp(-lambda1*y) - y * dexp(y, rate = lambda1)
    df2.dlambda2 <- exp(-lambda2*y) - y * dexp(y, rate = lambda2)
    dl.dphi <- (f1-f2) / pdf
    dl.dlambda1 <- phi * df1.dlambda1 / pdf
    dl.dlambda2 <- (1 - phi) * df2.dlambda2 / pdf

    c(w) * cbind(dl.dphi * dphi.deta,
                 dl.dlambda1 * dlambda1.deta,
                 dl.dlambda2 * dlambda2.deta)
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2 ))),
  weight = eval(substitute(expression({
    run.mean <- 0
    for (ii in 1:( .nsimEIM )) {
      ysim <- ifelse(runif(n) < phi, rexp(n, lambda1),
                                     rexp(n, lambda2))
      f1 <- dexp(x = ysim, rate=lambda1)
      f2 <- dexp(x = ysim, rate=lambda2)
      pdf <- phi*f1 + (1 - phi)*f2

      df1.dlambda1 <- exp(-lambda1*ysim) - ysim * dexp(ysim, rate = lambda1)
      df2.dlambda2 <- exp(-lambda2*ysim) - ysim * dexp(ysim, rate = lambda2)
      dl.dphi <- (f1-f2) / pdf
      dl.dlambda1 <- phi * df1.dlambda1 / pdf
      dl.dlambda2 <- (1 - phi) * df2.dlambda2 / pdf
      d2f1.dlambda12 <- ysim*(ysim*lambda1-2)*exp(-lambda1*ysim)
      d2f2.dlambda22 <- ysim*(ysim*lambda2-2)*exp(-lambda2*ysim)
      d2l.dphi2 <-  dl.dphi^2
      d2l.dlambda12 <- phi * (phi * df1.dlambda1^2 / pdf -
                       d2f1.dlambda12) / pdf
      d2l.dlambda22 <- (1 - phi) * ((1 - phi) * df2.dlambda2^2 / pdf -
                       d2f2.dlambda22) / pdf
      d2l.dlambda1lambda2 <- phi * (1 - phi) *
                             df1.dlambda1 * df2.dlambda2 / pdf^2
      d2l.dphilambda1 <- df1.dlambda1 * (phi*(f1-f2)/pdf - 1) / pdf
      d2l.dphilambda2 <- df2.dlambda2 * ((1 - phi)*(f1-f2)/pdf - 1) / pdf
      rm(ysim)

      temp3 <- matrix(0, n, dimm(M))
      temp3[, iam(1, 1, M = 3)] <- d2l.dphi2
      temp3[, iam(2, 2, M = 3)] <- d2l.dlambda12
      temp3[, iam(3, 3, M = 3)] <- d2l.dlambda22
      temp3[, iam(1, 2, M = 3)] <- d2l.dphilambda1
      temp3[, iam(1, 3, M = 3)] <- d2l.dphilambda2
      temp3[, iam(2, 3, M = 3)] <- d2l.dlambda1lambda2
      run.mean <- ((ii-1) * run.mean + temp3) / ii
    }
    wz <- if (intercept.only)
         matrix(colMeans(run.mean), n, dimm(M), byrow = TRUE) else
         run.mean

    dtheta.detas <- cbind(dphi.deta, dlambda1.deta, dlambda2.deta)
    index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    wz <- wz * dtheta.detas[, index0$row] *
               dtheta.detas[, index0$col]
    c(w) * wz
  }), list(.lphi = lphi, .llambda = llambda,
           .ephi = ephi, .el1 = el1, .el2 = el2,
           .nsimEIM = nsimEIM ))))
}




