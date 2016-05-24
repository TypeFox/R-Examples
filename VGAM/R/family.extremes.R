# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.

















rgev <- function(n, location = 0, scale = 1, shape = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
             stop("bad input for argument 'n'") else n

  if (!is.Numeric(location)) 
    stop("bad input for argument argument 'location'")
  if (!is.Numeric(shape))
    stop("bad input for argument argument 'shape'")

  ans <- numeric(use.n)
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n)
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)

  scase <- abs(shape) < sqrt( .Machine$double.eps )
  nscase <- sum(scase)
  if (use.n - nscase)
    ans[!scase] <- location[!scase] + scale[!scase] *
    ((-log(runif(use.n - nscase)))^(-shape[!scase]) -1) / shape[!scase]
  if (nscase)
    ans[scase] <- rgumbel(nscase, location = location[scase],
                          scale = scale[scase])
  ans[scale <= 0] <- NaN
  ans
}



 dgev <- function(x, location = 0, scale = 1, shape = 0, log = FALSE,
                  tolshape0 = sqrt( .Machine$double.eps ),
                  oobounds.log = -Inf, giveWarning = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (oobounds.log > 0)
    stop("bad input for argument 'oobounds.log'")

  if (!is.Numeric(tolshape0, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'tolshape0'")

  use.n <- max(length(x), length(location), length(scale), length(shape))
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n)
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)



  x         <- rep(x,          length.out = use.n)

  logdensity <- rep(log(0), length.out = use.n)
  scase <- (abs(shape) < tolshape0)
  nscase <- sum(scase)
  if (use.n - nscase) {
    zedd <- 1 + shape * (x - location) / scale # pmax(0, (1+shape*xc/scale))
    xok <- (!scase) & (zedd > 0)
    logdensity[xok] <- -log(scale[xok]) - zedd[xok]^(-1/shape[xok]) -
                       (1 + 1/shape[xok]) * log(zedd[xok])
    outofbounds <- (!scase) & (zedd <= 0)
    if (any(outofbounds)) {
      logdensity[outofbounds] <- oobounds.log
      no.oob <- sum(outofbounds)
      if (giveWarning)
        warning(no.oob, " observation",
                ifelse(no.oob > 1, "s are", " is"), " out of bounds")
    }
  }
  if (nscase) {
    logdensity[scase] <- dgumbel(x[scase], location = location[scase],
                                 scale = scale[scase], log = TRUE)
  }

  logdensity[scale <= 0] <- NaN

  logdensity[is.infinite(x)] <- log(0)  # 20141209 KaiH

  if (log.arg) logdensity else exp(logdensity)
}



pgev <- function(q, location = 0, scale = 1, shape = 0,
                 lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  

  use.n <- max(length(q), length(location), length(scale), length(shape))
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n)
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)
  if (length(q)        != use.n)
    q        <- rep(q,            length.out = use.n)

  scase0 <- abs(shape) < sqrt( .Machine$double.eps )  # Effectively 0
  zedd <- (q - location) / scale
  use.zedd <- pmax(0, 1 + shape * zedd)

  if (lower.tail) {
    if (log.p) {
      ans <- -use.zedd^(-1 / shape)
    } else {
      ans <- exp(-use.zedd^(-1 / shape))
    }
  } else {
    if (log.p) {
      ans <- log(-expm1(-use.zedd^(-1 / shape)))
    } else {
      ans <- -expm1(-use.zedd^(-1 / shape))
    }
  }

  if (any(scase0)) {
    ans[scase0] <- pgumbel(q[scase0], location = location[scase0],
                           scale = scale[scase0],
                           lower.tail = lower.tail, log.p = log.p)
  }

  ans[scale <= 0] <- NaN
  ans
}



qgev <- function(p, location = 0, scale = 1, shape = 0,
                 lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
 


  use.n <- max(length(p), length(location), length(scale), length(shape))
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n)
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)
  if (length(p)        != use.n)
    p        <- rep(p,            length.out = use.n)


  scase0 <- abs(shape) < sqrt( .Machine$double.eps )
  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- location + scale * ((-ln.p)^(-shape) - 1) / shape
      ans[ln.p > 0] <- NaN
    } else {
      ans <- location + scale * ((-log(p))^(-shape) - 1) / shape
      ans[p == 1] <-  Inf
      ans[p >  1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- location + scale * ((-log1p(-exp(ln.p)))^(-shape) - 1) / shape
      ans[ln.p > 0] <- NaN
    } else {
      ans <- location + scale * ((-log1p(-p))^(-shape) - 1) / shape
      ans[p == 1] <-  Inf
      ans[p >  1] <- NaN
      ans[p <  0] <- NaN
    }
  }

  if (any(scase0))
    ans[scase0] <- qgumbel(p[scase0], location = location[scase0],
                           scale = scale[scase0],
                           lower.tail = lower.tail, log.p = log.p)
  ans[scale <= 0] <- NaN
  ans
}





 gev <- function(
          llocation = "identitylink",
          lscale = "loge",
          lshape = logoff(offset = 0.5),
          percentiles = c(95, 99),
          iscale = NULL, ishape = NULL,
          imethod = 1, gshape = c(-0.45, 0.45),
          tolshape0 = 0.001,
          type.fitted = c("percentiles", "mean"),
          giveWarning = TRUE,
          zero = c("scale", "shape")) {



  type.fitted <- match.arg(type.fitted,
                           c("percentiles", "mean"))[1]

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")





  if (!is.logical(giveWarning) || length(giveWarning) != 1)
    stop("bad input for argument 'giveWarning'")

  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")



  if (length(percentiles) &&
     (!is.Numeric(percentiles, positive = TRUE) ||
      max(percentiles) >= 100))
    stop("bad input for argument 'percentiles'")

  if (!is.Numeric(imethod, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
     imethod > 2.5)
    stop("argument 'imethod' must be 1 or 2")
  if (length(ishape) && !is.Numeric(ishape))
      stop("bad input for argument 'ishape'")

  if (!is.Numeric(tolshape0, length.arg = 1, positive = TRUE) ||
      tolshape0 > 0.1)
    stop("bad input for argument 'tolshape0'")
  if (!is.Numeric(gshape, length.arg = 2) ||
      gshape[1] >= gshape[2])
    stop("bad input for argument 'gshape'")


  new("vglmff",
  blurb = c("Generalized extreme value distribution\n",
            "Links:    ",
            namesof("location", llocat, elocat), ", ", 
            namesof("scale",    lscale, escale), ", ",
            namesof("shape",    lshape, eshape)),
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
         parameters.names = c("location", "scale", "shape"),
         llocation = .llocat ,
         lscale    = .lscale ,
         lshape    = .lshape ,
         type.fitted = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocation, .lscale = lscale, .lshape = lshape,
           .type.fitted = type.fitted ))),


  initialize = eval(substitute(expression({
    M1 <- extra$M1 <- 3
    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$type.fitted <- .type.fitted
    extra$M1 <- M1


    mynames1  <- "location"
    mynames2  <- "scale"
    mynames3  <- "shape"
    llocat <- .llocat
    lscale <- .lscale
    lshape <- .lshape


    predictors.names <- c(
      namesof(mynames1, .llocat , .elocat , short = TRUE),
      namesof(mynames2, .lscale , .escale , short = TRUE),
      namesof(mynames3, .lshape , .eshape , short = TRUE))




    y <- as.matrix(y)





    if (ncol(y) > 1)
      y <- -t(apply(-y, 1, sort, na.last = TRUE))







    r.vec <- rowSums(cbind(!is.na(y)))


    if (any(r.vec == 0))
      stop("A row contains all missing values")

    extra$percentiles <- .percentiles
    if (!length(etastart)) {
      init.sig <- if (length( .iscale ))
        rep( .iscale, length.out = nrow(y)) else NULL
      init.xi <- if (length( .ishape ))
        rep( .ishape, length.out = nrow(y)) else NULL
      LIST.lshape <- .lshape

      if ( .lshape == "extlogit" && length(init.xi) &&
          (any(init.xi <= LIST.lshape$min |
               init.xi >= LIST.lshape$max)))
          stop("bad input for an argument in 'lshape'")

      if ( .imethod == 1) {
        nvector <- 4:10  # Arbitrary; could be made an argument
        ynvector <- quantile(y[, 1], probs = 1 - 1/nvector)
          objecFunction <- -Inf  # Actually the log-likelihood
          est.sigma <- !length(init.sig)
          gshape <- .gshape
          temp234 <- if (length(init.xi)) init.xi[1] else
                         seq(gshape[1], gshape[2], length.out = 12)
          for (shapeTry in temp234) {
              xvec <- if (abs(shapeTry) < .tolshape0) log(nvector) else
                      (nvector^shapeTry - 1) / shapeTry
              fit0 <- lsfit(x = xvec, y = ynvector, intercept = TRUE)
              sigmaTry <- if (est.sigma)
                rep(fit0$coef["X"], length.out = nrow(y)) else
                init.sig
              LocatTry <- rep(fit0$coef["Intercept"], length.out = nrow(y))
              llTry <- egev(giveWarning =
               FALSE)@loglikelihood(mu = NULL, y = y[, 1], w = w,
                                    residuals = FALSE,
                                    eta =
               cbind(theta2eta(LocatTry, .llocat , .elocat ),
                     theta2eta(sigmaTry, .lscale , .escale ),
                     theta2eta(shapeTry, .lshape , .eshape )))
                if (llTry >= objecFunction) {
                  if (est.sigma)
                    init.sig <- sigmaTry
                  init.mu <- rep(LocatTry, length.out = nrow(y))
                  objecFunction <- llTry
                  bestxi <- shapeTry
                }
            }
            if (!length(init.xi))
                init.xi <- rep(bestxi, length.out = nrow(y))
      } else {
        init.xi <- rep(0.05, length.out = nrow(y))
        if (!length(init.sig))
          init.sig <- rep(sqrt(6 * var(y[, 1]))/pi,
                          length.out = nrow(y))
        EulerM <- -digamma(1)
        init.mu <- rep(median(y[, 1]) - EulerM*init.sig,
                       length.out = nrow(y))
      }

      bad <- ((1 + init.xi*(y-init.mu)/init.sig) <= 0)
      if (fred <- sum(bad)) {
        warning(fred, "observations violating boundary constraints ",
                "while initializing. Taking corrective action")
        init.xi[bad] <- ifelse(y[bad] > init.mu[bad], 0.1, -0.1)
      }

      etastart <-
        cbind(theta2eta(init.mu,  .llocat , .elocat ),
              theta2eta(init.sig, .lscale , .escale ),
              theta2eta(init.xi,  .lshape , .eshape ))
    }
  }), list( 
            .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,
            .ishape = ishape, .iscale = iscale,

            .gshape = gshape, .type.fitted = type.fitted,
            .percentiles = percentiles,
            .tolshape0 = tolshape0,
            .imethod = imethod, .giveWarning = giveWarning ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Locat <- eta2theta(eta[, 1], .llocat , .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , .escale )
    shape <- eta2theta(eta[, 3], .lshape , .eshape )


    type.fitted <-
      if (length(extra$type.fitted)) {
        extra$type.fitted
      } else {
        warning("cannot find 'type.fitted'. Returning 'percentiles'.")
        "percentiles"
      }

    type.fitted <- match.arg(type.fitted,
                             c("percentiles", "mean"))[1]



    is.zero <- (abs(shape) < .tolshape0 )

    cent <- extra$percentiles
    LP <- length(cent)
    if (type.fitted == "percentiles" &&  # Upward compatibility:
        LP > 0) {
      fv <- matrix(NA_real_, nrow(eta), LP)
      for (ii in 1:LP) {
        yp <- -log(cent[ii] / 100)
        fv[!is.zero, ii] <- Locat[!is.zero] - sigma[!is.zero] *
                            (1 - yp^(-shape[!is.zero])) / shape[!is.zero]
        fv[ is.zero, ii] <- Locat[ is.zero] - sigma[ is.zero] *
                            log(yp)
      }
      dimnames(fv) <- list(dimnames(eta)[[1]],
                           paste(as.character(cent), "%", sep = ""))
    } else {
      EulerM <- -digamma(1)
      fv <- Locat + sigma * EulerM  # When shape = 0, is Gumbel
      fv[!is.zero] <- Locat[!is.zero] + sigma[!is.zero] *
                      (gamma(1-shape[!is.zero])-1) / shape[!is.zero]
      fv[shape >= 1] <- NA  # Mean exists only if shape < 1.
    }
    fv
  }, list(
            .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,

            .tolshape0 = tolshape0 ))),
  last = eval(substitute(expression({
    misc$earg <- vector("list", M)
    names(misc$earg) <- c(mynames1, mynames2, mynames3)
    misc$earg[[1]] <- .elocat
    misc$earg[[2]] <- .escale
    misc$earg[[3]] <- .eshape

    misc$link <-       c( .llocat , .lscale , .lshape )
    names(misc$link) <- c(mynames1, mynames2, mynames3)

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- FALSE





    misc$true.mu <- !length( .percentiles)  # @fitted is not a true mu
    misc$percentiles <- .percentiles
    misc$expected <- TRUE
    misc$tolshape0 <- .tolshape0
    if (ncol(y) == 1)
      y <- as.vector(y)
    if (any(shape < -0.5))
      warning("some values of the shape parameter are less than -0.5")
  }), list(
            .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,

            .tolshape0 = tolshape0, .percentiles = percentiles ))),
  loglikelihood = eval(substitute(
  function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Locat <- eta2theta(eta[, 1], .llocat , .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , .escale )
    shape <- eta2theta(eta[, 3], .lshape , .eshape )


    is.zero <- (abs(shape) < .tolshape0 )
    zedd <- (y-Locat) / sigma
    r.vec <- rowSums(cbind(!is.na(y)))
    A <- 1 + shape * (y-Locat)/sigma
    ii <- 1:nrow(eta)
    A1 <- A[cbind(ii, r.vec)]
    mytolerance <- 0  # .Machine$double.eps
    if (any(bad <- (A1 <= mytolerance), na.rm = TRUE)) {
      if ( .giveWarning )
        warning("There are", sum(bad), "range violations in @loglikelihood")
      cat("There are", sum(bad),
          "range violations in @loglikelihood\n")
      flush.console()
    }
    igev <- !is.zero &  !bad
    igum <-  is.zero &  !bad
    pow <- 1 + 1/shape[igev]
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {





      old.answer <-
        sum(bad) * (-1.0e10) +
        sum(w[igum] * (-r.vec[igum] * log(sigma[igum]) -
                  exp(-zedd[igum, r.vec[igum]]) -
                  rowSums(cbind(zedd)[igum, , drop = FALSE], na.rm = TRUE))) +
        sum(w[igev] * (-r.vec[igev] * log(sigma[igev]) -
                       pow * rowSums(cbind(log(A[igev])), na.rm = TRUE) -
                       A1[igev]^(-1/shape[igev])))

        new.answer <-
          sum(w * dgev(x = y, location = Locat, scale = sigma, shape = shape,
                       tolshape0 = .tolshape0 ,
                       giveWarning = .giveWarning,
                       log = TRUE, oobounds.log = -1.0e04))
        check0 <- old.answer - new.answer
      old.answer
      new.answer
    }
  }, list( 
            .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,

           .giveWarning = giveWarning, .tolshape0 = tolshape0 ))),
  vfamily = c("gev", "vextremes"),
  deriv = eval(substitute(expression({
    M1 <- 3
    r.vec <- rowSums(cbind(!is.na(y)))

    Locat <- eta2theta(eta[, 1], .llocat , .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , .escale )
    shape <- eta2theta(eta[, 3], .lshape , .eshape )

    dmu.deta <- dtheta.deta(Locat, .llocat , .elocat )
    dsi.deta <- dtheta.deta(sigma, .lscale , .escale )
    dxi.deta <- dtheta.deta(shape, .lshape , .eshape )

    is.zero <- (abs(shape) < .tolshape0 )
    ii <- 1:nrow(eta)
    zedd <- (y-Locat) / sigma
    A <- 1 + shape * zedd
    dA.dxi <- zedd                   # matrix
    dA.dmu <- -shape/sigma           # vector
    dA.dsigma <- -shape*zedd/sigma   # matrix
    pow <- 1 + 1/shape
    A1 <- A[cbind(ii, r.vec)]

    AAr1 <- dA.dmu/(shape * A1^pow) -
           pow * rowSums(cbind(dA.dmu/A), na.rm = TRUE)
    AAr2 <- dA.dsigma[cbind(ii,r.vec)] / (shape * A1^pow) -
           pow * rowSums(cbind(dA.dsigma/A), na.rm = TRUE)
    AAr3 <- 1/(shape * A1^pow) -
           pow * rowSums(cbind(dA.dsigma/A), na.rm = TRUE)
    dl.dmu <- AAr1
    dl.dsi <- AAr2 - r.vec/sigma
    dl.dxi <- rowSums(cbind(log(A)), na.rm = TRUE)/shape^2 -
              pow * rowSums(cbind(dA.dxi/A), na.rm = TRUE) -
              (log(A1) / shape^2 -
              dA.dxi[cbind(ii,r.vec)] / (shape*A1)) * A1^(-1/shape)

    if (any(is.zero)) {
      zorro <- c(zedd[cbind(1:n,r.vec)])
      zorro <- zorro[is.zero]
      ezedd <- exp(-zorro)
      dl.dmu[is.zero] <- (1-ezedd) / sigma[is.zero]
      dl.dsi[is.zero] <- (zorro *  (1 - ezedd) - 1) / sigma[is.zero]
      dl.dxi[is.zero] <-  zorro * ((1 - ezedd) * zorro / 2 - 1)
    }

    c(w) * cbind(dl.dmu * dmu.deta,
                 dl.dsi * dsi.deta,
                 dl.dxi * dxi.deta)
  }), list(
            .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,

            .tolshape0 = tolshape0 ))),

  weight = eval(substitute(expression({
    kay <- -shape
    dd <- digamma(r.vec-kay+1)
    ddd <- digamma(r.vec+1)  # Unnecessarily evaluated at each iteration
    temp13 <- -kay * dd + (kay^2 - kay + 1) / (1-kay)
    temp33 <- 1 - 2 * kay * ddd +
              kay^2 * (1 + trigamma(r.vec+1) + ddd^2)
    temp23 <- -kay * dd + (1+(1-kay)^2) / (1-kay)
    GR.gev <- function(j, ri, kay) gamma(ri - j*kay + 1) / gamma(ri)
    tmp2 <- (1-kay)^2 * GR.gev(2, r.vec, kay)  # Latter is GR2
    tmp1 <- (1-2*kay) * GR.gev(1, r.vec, kay)  # Latter is GR1
    k0 <- (1-2*kay)
    k1 <- k0 * kay
    k2 <- k1 * kay
    k3 <- k2 * kay  # kay^3 * (1-2*kay)

    wz <- matrix(NA_real_, n, 6)
    wz[, iam(1, 1, M)] <- tmp2 / (sigma^2 * k0)
    wz[, iam(1, 2, M)] <- (tmp2 - tmp1) / (sigma^2 * k1)
    wz[, iam(1, 3, M)] <- (tmp1 * temp13 - tmp2) / (sigma * k2)
    wz[, iam(2, 2, M)] <- (r.vec*k0 - 2*tmp1 + tmp2) / (sigma^2 * k2)
    wz[, iam(2, 3, M)] <- (r.vec*k1*ddd + tmp1 *
                       temp23 - tmp2 - r.vec*k0) / (sigma * k3)
    wz[, iam(3, 3, M)] <- (2*tmp1*(-temp13) + tmp2 +
                         r.vec*k0*temp33)/(k3*kay)

    if (any(is.zero)) {
      if (ncol(y) > 1)
        stop("cannot handle shape == 0 with a multivariate response")

      EulerM <- -digamma(1)
      wz[is.zero, iam(2, 2, M)] <- (pi^2/6 + (1-EulerM)^2)/sigma[is.zero]^2
      wz[is.zero, iam(3, 3, M)] <- 2.4236
      wz[is.zero, iam(1, 2, M)] <-
        (digamma(2) + 2 * (EulerM-1)) / sigma[is.zero]^2
      wz[is.zero, iam(1, 3, M)] <-
       -(trigamma(1) / 2 + digamma(1) * (digamma(1)/2+1)) / sigma[is.zero]
      wz[is.zero, iam(2, 3, M)] <-
        (-dgammadx(2, 3)/6 + dgammadx(1, 1) +
                           2*dgammadx(1, 2) +
                           2*dgammadx(1, 3)/3) / sigma[is.zero]

      if (FALSE ) {
        wz[, iam(1, 2, M)] <- 2 * r.vec / sigma^2
        wz[, iam(2, 2, M)] <- -4 * r.vec * digamma(r.vec+1) + 2 * r.vec +
                              (4 * dgammadx(r.vec+1, deriv.arg = 1) -
          3 * dgammadx(r.vec+1,
                       deriv.arg = 2)) / gamma(r.vec)  # Not checked
      }
    }

    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] * dmu.deta^2
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] * dsi.deta^2
    wz[, iam(3, 3, M)] <- wz[, iam(3, 3, M)] * dxi.deta^2
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)] * dmu.deta *   dsi.deta
    wz[, iam(1, 3, M)] <- wz[, iam(1, 3, M)] * dmu.deta * (-dxi.deta)
    wz[, iam(2, 3, M)] <- wz[, iam(2, 3, M)] * dsi.deta * (-dxi.deta)
    c(w) * wz
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape ))))
}




dgammadx <- function(x, deriv.arg = 1) {
  if (deriv.arg == 0) {
    gamma(x)
  } else if (deriv.arg == 1) {
    digamma(x) * gamma(x)
  } else if (deriv.arg == 2) {
    gamma(x) * (trigamma(x) + digamma(x)^2)
  } else if (deriv.arg == 3) {
    gamma(x) * (psigamma(x, deriv = 2) +
    2 * digamma(x) * trigamma(x)) +
    Recall(x, deriv.arg = 1) * (trigamma(x) + digamma(x)^2)
  } else if (deriv.arg == 4) {
      Recall(x, deriv.arg = 2) * (trigamma(x) + digamma(x)^2) +
  2 * Recall(x, deriv.arg = 1) * (psigamma(x, deriv = 2) +
      2*digamma(x) * trigamma(x)) +
      gamma(x) * (psigamma(x, deriv = 3) + 2*trigamma(x)^2 +
               2 * digamma(x) * psigamma(x, deriv = 2))
  } else {
    stop("cannot handle 'deriv' > 4")
  }
}






 egev <- function(llocation = "identitylink",
                  lscale = "loge",
                  lshape = logoff(offset = 0.5),
                  percentiles = c(95, 99),
                  iscale = NULL, ishape = NULL,
                  imethod = 1, gshape = c(-0.45, 0.45),
                  tolshape0 = 0.001,
                  type.fitted = c("percentiles", "mean"),
                  giveWarning = TRUE,
                  zero = c("scale", "shape")) {
  if (!is.logical(giveWarning) || length(giveWarning) != 1)
    stop("bad input for argument 'giveWarning'")
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  type.fitted <- match.arg(type.fitted,
                           c("percentiles", "mean"))[1]

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")




  if (!is.Numeric(gshape, length.arg = 2) ||
      gshape[1] >= gshape[2])
      stop("bad input for argument 'gshape'")
    if (length(percentiles) && 
      (!is.Numeric(percentiles, positive = TRUE) ||
       max(percentiles) >= 100))
      stop("bad input for argument 'percentiles'")
    if (!is.Numeric(imethod, length.arg = 1,
                    positive = TRUE, integer.valued = TRUE) ||
       imethod > 2.5)
      stop("argument 'imethod' must be 1 or 2")
    if (length(ishape) && !is.Numeric(ishape))
      stop("bad input for argument 'ishape'")
    if (!is.Numeric(tolshape0, length.arg = 1,
                    positive = TRUE) ||
        tolshape0 > 0.1)
      stop("bad input for argument 'tolshape0'")


  new("vglmff",
  blurb = c("Generalized extreme value distribution\n",
          "Links:    ",
          namesof("location", link = llocat, earg = elocat), ", ", 
          namesof("scale",    link = lscale, earg = escale), ", ",
          namesof("shape",    link = lshape, earg = eshape)),
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
         parameters.names = c("location", "scale", "shape"),
         llocation = .llocat ,
         lscale    = .lscale ,
         lshape    = .lshape ,
         type.fitted = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocation, .lscale = lscale, .lshape = lshape,
           .type.fitted = type.fitted ))),


  initialize = eval(substitute(expression({
    M1 <- extra$M1 <- 3
    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$type.fitted <- .type.fitted
    extra$M1 <- M1

    predictors.names <-
     c(namesof("location", .llocat , earg = .elocat , short = TRUE),
       namesof("scale",    .lscale , earg = .escale , short = TRUE),
       namesof("shape",    .lshape , earg = .eshape , short = TRUE))



    if (ncol(as.matrix(y)) != 1)
      stop("response must be a vector or one-column matrix")



    if (!length(etastart)) {
      init.sig <- if (length( .iscale ))
                  rep( .iscale , length.out = length(y)) else NULL
      init.xi  <- if (length( .ishape ))
                  rep( .ishape , length.out = length(y)) else NULL
      eshape <- .eshape
      if ( .lshape == "extlogit" && length(init.xi) && 
         (any(init.xi <= eshape$min | init.xi >= eshape$max)))
        stop("bad input for argument 'eshape'")
      if ( .imethod == 1) {
        nvector <- 4:10   # Arbitrary; could be made an argument
        ynvector <- quantile(y, probs = 1-1/nvector)
        objecFunction <- -Inf   # Actually the log-likelihood
        est.sigma <- !length(init.sig)
        gshape <- .gshape
        temp234 <- if (length(init.xi)) init.xi[1] else
                   seq(gshape[1], gshape[2], length.out = 12)
        for (xi.try in temp234) {
          xvec <- if (abs(xi.try) < .tolshape0 ) log(nvector) else
                  (nvector^xi.try - 1) / xi.try
          fit0 <- lsfit(x = xvec, y=ynvector, intercept = TRUE)
          if (est.sigma) {
            sigmaTry <- rep(fit0$coef["X"], length.out = length(y))
          } else { 
            sigmaTry <- init.sig
          }
          muTry <- rep(fit0$coef["Intercept"], length.out = length(y))
          llTry <- egev(giveWarning = FALSE)@loglikelihood(mu = NULL,
                                             y = y, w = w,
                                             residuals = FALSE,
          eta <- cbind(theta2eta(muTry,    .llocat , earg = .elocat ),
                       theta2eta(sigmaTry, .lscale , earg = .escale ), 
                       theta2eta(xi.try,   .lshape , earg = .eshape )))
          if (llTry >= objecFunction) {
            if (est.sigma)
              init.sig <- sigmaTry
            init.mu <- rep(muTry, length.out = length(y))
            objecFunction <- llTry
            bestxi <- xi.try
          }
        }
        if (!length(init.xi))
          init.xi <- rep(bestxi, length.out = length(y))

        } else {
          init.xi <- rep(if (length(init.xi)) init.xi else 0.05,
                         length.out = length(y))
          if (!length(init.sig))
            init.sig <- rep(sqrt(6 * var(y)) / pi, length.out = length(y))
          EulerM <- -digamma(1)
          init.mu <- rep(median(y) - EulerM * init.sig,
                         length.out = length(y))
        }
        bad <- (1 + init.xi * (y - init.mu) / init.sig <= 0)
        if (fred <- sum(bad, na.rm = TRUE)) {
          warning(fred, "observations violating boundary constraints ",
                  "while initializing. Taking corrective action")
          init.xi[bad] <- ifelse(y[bad] > init.mu[bad], 0.01, -0.01)
        }

      extra$percentiles <- .percentiles

      etastart <-
        cbind(theta2eta(init.mu,  .llocat ,    earg = .elocat ),
              theta2eta(init.sig, .lscale ,    earg = .escale ), 
              theta2eta(init.xi,  .lshape ,    earg = .eshape ))
    }
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,
            .percentiles = percentiles, .tolshape0 = tolshape0,
            .imethod = imethod, .type.fitted = type.fitted,
            .giveWarning= giveWarning,
            .iscale = iscale, .ishape = ishape, .gshape = gshape ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    loc   <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )
    xi    <- eta2theta(eta[, 3], .lshape , earg = .eshape )
    type.fitted <-
      if (length(extra$type.fitted)) {
        extra$type.fitted
      } else {
        warning("cannot find 'type.fitted'. Returning 'percentiles'.")
        "percentiles"
      }

    type.fitted <- match.arg(type.fitted,
                             c("percentiles", "mean"))[1]



    is.zero <- (abs(xi) < .tolshape0 )
    cent <- extra$percentiles
    LP <- length(cent)
    if (type.fitted == "percentiles" &&  # Upward compatibility:
        LP > 0) {
      fv <- matrix(NA_real_, nrow(eta), LP)
      for (ii in 1:LP) {
        yp <- -log(cent[ii] / 100)
        fv[!is.zero, ii] <- loc[!is.zero] - sigma[!is.zero] *
                            (1 - yp^(-xi[!is.zero])) / xi[!is.zero]
        fv[is.zero, ii] <- loc[is.zero] - sigma[is.zero] * log(yp)
      }
      dimnames(fv) <- list(dimnames(eta)[[1]],
                           paste(as.character(cent), "%", sep = ""))
    } else {
      EulerM <- -digamma(1)
      fv <- loc + sigma * EulerM  # When xi = 0, is Gumbel
      fv[!is.zero] <- loc[!is.zero] + sigma[!is.zero] *
                    (gamma(1-xi[!is.zero])-1) / xi[!is.zero]
      fv[xi >= 1] <- NA  # Mean exists only if xi < 1.
    }
    fv
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape,
           .type.fitted = type.fitted, .tolshape0 = tolshape0 ))),
  last = eval(substitute(expression({
    misc$links <-   c(location = .llocat, 
                      scale    = .lscale ,
                      shape    = .lshape)

    misc$earg <- list(location = .elocat,
                      scale    = .escale,
                      shape    = .eshape)


    misc$true.mu <- !length( .percentiles)  # @fitted is not a true mu
    misc$percentiles <- .percentiles
    misc$tolshape0 <- .tolshape0
    misc$expected <- TRUE 
    if (any(xi < -0.5))
      warning("some values of the shape parameter are less than -0.5")
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,
            .tolshape0 = tolshape0,  .percentiles = percentiles ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    mmu   <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )
    xi    <- eta2theta(eta[, 3], .lshape , earg = .eshape )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dgev(x = y, location = mmu, scale = sigma,
                             shape = xi, tolshape0 = .tolshape0 ,
                             log = TRUE, oobounds.log = -1.0e04,
                             giveWarning = .giveWarning )
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
           .elocat = elocat, .escale = escale, .eshape = eshape,
           .giveWarning = giveWarning, .tolshape0 = tolshape0 ))),
  vfamily = c("egev", "vextremes"),
  deriv = eval(substitute(expression({
    Locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )
    xi    <- eta2theta(eta[, 3], .lshape , earg = .eshape)
    is.zero <- (abs(xi) < .tolshape0)
    zedd <- (y - Locat) / sigma
    A <- 1 + xi * zedd
    dA.dxi <- zedd
    dA.dmu <- -xi / sigma
    dA.dsigma <- -xi * zedd / sigma
    pow <- 1 + 1/xi
    if (any(bad <- A <= 0, na.rm = TRUE))
      stop(sum(bad, na.rm = TRUE),
           " observations violating boundary constraints in '@deriv'")
    AA <- 1/(xi*A^pow)- pow/A 
    dl.dmu <- dA.dmu * AA
    dl.dsi <- dA.dsigma * AA - 1/sigma
    dl.dxi <-  log(A)/xi^2 - pow * dA.dxi / A -
           (log(A)/xi^2 - dA.dxi /(xi*A)) * A^(-1/xi)
    if (any(is.zero)) {
      ezedd <- exp(-zedd[is.zero])
      dl.dmu[is.zero] <- (1 - ezedd) / sigma[is.zero]
      dl.dsi[is.zero] <- (zedd[is.zero] *
                        (1 - ezedd) - 1) / sigma[is.zero]
      dl.dxi[is.zero] <- zedd[is.zero] *
                        ((1 - ezedd) * zedd[is.zero] / 2 - 1)
    }
    dmu.deta <- dtheta.deta(Locat, .llocat , earg = .elocat )
    dsi.deta <- dtheta.deta(sigma, .lscale , earg = .escale )
    dxi.deta <- dtheta.deta(xi, .lshape , earg = .eshape)
    c(w) * cbind(dl.dmu * dmu.deta,
                 dl.dsi * dsi.deta,
                 dl.dxi * dxi.deta)
  }), list( .llocat = llocat, .lscale = lscale, .lshape = lshape,
            .elocat = elocat, .escale = escale, .eshape = eshape,
            .tolshape0 = tolshape0 ))),
  weight = eval(substitute(expression({
    bad <- A <= 0
    if (any(bad, na.rm = TRUE))
      stop(sum(bad, na.rm = TRUE),
           " observations violating boundary constraints in '@weight'")
    kay <- -xi  # for the formulae 
    kay[abs(kay-0.5) < .tolshape0] <- 0.501
    temp100 <- gamma(2-kay)
    pp <- (1-kay)^2 * gamma(1-2*kay)  # gamma(0) is undefined so kay != 0.5
    qq <- temp100 * (digamma(1-kay) - (1-kay)/kay)
    wz <- matrix(NA_real_, n, 6)
    wz[, iam(1, 1, M)] <- pp / sigma^2
    wz[, iam(2, 2, M)] <- (1 - 2*temp100 + pp) / (sigma * kay)^2
    EulerM <- -digamma(1)
    wz[, iam(3, 3, M)] <- (pi^2 / 6 + (1-EulerM-1/kay)^2 +
                          (2*qq + pp/kay)/kay) / kay^2 
    wz[, iam(1, 2, M)] <- (pp - temp100) / (sigma^2 * kay)
    wz[, iam(1, 3, M)] <- -(qq + pp/kay) / (sigma * kay)
    wz[, iam(2, 3, M)] <- (1-EulerM - (1-temp100)/kay - qq -
                        pp/kay) / (sigma * kay^2)
    if (any(is.zero)) {
      wz[is.zero, iam(2, 2, M)] <- (pi^2/6 + (1-EulerM)^2) / sigma[is.zero]^2
      wz[is.zero, iam(3, 3, M)] <- 2.4236
      wz[is.zero, iam(1, 2, M)] <- (digamma(2) + 2*(EulerM-1)) / sigma[is.zero]^2
      wz[is.zero, iam(1, 3, M)] <- -(trigamma(1)/2 + digamma(1)*
                                     (digamma(1)/2+1))/sigma[is.zero]
      wz[is.zero, iam(2, 3, M)] <- (-dgammadx(2, 3)/6 + dgammadx(1, 1) +
                                   2*dgammadx(1, 2) +
                                   2*dgammadx(1, 3) / 3) / sigma[is.zero]
    }
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] * dmu.deta^2
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] * dsi.deta^2
    wz[, iam(3, 3, M)] <- wz[, iam(3, 3, M)] * dxi.deta^2
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)] * dmu.deta * dsi.deta
    wz[, iam(1, 3, M)] <- wz[, iam(1, 3, M)] * dmu.deta * (-dxi.deta)
    wz[, iam(2, 3, M)] <- wz[, iam(2, 3, M)] * dsi.deta * (-dxi.deta)
    c(w) * wz
  }), list( .eshape = eshape, .tolshape0 = tolshape0 ))))
}





rgumbel <- function(n, location = 0, scale = 1) {
  answer <- location - scale * log(-log(runif(n)))
  answer[scale <= 0] <- NaN
  answer
}


dgumbel <- function(x, location = 0, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  zedd <- (x - location) / scale
  logdensity <- -zedd - exp(-zedd) - log(scale)
  logdensity[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) logdensity else exp(logdensity)
}



qgumbel <- function(p, location = 0, scale = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- location - scale * log(-ln.p)
    } else {
      ans <- location - scale * log(-log(p))
      ans[p == 0] <- -Inf
      ans[p == 1] <-  Inf
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- location - scale * log(-log(-expm1(ln.p)))
      ans[ln.p > 0] <- NaN
    } else { 
      ans <- location - scale * log(-log1p(-p))
      ans[p == 0] <-  Inf
      ans[p == 1] <- -Inf
    }
  }
  ans[scale <= 0] <- NaN
  ans
}



pgumbel <- function(q, location = 0, scale = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ans <- -exp(-(q - location) / scale)
      ans[q <= -Inf] <- -Inf
      ans[q ==  Inf] <- 0
    } else {
      ans <- exp(-exp(-(q - location) / scale))
      ans[q <= -Inf] <- 0
      ans[q ==  Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log(-expm1(-exp(-(q - location) / scale)))
      ans[q <= -Inf] <- 0
      ans[q ==  Inf] <- -Inf
    } else {
      ans <- -expm1(-exp(-(q - location) / scale))
      ans[q <= -Inf] <- 1
      ans[q ==  Inf] <- 0
    }
  } 

  ans[scale <= 0] <- NaN
  ans
}




 gumbel <- function(llocation = "identitylink",
                    lscale = "loge",
                    iscale = NULL,
                    R = NA, percentiles = c(95, 99),
                    mpv = FALSE, zero = NULL) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (!is.logical(mpv) || length(mpv) != 1)
    stop("bad input for argument 'mpv'")

  if (length(percentiles) &&
     (!is.Numeric(percentiles, positive = TRUE) ||
      max(percentiles) >= 100))
    stop("bad input for argument 'percentiles'")


  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")



  new("vglmff",
  blurb = c("Gumbel distribution for extreme value regression\n",
            "Links:    ",
            namesof("location", llocat,  earg = elocat ), ", ",
            namesof("scale",    lscale, earg = escale )),
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
         llocation = .llocat ,
         lscale    = .lscale ,
         mpv = .mpv ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocation, .lscale = lscale,
           .mpv = mpv ))),


  initialize = eval(substitute(expression({

    predictors.names <-
    c(namesof("location", .llocat , earg = .elocat , short = TRUE),
      namesof("scale",    .lscale , earg = .escale , short = TRUE))


    y <- as.matrix(y)
    if (ncol(y) > 1)
      y <- -t(apply(-y, 1, sort, na.last = TRUE))



    r.vec <- rowSums(cbind(!is.na(y)))
    if (any(r.vec == 0))
      stop("There is at least one row of the response containing all NAs")


    if (ncol(y) > 1) {
      yiri <- y[cbind(1:nrow(y), r.vec)]
      sc.init <- if (is.Numeric( .iscale, positive = TRUE))
                .iscale else {3 * (rowMeans(y, na.rm = TRUE) - yiri)}
      sc.init <- rep(sc.init, length = nrow(y))
      sc.init[sc.init <= 0.0001] <- 1  # Used to be .iscale
      loc.init <- yiri + sc.init * log(r.vec)
    } else {
      sc.init <-  if (is.Numeric( .iscale, positive = TRUE))
                     .iscale else 1.1 * (0.01+sqrt(var(y)*6)) / pi
      sc.init <- rep(sc.init, length.out = n)
      EulerM <- -digamma(1)
      loc.init <- (y - sc.init * EulerM)
      loc.init[loc.init <= 0] <- min(y)
    }

    extra$R <- .R
    extra$mpv <- .mpv
    extra$percentiles <- .percentiles

    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(loc.init, .llocat , earg = .elocat ),
              theta2eta( sc.init, .lscale , earg = .escale ))
    }

  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
                              .iscale = iscale,
            .R = R, .mpv = mpv, .percentiles = percentiles ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    loc   <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )

    Percentiles <- extra$percentiles
    LP <- length(Percentiles)  # may be 0
    if (LP > 0) {
      mpv <- extra$mpv
      mu <- matrix(NA_real_, nrow(eta), LP + mpv)  # LP may be 0
      Rvec <- extra$R
      for (ii in 1:LP) {
        ci <- if (is.Numeric(Rvec))
              Rvec * (1 - Percentiles[ii] / 100) else
              -log(Percentiles[ii] / 100)
        mu[, ii] <- loc - sigma * log(ci)
      }
      if (mpv) 
        mu[, ncol(mu)] <- loc - sigma * log(log(2))


    dmn2 <- paste(as.character(Percentiles), "%", sep = "")
    if (mpv) 
      dmn2 <- c(dmn2, "MPV")
    dimnames(mu) <- list(dimnames(eta)[[1]], dmn2)
  } else {
    EulerM <- -digamma(1)
    mu <- loc + sigma * EulerM
  }
  mu
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale ))),

  last = eval(substitute(expression({
    misc$links <-   c(location = .llocat , scale = .lscale )

    misc$earg <- list(location = .elocat , scale = .escale )

    misc$R <- .R
    misc$mpv <- .mpv
    misc$true.mu <- !length( .percentiles )  # @fitted is not a true mu
    misc$percentiles <- .percentiles
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .percentiles = percentiles,
            .mpv = mpv, .R = R ))),
  vfamily = c("gumbel", "vextremes"),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    loc   <- eta2theta(eta[, 1], .llocat,  earg = .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )

    r.vec <- rowSums(cbind(!is.na(y)))
    yiri <- y[cbind(1:nrow(y), r.vec)]
    ans <- -r.vec * log(sigma) - exp( -(yiri-loc)/sigma )
    max.r.vec <- max(r.vec)
    for (jay in 1:max.r.vec) {
      index <- (jay <= r.vec)
      ans[index] <- ans[index] - (y[index,jay] - loc[index]) / sigma[index]
    }


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * ans
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale ))),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )

    r.vec <- rowSums(cbind(!is.na(y)))
    yiri <- y[cbind(1:nrow(y), r.vec)]
    yi.bar <- rowMeans(y, na.rm = TRUE)
    temp2 <- (yiri - locat) / sigma
    term2 <- exp(-temp2)

    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )
    dsigma.deta <- dtheta.deta(sigma, .lscale , earg = .escale )

    dl.dlocat <- (r.vec - term2) / sigma
    dl.dsigma <- (rowSums((y - locat) / sigma, na.rm = TRUE) -
                 r.vec - temp2 * term2) / sigma

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dsigma * dsigma.deta)
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  weight = eval(substitute(expression({
    temp6 <- digamma(r.vec)  # , integer = T
    temp5 <- digamma(1:max(r.vec))  # , integer=T
    temp5 <- matrix(temp5, n, max(r.vec), byrow = TRUE)
    temp5[col(temp5) > r.vec] <- 0
    temp5 <- temp5 %*% rep(1, ncol(temp5))

    wz <- matrix(NA_real_, n, dimm(M = 2))  # 3=dimm(M = 2)
    wz[, iam(1, 1, M)] <- r.vec / sigma^2
    wz[, iam(2, 1, M)] <- -(1 + r.vec * temp6) / sigma^2
    wz[, iam(2, 2, M)] <- (2*(r.vec+1)*temp6 + r.vec*(trigamma(r.vec) +
                          temp6^2) + 2 - r.vec - 2*temp5) / sigma^2

    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] * dlocat.deta^2
    wz[, iam(2, 1, M)] <- wz[, iam(2, 1, M)] * dsigma.deta * dlocat.deta
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] * dsigma.deta^2

    c(w) * wz
  }), list( .lscale = lscale ))))
}




rgpd <- function(n, location = 0, scale = 1, shape = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
             stop("bad input for argument 'n'") else n

  if (!is.Numeric(location))
    stop("bad input for argument 'location'")
  if (!is.Numeric(shape))
    stop("bad input for argument 'shape'")

  ans <- numeric(use.n)
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n); 
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)


  scase <- abs(shape) < sqrt( .Machine$double.eps )
  nscase <- sum(scase)
  if (use.n - nscase)
    ans[!scase] <- location[!scase] +
                   scale[!scase] *
       ((runif(use.n - nscase))^(-shape[!scase])-1) / shape[!scase]
  if (nscase)
    ans[scase] <- location[scase] - scale[scase] * log(runif(nscase))
  ans[scale <= 0] <- NaN
  ans
}



dgpd <- function(x, location = 0, scale = 1, shape = 0, log = FALSE,
                tolshape0 = sqrt( .Machine$double.eps ),
                oobounds.log = -Inf, giveWarning = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  if (oobounds.log > 0)
    stop("bad input for argument 'oobounds.log'")

  if (!is.Numeric(tolshape0, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'tolshape0'")


  L <- max(length(x), length(location), length(scale), length(shape))
  if (length(shape)    != L)
    shape    <- rep(shape,        length.out = L)
  if (length(location) != L)
    location <- rep(location,     length.out = L)
  if (length(scale)    != L)
    scale    <- rep(scale,        length.out = L)
  if (length(x)        != L)
    x        <- rep(x,            length.out = L)




  logdensity <- rep(log(0), length.out = L)
  scase <- abs(shape) < tolshape0
  nscase <- sum(scase)
  if (L - nscase) {
    zedd <- (x-location) / scale
    xok <- (!scase) & (zedd > 0) & (1 + shape*zedd > 0)
    logdensity[xok] <- -(1 + 1/shape[xok])*log1p(shape[xok]*zedd[xok]) -
                      log(scale[xok])
    outofbounds <- (!scase) & ((zedd <= 0) | (1 + shape*zedd <= 0))
    if (any(outofbounds)) {
      logdensity[outofbounds] <- oobounds.log
      no.oob <- sum(outofbounds)
      if (giveWarning)
        warning(no.oob, " observation",
                ifelse(no.oob > 1, "s are", " is"), " out of bounds")
    }
  }
  if (nscase) {
    xok <- scase & (x > location)
    logdensity[xok] <- -(x[xok] - location[xok]) / scale[xok] -
                      log(scale[xok])
    outofbounds <- scase & (x <= location)
    if (any(outofbounds)) {
        logdensity[outofbounds] <- oobounds.log
        no.oob <- sum(outofbounds)
        if (giveWarning)
          warning(no.oob, " observation",
                  ifelse(no.oob > 1, "s are", " is"), " out of bounds")
    }
  }

  logdensity[scale <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}



pgpd <- function(q, location = 0, scale = 1, shape = 0,
                 lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  use.n <- max(length(q), length(location), length(scale), length(shape))

  ans <- numeric(use.n)
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n)
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)
  if (length(q)        != use.n)
    q        <- rep(q,            length.out = use.n)

  zedd <- (q - location) / scale
  use.zedd <- pmax(zedd, 0)


  scase0 <- abs(shape) < sqrt( .Machine$double.eps )
  nscase0 <- sum(scase0)
  if (use.n - nscase0) {
    ans <- 1 - pmax(1 + shape * use.zedd, 0)^(-1/shape)
  }
  if (nscase0) {
    pos <- (zedd >= 0)
    ind9 <- ( pos & scase0)
    ans[ind9] <-  -expm1(-use.zedd[ind9])
    ind9 <- (!pos & scase0)
    ans[ind9] <- 0
  }
  ans[scale <= 0] <- NaN

  if (lower.tail) {
    if (log.p) log(ans) else ans
  } else {
    if (log.p) log1p(-ans) else 1-ans
  }
}



qgpd <- function(p, location = 0, scale = 1, shape = 0,
                 lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")

  if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  rm(log.p)
  
  if (lower.tail) {
    if (log.arg) p <- exp(p)
  } else {
    p <- if (log.arg) -expm1(p) else 1 - p
  }

  use.n <- max(length(p), length(location), length(scale), length(shape))

  ans <- numeric(use.n)
  if (length(shape)    != use.n)
    shape    <- rep(shape,        length.out = use.n)
  if (length(location) != use.n)
    location <- rep(location,     length.out = use.n)
  if (length(scale)    != use.n)
    scale    <- rep(scale,        length.out = use.n)
  if (length(p)        != use.n)
    p        <- rep(p,            length.out = use.n)



  scase <- abs(shape) < sqrt( .Machine$double.eps )
  nscase <- sum(scase)
  if (use.n - nscase) {
    ans[!scase] <- location[!scase] + scale[!scase] *
        ((1-p[!scase])^(-shape[!scase]) - 1) / shape[!scase]
  }
  if (nscase) {
    ans[scase] <- location[scase] - scale[scase] * log1p(-p[scase])
  }

  ans[p <  0] <- NaN
  ans[p >  1] <- NaN
  ans[(p == 0)] <- location[p == 0]
  ans[(p == 1) & (shape >= 0)] <- Inf
  ind5 <- (p == 1) & (shape < 0)
  ans[ind5] <- location[ind5] - scale[ind5] / shape[ind5]

  ans[scale <= 0] <- NaN
  ans
}







 gpd <- function(threshold = 0,
          lscale = "loge",
          lshape = logoff(offset = 0.5),
          percentiles = c(90, 95),
          iscale = NULL,
          ishape = NULL, 
          tolshape0 = 0.001,
          type.fitted = c("percentiles", "mean"),
          giveWarning = TRUE,
          imethod = 1,
          zero = "shape") {

  type.fitted <- match.arg(type.fitted,
                           c("percentiles", "mean"))[1]

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (!is.logical(giveWarning) || length(giveWarning) != 1)
    stop("bad input for argument 'giveWarning'")
  if (!is.Numeric(threshold)) 
    stop("bad input for argument 'threshold'")
  if (!is.Numeric(imethod, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
     imethod > 2.5)
    stop("argument 'imethod' must be 1 or 2")

  if (length(percentiles) && 
    (!is.Numeric(percentiles, positive = TRUE) ||
     max(percentiles) >= 100))
    stop("bad input for argument 'percentiles'")
  if (!is.Numeric(tolshape0, length.arg = 1, positive = TRUE) ||
      tolshape0 > 0.1)
    stop("bad input for argument 'tolshape0'")



  new("vglmff",
  blurb = c("Generalized Pareto distribution\n",
          "Links:    ",
          namesof("scale", link = lscale, earg = escale ), ", ",
          namesof("shape", link = lshape, earg = eshape )),
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
         parameters.names = c("scale", "shape"),
         lscale    = .lscale ,
         lshape    = .lshape ,
         type.fitted = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero, .type.fitted = type.fitted,
           .lscale = lscale, .lshape = lshape
         ))),



  initialize = eval(substitute(expression({


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
    y.names <- dimnames(y)[[2]]
    if (length(y.names) != ncoly)
      y.names <- paste("Y", 1:ncoly, sep = "")
    extra$y.names <- y.names
    extra$type.fitted <- .type.fitted
    extra$percentiles <- .percentiles



    Threshold <- if (is.Numeric( .threshold )) .threshold else 0
    Threshold <- matrix(Threshold, n, ncoly, byrow = TRUE)
    if (is.Numeric(  .threshold )) {
      orig.y <- y
    }
    ystar <- as.matrix(y - Threshold)  # Operate on ystar
    if (min(ystar, na.rm = TRUE) < 0)
      stop("some response values, after subtracting ",
           "argument 'threshold', are negative. ",
           "Maybe argument 'subset' should be used. ",
           "A threshold value no more than ", min(orig.y, na.rm = TRUE),
           " is needed.")
    extra$threshold <- Threshold


    mynames1 <- param.names("scale", ncoly)
    mynames2 <- param.names("shape", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lscale , earg = .escale , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]



    if (!length(etastart)) {
      meany <- colSums(ystar * w) / colSums(w)
      vary <- apply(ystar, 2, var)
      mediany <- apply(ystar, 2, median)


      init.xii <- if (length( .ishape )) .ishape else {
        if ( .imethod == 1)
          -0.5 * (meany^2 / vary - 1) else
           0.5 * (1 - mediany^2 / vary)
      }
      init.sig <- if (length( .iscale )) .iscale else {
        if (.imethod == 1)
          0.5 * meany * (meany^2 / vary + 1) else
          abs(1 - init.xii) * mediany
      }


      init.xii <- matrix(init.xii, n, ncoly, byrow = TRUE)
      init.sig <- matrix(init.sig, n, ncoly, byrow = TRUE)


      init.sig[init.sig <=  0.0] <-  0.01 # sigma > 0
      init.xii[init.xii <= -0.5] <- -0.40 # Fisher scoring works if xi > -0.5
      init.xii[init.xii >=  1.0] <-  0.90 # Mean/var exists if xi < 1 / 0.5
      if ( .lshape == "loge")
        init.xii[init.xii <= 0.0] <-  0.05



      etastart <-
        cbind(theta2eta(init.sig, .lscale , earg = .escale ),
              theta2eta(init.xii, .lshape , earg = .eshape ))[,
              interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lscale = lscale, .lshape = lshape,
            .iscale = iscale, .ishape = ishape,
            .escale = escale, .eshape = eshape,
            .percentiles = percentiles,
            .threshold = threshold, .type.fitted = type.fitted,
            .imethod = imethod ))),



  linkinv = eval(substitute(function(eta, extra = NULL) {
    sigma <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , earg = .escale )
    shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    if (!is.matrix(sigma))
      sigma <- as.matrix(sigma)
    if (!is.matrix(shape))
      shape <- as.matrix(shape)


    type.fitted <-
      if (length(extra$type.fitted)) {
        extra$type.fitted
      } else {
        warning("cannot find 'type.fitted'. Returning 'percentiles'.")
        "percentiles"
      }

    type.fitted <- match.arg(type.fitted,
                             c("percentiles", "mean"))[1]


    M1 <- 2
    pcent <- extra$percentiles  # Post-20140912


    LP <- length(pcent)  # NULL means LP == 0 and the mean is returned
    ncoly <- ncol(eta) / M1
    if (!length(y.names <- extra$y.names))
      y.names <- paste("Y", 1:ncoly, sep = "")

    Threshold <- extra$threshold



    if (type.fitted == "percentiles" &&  # Upward compatibility:
        LP > 0) {




    do.one <- function(yvec, shape, scale, 
                       threshold,
                       percentiles = c(90, 95),
                       y.name = NULL,
                       tolshape0 = 0.001) {
      is.zero <- (abs(shape) < tolshape0 )  # A matrix

      LP <- length(percentiles)
      fv <- matrix(NA_real_, length(shape), LP)
      is.zero <- (abs(shape) < tolshape0)
      for (ii in 1:LP) {
        temp <- 1 - percentiles[ii] / 100
        fv[!is.zero, ii] <- threshold[!is.zero] +
                           (temp^(-shape[!is.zero]) - 1) *
                           scale[!is.zero] / shape[!is.zero]
        fv[ is.zero, ii] <- threshold[is.zero] - scale[is.zero] * log(temp)
      }

      post.name <- paste(as.character(percentiles), "%", sep = "")

      dimnames(fv) <-
        list(dimnames(shape)[[1]],
             if (is.null(y.name))
               post.name else
               paste(y.name, post.name, sep = " "))
      fv
    }




      fv <- matrix(-1, nrow(sigma),  LP * ncoly)
      colnames.cumsum.fv <- NULL
      for (jlocal in 1:ncoly) {
        block.mat.fv <-
          do.one(yvec = y[, jlocal],
                 shape = shape[, jlocal],
                 scale = sigma[, jlocal],
                 threshold = Threshold[, jlocal],
                 percentiles = pcent,
                 y.name = if (ncoly > 1) y.names[jlocal] else NULL,
                 tolshape0 = .tolshape0 )
        fv[, (jlocal - 1) *  LP + (1:LP)] <- block.mat.fv
        colnames.cumsum.fv <- c(colnames.cumsum.fv,
                                colnames(block.mat.fv))
      }
      colnames(fv) <- colnames.cumsum.fv
    } else {
      fv <- Threshold + sigma / (1 - shape)
      fv[shape >= 1] <- Inf  # Mean exists only if shape < 1.
      dimnames(fv) <- list(dimnames(eta)[[1]], y.names)
    }

    fv
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape,
           .threshold = threshold,
           .tolshape0 = tolshape0 ))),




  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lscale , length = ncoly),
        rep( .lshape , length = ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .escale
      misc$earg[[M1*ii  ]] <- .eshape
    }

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

    misc$true.mu <- FALSE    # @fitted is not a true mu
    misc$percentiles <- .percentiles
    misc$tolshape0 <- .tolshape0
      if (any(Shape < -0.5))
        warning("some values of the shape parameter are less than -0.5")
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .threshold = threshold,
            .tolshape0 = tolshape0, .percentiles = percentiles ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    sigma <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , earg = .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )
    Threshold <- extra$threshold

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- 
        c(w) * dgpd(x = y, location = Threshold, scale = sigma,
                    shape = Shape, tolshape0 = .tolshape0,
                    giveWarning = .giveWarning,
                    log = TRUE, oobounds.log = -1.0e04)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .tolshape0 = tolshape0, .giveWarning= giveWarning,
           .escale = escale, .eshape = eshape,
           .lscale = lscale, .lshape = lshape ))),
  vfamily = c("gpd", "vextremes"),
  deriv = eval(substitute(expression({
    M1 <- 2
    sigma <- eta2theta(eta[, c(TRUE, FALSE)], .lscale , earg = .escale )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    Threshold <- extra$threshold
    ystar <- y - Threshold # Operate on ystar
    A <- 1 + Shape * ystar / sigma

    mytolerance <- .Machine$double.eps
    bad <- (A <= mytolerance)
    if (any(bad) && any(w[bad] != 0)) {
      cat(sum(w[bad],na.rm = TRUE),  # "; ignoring them"
          "observations violating boundary constraints\n")
      flush.console()
    }
    if (any(is.zero <- (abs(Shape) < .tolshape0))) {
    }
    igpd <- !is.zero &  !bad
    iexp <-  is.zero &  !bad

    dl.dShape <- dl.dsigma <- rep(0, length.out = length(y))
    dl.dsigma[igpd] <- ((1 + Shape[igpd]) * ystar[igpd]     / (sigma[igpd] +
                             Shape[igpd]  * ystar[igpd])-1) /  sigma[igpd]

    dl.dShape[igpd] <- log(A[igpd])/Shape[igpd]^2 - (1 + 1/Shape[igpd]) *
                       ystar[igpd] / (A[igpd] * sigma[igpd])
    dl.dShape[iexp] <- ystar[iexp] *
                       (0.5*ystar[iexp]/sigma[iexp] - 1) / sigma[iexp]

    dsigma.deta <- dtheta.deta(sigma, .lscale , earg = .escale )
    dShape.deta <- dtheta.deta(Shape, .lshape , earg = .eshape )

    myderiv <- 
    c(w) * cbind(dl.dsigma * dsigma.deta,
                 dl.dShape * dShape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .tolshape0 = tolshape0,
            .lscale = lscale, .escale = escale,
            .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({


    ned2l.dscale2 <- 1 / ((1+2*Shape) * sigma^2)
    ned2l.dshape2 <- 2 / ((1+2*Shape) * (1+Shape))
    ned2l.dshapescale <- 1 / ((1+2*Shape) * (1+Shape) * sigma)  # > 0 !

    S <- M / M1

    wz <- array(c(c(w) * ned2l.dscale2 * dsigma.deta^2,
                  c(w) * ned2l.dshape2 * dShape.deta^2,
                  c(w) * ned2l.dshapescale * dsigma.deta * dShape.deta),
                dim = c(n, S, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)

    wz
  }), list( .lscale = lscale ))))
}










meplot.default <- function(y, main = "Mean Excess Plot",
    xlab = "Threshold", ylab = "Mean Excess", lty = c(2, 1:2), 
    conf = 0.95, col = c("blue", "black", "blue"), type = "l", ...) {


  if (!is.Numeric(y))
    stop("bad input for argument 'y'")

  n <- length(y)
  sy <- sort(y)
  dsy <- rev(sy)  # decreasing sequence
  me <- rev(cumsum(dsy)) / (n:1) - sy
  me2 <- rev(cumsum(dsy^2))
  var <- (me2 - (n:1) * (me+sy)^2) / (n:1)
  ci <- qnorm((1+conf)/2) * sqrt(abs(var)) / sqrt(n:1)

  ci[length(ci)] <- NA

  mymat <- cbind(me - ci, me, me + ci)
  sy <- sy - sqrt( .Machine$double.eps )

  matplot(sy, mymat, main = main,
          xlab = xlab, ylab = ylab, 
          lty = lty, col = col, type = type, ...)
  invisible(list(threshold = sy,
                 meanExcess = me,
                 plusminus = ci))
}



meplot.vlm <- function(object, ...) {
   if (!length(y <- object@y))
     stop("y slot is empty")
   ans <- meplot(as.numeric(y), ...) 
   invisible(ans)
}



if (!isGeneric("meplot"))
    setGeneric("meplot",
         function(object, ...)
         standardGeneric("meplot"))


setMethod("meplot", "numeric",
         function(object, ...)
         meplot.default(y=object, ...))


setMethod("meplot", "vlm",
         function(object, ...)
         meplot.vlm(object, ...))




guplot.default <-
  function(y, main = "Gumbel Plot",
           xlab = "Reduced data",
           ylab = "Observed data", type = "p", ...) {

    if (!is.Numeric(y))
      stop("bad input for argument 'y'")

    n <- length(y)
    sy <- sort(y)
    x <- -log(-log(((1:n) - 0.5) / n))
    plot(x, sy, main = main, xlab = xlab, ylab = ylab,
         type = type, ...)
    invisible(list(x = x, y = sy))
}



guplot.vlm <- function(object, ...) {
    if (!length(y <- object@y))
      stop("y slot is empty")
    ans <- guplot(as.numeric(y), ...) 
    invisible(ans)
}



if (!isGeneric("guplot"))
    setGeneric("guplot", function(object, ...)
    standardGeneric("guplot"))


setMethod("guplot", "numeric",
         function(object, ...)
         guplot.default(y=object, ...))


setMethod("guplot", "vlm",
         function(object, ...)
         guplot.vlm(object, ...))







 egumbel <- function(llocation = "identitylink",
                     lscale = "loge",
                     iscale = NULL,
                     R = NA, percentiles = c(95, 99),
                     mpv = FALSE, zero = NULL) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (!is.logical(mpv) || length(mpv) != 1)
    stop("bad input for argument 'mpv'")
  if (length(percentiles) &&
     (!is.Numeric(percentiles, positive = TRUE) ||
      max(percentiles) >= 100))
    stop("bad input for argument 'percentiles'")


  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")

  new("vglmff",
  blurb = c("Gumbel distribution (univariate response)\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat, tag = TRUE), ", ",
            namesof("scale",    lscale, earg = escale, tag = TRUE), "\n",
            "Mean:     location + scale*0.5772..\n",
            "Variance: pi^2 * scale^2 / 6"),
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
         llocation = .llocat ,
         lscale    = .lscale ,
         mpv = .mpv ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocation, .lscale = lscale,
           .mpv = mpv ))),



  initialize = eval(substitute(expression({
    y <- cbind(y)
    if (ncol(y) > 1)
      stop("Use gumbel() to handle multivariate responses")
    if (min(y) <= 0)
      stop("all response values must be positive")






    predictors.names <-
      c(namesof("location", .llocat , earg = .elocat , tag = FALSE),
        namesof("scale",    .lscale , earg = .escale , tag = FALSE))


    extra$R <- .R
    extra$mpv <- .mpv
    extra$percentiles <- .percentiles

    if (!length(etastart)) {
      sca.init <-  if (is.Numeric( .iscale, positive = TRUE)) 
                     .iscale else 1.5 * (0.01+sqrt(var(y)*6)) / pi
      sca.init <- rep(sca.init, length.out = n)
      EulerM <- -digamma(1)
      loc.init <- (y - sca.init * EulerM)
      etastart <-
        cbind(theta2eta(loc.init, .llocat , earg = .elocat ),
              theta2eta(sca.init, .lscale , earg = .escale ))
    }
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
                              .iscale = iscale, 
            .R = R, .mpv = mpv, .percentiles = percentiles ))),
  linkinv = eval(substitute( function(eta, extra = NULL) {
    locat <- eta2theta(eta[, 1], .llocat,  earg = .elocat)
    sigma <- eta2theta(eta[, 2], .lscale , earg = .escale )
    EulerM <- -digamma(1)
    Percentiles <- extra$percentiles
    mpv <- extra$mpv
    LP <- length(Percentiles)  # may be 0
    if (!LP) return(locat + sigma * EulerM)
    mu <- matrix(NA_real_, nrow(eta), LP + mpv)
    Rvec <- extra$R
    if (1 <= LP)
    for (ii in 1:LP) {
      ci <- if (is.Numeric(Rvec)) Rvec * (1 - Percentiles[ii] / 100) else
          -log(Percentiles[ii] / 100)
      mu[,ii] <- locat - sigma * log(ci)
    }
    if (mpv)
      mu[, ncol(mu)] <- locat - sigma * log(log(2))
    dmn2 <- if (LP >= 1) paste(as.character(Percentiles), "%",
                               sep = "") else NULL
    if (mpv)
      dmn2 <- c(dmn2, "MPV")
    dimnames(mu) <- list(dimnames(eta)[[1]], dmn2)
    mu
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat,  scale = .lscale) 
    misc$earg <- list(location = .elocat, scale = .escale)
    misc$true.mu <- !length( .percentiles)  # @fitted is not a true mu
    misc$R <- .R
    misc$mpv <- .mpv
    misc$percentiles = .percentiles
  }), list( .llocat = llocat, .lscale = lscale, .mpv = mpv,
            .elocat = elocat, .escale = escale,
            .R = R, .percentiles = percentiles ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    loc <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sca <- eta2theta(eta[, 2], .lscale , earg = .escale )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
                 dgumbel(x = y, location = loc, scale = sca, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale ))),
  vfamily = "egumbel",
  deriv = eval(substitute(expression({
    loc <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sca <- eta2theta(eta[, 2], .lscale , earg = .escale )
    zedd <- (y-loc) / sca
    temp2 <- -expm1(-zedd)
    dl.dloc <- temp2 / sca
    dl.dsca <- -1/sca + temp2 * zedd / sca
    dloc.deta <- dtheta.deta(loc, .llocat , earg = .elocat)
    dsca.deta <- dtheta.deta(sca, .lscale , earg = .escale )
    c(w) * cbind(dl.dloc * dloc.deta,
                 dl.dsca * dsca.deta)
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale ))),
  weight=expression({
    digamma1 <- digamma(1)
    ned2l.dsca2 <- ((2 + digamma1) * digamma1 + trigamma(1) + 1) / sca^2
    ned2l.dloc2 <- 1 / sca^2
    ned2l.dscaloc <- -(1 + digamma1) / sca^2 

    wz = matrix(NA_real_, n, dimm(M = 2))
    wz[, iam(1, 1, M)] <- ned2l.dloc2 * dloc.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dsca2 * dsca.deta^2
    wz[, iam(1, 2, M)] <- ned2l.dscaloc * dloc.deta * dsca.deta

    c(w) * wz
  }))
}




 cens.gumbel <- function(llocation = "identitylink",
                         lscale = "loge",
                         iscale = NULL,
                         mean = TRUE, percentiles = NULL,
                         zero = "scale") {
  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")


  if (!is.logical(mean) || length(mean) != 1)
      stop("mean must be a single logical value")
  if (!mean && (!is.Numeric(percentiles, positive = TRUE) ||
               any(percentiles >= 100)))
    stop("valid percentiles values must be given when mean = FALSE")



  new("vglmff",
  blurb = c("Censored Gumbel distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat, tag = TRUE), ", ", 
            namesof("scale",    lscale, earg = escale, tag = TRUE), "\n",
            "Mean:     location + scale*0.5772..\n",
            "Variance: pi^2 * scale^2 / 6"),
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
         llocation = .llocat ,
         lscale    = .lscale ,
         percentiles = .percentiles ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocation, .lscale = lscale,
           .percentiles = percentiles ))),

  initialize = eval(substitute(expression({
    y <- cbind(y)
    if (ncol(y) > 1)
      stop("Use gumbel.block() to handle multivariate responses")
    if (any(y) <= 0)
      stop("all response values must be positive")



    if (!length(extra$leftcensored))
      extra$leftcensored <- rep(FALSE, length.out = n)
    if (!length(extra$rightcensored))
      extra$rightcensored <- rep(FALSE, length.out = n)
    if (any(extra$rightcensored & extra$leftcensored))
      stop("some observations are both right and left censored!")

    predictors.names <-
    c(namesof("location", .llocat,  earg = .elocat, tag = FALSE),
      namesof("scale",    .lscale ,    earg = .escale   , tag = FALSE))

    if (!length(etastart)) {
      sca.init <-  if (is.Numeric( .iscale, positive = TRUE)) 
                      .iscale else 1.1 * sqrt(var(y) * 6 ) / pi
      sca.init <- rep(sca.init, length.out = n)
      EulerM <- -digamma(1)
      loc.init <- (y - sca.init * EulerM)
      loc.init[loc.init <= 0] = min(y)
      etastart <-
        cbind(theta2eta(loc.init, .llocat , earg = .elocat ),
              theta2eta(sca.init, .lscale , earg = .escale ))
    }
  }), list( .lscale = lscale, .iscale = iscale,
            .llocat = llocat, 
            .elocat = elocat, .escale = escale ))), 
  linkinv = eval(substitute( function(eta, extra = NULL) {
    loc  <- eta2theta(eta[, 1], .llocat)
    sc   <- eta2theta(eta[, 2], .lscale)
    EulerM <- -digamma(1)
    if (.mean) loc + sc * EulerM else {
      LP <- length(.percentiles)  # 0 if NULL
      mu <- matrix(NA_real_, nrow(eta), LP)
      for (ii in 1:LP) {
          ci <- -log( .percentiles[ii] / 100)
          mu[, ii] <- loc - sc * log(ci)
      }
      dmn2 <- paste(as.character(.percentiles), "%", sep = "")
      dimnames(mu) <- list(dimnames(eta)[[1]], dmn2)
      mu
    }
  }, list( .lscale = lscale, .percentiles = percentiles,
           .llocat = llocat, 
           .elocat = elocat, .escale = escale ,
           .mean=mean ))), 
  last = eval(substitute(expression({
        misc$link <- c(location= .llocat,  scale = .lscale) 
        misc$earg <- list(location= .elocat, scale= .escale )
        misc$true.mu <- .mean    # if FALSE then @fitted is not a true mu 
        misc$percentiles = .percentiles
  }), list( .lscale = lscale, .mean=mean,
            .llocat = llocat, 
            .elocat = elocat, .escale = escale ,
            .percentiles = percentiles ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    loc <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    sc  <- eta2theta(eta[, 2], .lscale , earg = .escale )
    zedd <- (y-loc) / sc

    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cen0 <- !cenL & !cenU   # uncensored obsns
    Fy <- exp(-exp(-zedd))
    ell1 <- -log(sc[cen0]) - zedd[cen0] - exp(-zedd[cen0])
    ell2 <- log(Fy[cenL])
    ell3 <- log1p(-Fy[cenU])
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(w[cen0] * ell1) + sum(w[cenL] * ell2) + sum(w[cenU] * ell3)
  }, list( .lscale = lscale,
           .llocat = llocat, 
           .elocat = elocat, .escale = escale ))),
  vfamily = "cens.gumbel",
  deriv = eval(substitute(expression({
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cen0 <- !cenL & !cenU   # uncensored obsns

    loc <- eta2theta(eta[, 1], .llocat,  earg = .elocat )
    sc  <- eta2theta(eta[, 2], .lscale , earg = .escale )
    zedd <- (y-loc) / sc
    temp2 <- -expm1(-zedd)
    dl.dloc <- temp2 / sc
    dl.dsc <- -1/sc + temp2 * zedd / sc
    dloc.deta <- dtheta.deta(loc, .llocat,  earg = .elocat )
    dsc.deta <- dtheta.deta(sc, .lscale , earg = .escale )

    ezedd <- exp(-zedd)
    Fy <- exp(-ezedd)
    dFy.dloc <- -ezedd * Fy / sc
    dFy.dsc <- zedd * dFy.dloc # -zedd * exp(-zedd) * Fy / sc
    if (any(cenL)) {
      dl.dloc[cenL] <- -ezedd[cenL] / sc[cenL]
      dl.dsc[cenL] <- -zedd[cenL] * ezedd[cenL] / sc[cenL]
    }
    if (any(cenU)) {
      dl.dloc[cenU] <- -dFy.dloc[cenU] / (1-Fy[cenU])
      dl.dsc[cenU] <- -dFy.dsc[cenU] / (1-Fy[cenU])
    }
    c(w) * cbind(dl.dloc * dloc.deta,
                 dl.dsc * dsc.deta)
  }), list( .lscale = lscale,
            .llocat = llocat, 
            .elocat = elocat, .escale = escale ))),
  weight = expression({
    A1 <- ifelse(cenL, Fy, 0)
    A3 <- ifelse(cenU, 1-Fy, 0)
    A2 <- 1 - A1 - A3   # Middle; uncensored
    digamma1 <- digamma(1)
    ed2l.dsc2 <- ((2+digamma1)*digamma1 + trigamma(1) + 1) / sc^2
    ed2l.dloc2 <- 1 / sc^2
    ed2l.dlocsc <- -(1 + digamma1) / sc^2 
    wz <- matrix(NA_real_, n, dimm(M = 2))
    wz[, iam(1, 1, M)] <- A2 * ed2l.dloc2 * dloc.deta^2
    wz[, iam(2, 2, M)] <- A2 * ed2l.dsc2 * dsc.deta^2
    wz[, iam(1, 2, M)] <- A2 * ed2l.dlocsc * dloc.deta * dsc.deta
    d2l.dloc2 <- -ezedd / sc^2
    d2l.dsc2 <- (2 - zedd) * zedd * ezedd / sc^2
    d2l.dlocsc <- (1 - zedd) * ezedd / sc^2
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)]-A1^2 * d2l.dloc2 * dloc.deta^2
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)]-A1^2 * d2l.dsc2 * dsc.deta^2
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)]-A1^2 * d2l.dlocsc *
                        dloc.deta * dsc.deta
    d2Fy.dloc2 <- dFy.dloc * dl.dloc + Fy * d2l.dloc2
    d2Fy.dsc2 <- dFy.dsc * dl.dsc + Fy * d2l.dsc2
    d2Fy.dlocsc <- dFy.dsc * dl.dloc + Fy * d2l.dlocsc
    d2l.dloc2 <- -((1-Fy) * d2Fy.dloc2 - dFy.dloc^2) / (1-Fy)^2
    d2l.dsc2 <- -((1-Fy) * d2Fy.dsc2 - dFy.dsc^2) / (1-Fy)^2
    d2l.dlocsc  <- -((1-Fy) * d2Fy.dlocsc - dFy.dloc * dFy.dsc) / (1-Fy)^2
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)]-A3^2 * d2l.dloc2 * dloc.deta^2
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)]-A3^2 * d2l.dsc2 * dsc.deta^2
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)]-A3^2 * d2l.dlocsc *
                          dloc.deta * dsc.deta
    c(w) * wz
  }))
}




dfrechet <- function(x, location = 0, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(scale), length(shape), length(location))
  x        <- rep(x,        length.out = L)
  scale    <- rep(scale,    length.out = L)
  shape    <- rep(shape,    length.out = L)
  location <- rep(location, length.out = L)

  logdensity <- rep(log(0), length.out = L)
  xok <- (x > location)
  rzedd <- scale / (x - location)
  logdensity[xok] <- log(shape[xok]) - (rzedd[xok]^shape[xok]) +
                    (shape[xok]+1) * log(rzedd[xok]) - log(scale[xok])
  logdensity[shape <= 0] <- NaN
  logdensity[scale <= 0] <- NaN

  if (log.arg) logdensity else exp(logdensity)
}



pfrechet <- function(q, location = 0, scale = 1, shape,
                     lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  rzedd <- scale / (q - location)


  if (lower.tail) {
    if (log.p) {
      ans <- -(rzedd^shape)
      ans[q <= location] <- -Inf
    } else {  
      ans <- exp(-(rzedd^shape))
      ans[q <= location] <- 0 
      }
  } else {
    if (log.p) {
      ans <- log(-expm1(-(rzedd^shape)))
      ans[q <= location] <- 0
    } else {
      ans <- -expm1(-(rzedd^shape)) 
      ans[q <= location]  <- 1
    }
  }
  ans
}



qfrechet <- function(p, location = 0, scale = 1, shape,
                     lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- location + scale * (-ln.p)^(-1 / shape)
      ans[ln.p > 0] <- NaN
    } else {
      ans <- location + scale * (-log(p))^(-1 / shape)
      ans[p < 0] <- NaN
      ans[p == 0] <- location
      ans[p == 1] <- Inf
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- location + scale * (-log(-expm1(ln.p)))^(-1 / shape)
      ans[ln.p > 0] <- NaN
    } else { 
      ans <- location + scale * (-log1p(-p))^(-1 / shape)
      ans[p < 0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- location
      ans[p > 1] <- NaN
    }
  }
  ans
}



rfrechet <- function(n, location = 0, scale = 1, shape) {
  if (!is.Numeric(scale, positive = TRUE))
    stop("scale must be positive")
  if (!is.Numeric(shape, positive = TRUE))
    stop("shape must be positive")

  location + scale * (-log(runif(n)))^(-1/shape)
}








frechet.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}



 frechet <- function(location = 0,
                     lscale = "loge",
                     lshape = logoff(offset = -2),
                     iscale = NULL, ishape = NULL,
                     nsimEIM = 250,
                     zero = NULL) {

  if (!is.Numeric(location))
    stop("bad input for argument 'location'")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")



  stopifnot(nsimEIM > 10, length(nsimEIM) == 1, nsimEIM == round(nsimEIM))


  new("vglmff",
  blurb = c("2-parameter Frechet distribution\n",
            "Links:    ",
            namesof("scale", link = lscale, earg = escale ), ", ",
            namesof("shape", link = lshape, earg = eshape )),
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
         parameters.names = c("scale", "shape"),
         lscale  = .lscale ,
         lshape  = .lshape ,
         nsimEIM = .nsimEIM ,
         zero = .zero )
  }, list( .zero = zero,
           .lscale = lscale,
           .lshape = lshape,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({


      temp5 <-
      w.y.check(w = w, y = y,
                ncol.w.max = 1,
                ncol.y.max = 1,
                out.wy = TRUE,
                colsyperw = 1,
                maximize = TRUE)
      w <- temp5$w
      y <- temp5$y



    predictors.names <-
      c(namesof("scale", .lscale , earg = .escale, short = TRUE),
        namesof("shape", .lshape , earg = .eshape, short = TRUE))


    extra$location <- rep( .location , length.out = n)  # stored here


    if (!length(etastart)) {
      locinit = extra$location
      if (any(y <= locinit))
        stop("initial values for 'location' are out of range")


      frech.aux <- function(shapeval, y, x, w, extraargs) {
        myprobs <- c(0.25, 0.5, 0.75)
        myobsns <- quantile(y, probs = myprobs)
        myquant <- (-log(myprobs))^(-1/shapeval)
        myfit <- lsfit(x = myquant, y = myobsns, intercept = TRUE)
        sum(myfit$resid^2)
      } 

      shape.grid <- c(100, 70, 40, 20, 12, 8, 4, 2, 1.5)
      shape.grid <- c(1 / shape.grid, 1, shape.grid)
      try.this <- grid.search(shape.grid, objfun = frech.aux,
                              y = y,  x = x, w = w, maximize = FALSE,
                              abs.arg = TRUE)

      shape.init <- if (length( .ishape ))
        rep( .ishape , length.out = n) else {
        rep(try.this , length.out = n)  # variance exists if shape > 2
      }


      myprobs <- c(0.25, 0.5, 0.75)
      myobsns <- quantile(y, probs = myprobs)
      myquant <- (-log(myprobs))^(-1/shape.init[1])
      myfit <- lsfit(x = myquant, y = myobsns)

    Scale.init <- if (length( .iscale ))
                 rep( .iscale , length.out = n) else {
      if (all(shape.init > 1)) {
        myfit$coef[2]
      } else {
        rep(1.0, length.out = n)
      }
    }

    etastart <-
      cbind(theta2eta(Scale.init, .lscale , earg = .escale ),
            theta2eta(shape.init, .lshape , earg = .eshape ))
    }
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .iscale = iscale, .ishape = ishape,
            .location = location ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    loc <- extra$location
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    ans <- rep(NA_real_, length.out = length(shape))
    ok <- shape > 1
    ans[ok] <- loc[ok] + Scale[ok] * gamma(1 - 1/shape[ok])
    ans
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),
  last = eval(substitute(expression({
    misc$links <-   c("scale" = .lscale , "shape" = .lshape )

    misc$earg <- list("scale" = .escale , "shape" = .eshape )

    misc$nsimEIM <- .nsimEIM
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    loctn <- extra$location
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
                 dfrechet(x = y, location = loctn, scale = Scale,
                          shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),
  vfamily = c("frechet", "vextremes"),
  deriv = eval(substitute(expression({
    loctn <- extra$location
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    rzedd <- Scale / (y - loctn)  # reciprocial of zedd
    dl.dloctn <- (shape + 1) / (y - loctn) -
                (shape / (y - loctn)) * (rzedd)^shape
    dl.dScale <- shape * (1 - rzedd^shape) / Scale
    dl.dshape <- 1 / shape + log(rzedd) * (1 -  rzedd^shape)

    dthetas.detas <- cbind(
      dScale.deta <- dtheta.deta(Scale, .lscale , earg = .escale ),
      dShape.deta <- dtheta.deta(shape, .lshape , earg = .eshape ))

    c(w) * cbind(dl.dScale,
                 dl.dshape) * dthetas.detas
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape ))),
  weight = eval(substitute(expression({

    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

    if (length( .nsimEIM )) {
      for (ii in 1:( .nsimEIM )) {
        ysim <- rfrechet(n, location = loctn, scale = Scale, shape = shape)

          rzedd <- Scale / (ysim - loctn)  # reciprocial of zedd
          dl.dloctn <- (shape + 1) / (ysim - loctn) -
                      (shape / (ysim - loctn)) * (rzedd)^shape
          dl.dScale <- shape * (1 - rzedd^shape) / Scale
          dl.dshape <- 1 / shape + log(rzedd) * (1 -  rzedd^shape)

          rm(ysim)
          temp3 <- cbind(dl.dScale, dl.dshape)
          run.varcov <- run.varcov +
                       temp3[, ind1$row.index] *
                       temp3[, ind1$col.index]
      }
      run.varcov <- run.varcov / .nsimEIM

      wz = if (intercept.only)
          matrix(colMeans(run.varcov),
                 n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz = c(w) * wz * dthetas.detas[, ind1$row.index] *
                       dthetas.detas[, ind1$col.index]
    } else {
      stop("argument 'nsimEIM' must be numeric")
    }

    wz
  }), list( .nsimEIM = nsimEIM ))))
}








rec.normal.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}


 rec.normal <- function(lmean = "identitylink", lsd = "loge",
                        imean = NULL, isd = NULL, imethod = 1,
                        zero = NULL) {
  lmean <- as.list(substitute(lmean))
  emean <- link2list(lmean)
  lmean <- attr(emean, "function.name")

  lsdev <- as.list(substitute(lsd))
  esdev <- link2list(lsdev)
  lsdev <- attr(esdev, "function.name")

  isdev <- isd


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3.5)
    stop("argument 'imethod' must be 1 or 2 or 3")



  new("vglmff",
  blurb = c("Upper record values from a univariate normal distribution\n\n",
            "Links:    ",
            namesof("mean", lmean, emean, tag = TRUE), "; ",
            namesof("sd",   lsdev, esdev, tag = TRUE),
            "\n",
            "Variance: sd^2"),
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
         parameters.names = c("mean", "sd"),
         lmean  = .lmean ,
         lsd    = .lsd ,
         imethod = .imethod ,
         zero = .zero )
  }, list( .zero = zero,
           .lmean = lmean,
           .lsd   = lsd,
           .imethod = imethod ))),



  initialize = eval(substitute(expression({



    predictors.names <-
      c(namesof("mean", .lmean, .emean, tag = FALSE),
        namesof("sd",   .lsdev, .esdev, tag = FALSE))

    if (ncol(y <- cbind(y)) != 1)
        stop("response must be a vector or a one-column matrix")

    if (any(diff(y) <= 0))
        stop("response must have increasingly larger and larger values")
    if (any(w != 1))
        warning("weights should have unit values only")


    if (!length(etastart)) {
        mean.init <- if (length( .imean )) rep( .imean ,
                                               length.out = n) else {
            if (.lmean == "loge") pmax(1/1024, min(y)) else min(y)}
        sd.init <- if (length( .isdev)) rep( .isdev , length.out = n) else {
            if (.imethod == 1)  1*(sd(c(y))) else
            if (.imethod == 2)  5*(sd(c(y))) else
                                  .5*(sd(c(y)))
            }
        etastart <-
          cbind(theta2eta(rep(mean.init, len = n), .lmean, .emean ),
                theta2eta(rep(sd.init,   len = n), .lsdev, .esdev ))
    }
  }), list( .lmean = lmean, .lsdev = lsdev,
            .emean = emean, .esdev = esdev,
            .imean = imean, .isdev = isdev,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmean, .emean )
  }, list( .lmean = lmean, .emean = emean ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmean , "sd" = .lsdev )
    misc$earg <- list("mu" = .emean , "sd" = .esdev )


    misc$expected = FALSE
  }), list( .lmean = lmean, .lsdev = lsdev,
            .emean = emean, .esdev = esdev ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    sdev <- eta2theta(eta[, 2], .lsdev )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      zedd <- (y - mu) / sdev
      NN <- nrow(eta)
      if (summation) {
        sum(w * (-log(sdev) - 0.5 * zedd^2)) -
        sum(w[-NN] * pnorm(zedd[-NN], lower.tail = FALSE, log.p = TRUE))
      } else {
        stop("cannot handle 'summation = FALSE' yet")
      }
    }
  }, list( .lsdev = lsdev, .esdev = esdev ))),
  vfamily = c("rec.normal"),
  deriv = eval(substitute(expression({
    NN <- nrow(eta)
    mymu <- eta2theta(eta[, 1], .lmean)
    sdev <- eta2theta(eta[, 2], .lsdev)
    zedd <- (y - mymu) / sdev
    temp200 <- dnorm(zedd) / (1-pnorm(zedd))
    dl.dmu <- (zedd - temp200) / sdev
    dl.dmu[NN] <- zedd[NN] / sdev[NN]
    dl.dsd <- (-1 + zedd^2 - zedd * temp200)  / sdev
    dl.dsd[NN] <- (-1 + zedd[NN]^2)  / sdev[NN]

    dmu.deta <- dtheta.deta(mymu, .lmean, .emean )
    dsd.deta <- dtheta.deta(sdev, .lsdev, .esdev )

    if (iter == 1) {
      etanew <- eta
    } else {
      derivold <- derivnew
      etaold <- etanew
      etanew <- eta
    }
    derivnew <- c(w) * cbind(dl.dmu * dmu.deta,
                            dl.dsd * dsd.deta)
    derivnew
    }), list( .lmean = lmean, .lsdev = lsdev,
              .emean = emean, .esdev = esdev ))),
    weight = expression({
      if (iter == 1) {
          wznew <- cbind(matrix(w, n, M), matrix(0, n, dimm(M)-M))
      } else {
        wzold <- wznew
        wznew <- qnupdate(w = w, wzold = wzold,
                          dderiv = (derivold - derivnew),
                          deta = etanew-etaold, M = M,
                          trace = trace)  # weights incorporated in args
    }
    wznew
  }))
}



rec.exp1.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 rec.exp1 <- function(lrate = "loge", irate = NULL, imethod = 1) {
  lrate <- as.list(substitute(lrate))
  erate <- link2list(lrate)
  lrate <- attr(erate, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                    integer.valued = TRUE, positive = TRUE) ||
      imethod > 3.5)
    stop("argument 'imethod' must be 1 or 2 or 3")



  new("vglmff",
  blurb = c("Upper record values from a ",
            "1-parameter exponential distribution\n\n",
            "Links:    ",
            namesof("rate", lrate, erate, tag = TRUE),
            "\n",
            "Variance: 1/rate^2"),
  initialize = eval(substitute(expression({
    predictors.names <-
      c(namesof("rate", .lrate , .erate , tag = FALSE))

    if (ncol(y <- cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")
    if (any(diff(y) <= 0))
      stop("response must have increasingly larger and larger values")
    if (any(w != 1))
      warning("weights should have unit values only")


    if (!length(etastart)) {
      rate.init <- if (length( .irate ))
                  rep( .irate , len = n) else {
          init.rate <-
              if (.imethod == 1) length(y) / y[length(y), 1] else
              if (.imethod == 2) 1/mean(y) else 1/median(y)
          if (.lrate == "loge") pmax(1/1024, init.rate) else
            init.rate}

      etastart <-
        cbind(theta2eta(rep(rate.init, len = n), .lrate , .erate ))
      }
  }), list( .lrate = lrate,
            .erate = erate,
            .irate = irate, .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .lrate , .erate )
  }, list( .lrate = lrate, .erate = erate ))),
  last = eval(substitute(expression({
    misc$link <-    c("rate" = .lrate)
    misc$earg <- list("rate" = .erate)

    misc$expected = TRUE
  }), list( .lrate = lrate, .erate = erate ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    rate <- eta2theta(eta, .lrate , .erate )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      NN <- length(eta)
      y <- cbind(y)
      if (summation) {
        sum(w * log(rate)) - w[NN] * rate[NN] * y[NN, 1]
      } else {
        stop("cannot handle 'summation = FALSE' yet")
      }
    }
  }, list( .lrate = lrate, .erate = erate ))),
  vfamily = c("rec.exp1"),
  deriv = eval(substitute(expression({
    NN <- length(eta)
    rate <- c(eta2theta(eta, .lrate , .erate ))

    dl.drate <- 1 / rate 
    dl.drate[NN] <- 1/ rate[NN] - y[NN, 1]

    drate.deta <- dtheta.deta(rate, .lrate , .erate )

    c(w) * cbind(dl.drate * drate.deta)
  }), list( .lrate = lrate, .erate = erate ))),
  weight = expression({
    ed2l.drate2 <- 1 / rate^2
    wz <- drate.deta^2 * ed2l.drate2
    c(w) * wz
  }))
}









dpois.points <- function(x, lambda, ostatistic,
                         dimension = 2, log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  L <- max(length(x), length(lambda),
           length(ostatistic), length(dimension))
  if (length(x) != L)
    x          <- rep(x,          length.out = L)
  if (length(lambda) != L)
    lambda     <- rep(lambda,     length.out = L)
  if (length(ostatistic) != L)
    ostatistic <- rep(ostatistic, length.out = L)
  if (length(dimension) != L)
    dimension  <- rep(dimension,  length.out = L)

  if (!all(dimension %in% c(2, 3)))
    stop("argument 'dimension' must have values 2 and/or 3")


  ans2 <- log(2) + ostatistic * log(pi * lambda) -
          lgamma(ostatistic) + (2 * ostatistic - 1) * log(x) -
          lambda * pi * x^2
  ans2[x < 0 | is.infinite(x)] <- log(0)  # 20141209 KaiH

  ans3 <- log(3) + ostatistic * log(4 * pi * lambda / 3) -
          lgamma(ostatistic) + (3 * ostatistic - 1) * log(x) -
          (4/3) * lambda * pi * x^3
  ans3[x < 0 | is.infinite(x)] <- log(0)  # 20141209 KaiH

  ans <- ifelse(dimension == 2, ans2, ans3)


  if (log.arg) ans else exp(ans)
}



 poisson.points <-
  function(ostatistic, dimension = 2,
           link = "loge",
           idensity = NULL, imethod = 1) {


  if (!is.Numeric(ostatistic,
                  length.arg = 1,
                  positive = TRUE))
    stop("argument 'ostatistic' must be a single positive integer")
  if (!is.Numeric(dimension, positive = TRUE,
                  length.arg = 1, integer.valued = TRUE) ||
      dimension > 3)
    stop("argument 'dimension' must be 2 or 3")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE) ||
      imethod > 2.5)
    stop("argument 'imethod' must be 1 or 2")
  if (length(idensity) &&
      !is.Numeric(idensity, positive = TRUE))
    stop("bad input for argument 'idensity'")

  new("vglmff",
  blurb = c(if (dimension == 2)
          "Poisson-points-on-a-plane distances distribution\n" else
          "Poisson-points-on-a-volume distances distribution\n",
          "Link:    ",
          namesof("density", link, earg = earg), "\n\n",
          if (dimension == 2)
            "Mean:    gamma(s+0.5) / (gamma(s) * sqrt(density * pi))" else
            "Mean:    gamma(s+1/3) / (gamma(s) * (4*density*pi/3)^(1/3))"),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")
    if (any(y <= 0))
      stop("response must contain positive values only")



    predictors.names <-
      namesof("density", .link, earg = .earg, tag = FALSE) 



    if (!length(etastart)) {
      use.this <- if ( .imethod == 1) median(y) + 1/8 else
                  weighted.mean(y,w)
      if ( .dimension == 2) {
        myratio <- exp(lgamma( .ostatistic + 0.5) -
                       lgamma( .ostatistic ))
        density.init <- if (is.Numeric( .idensity ))
            rep( .idensity , len = n) else
            rep(myratio^2 / (pi * use.this^2), len = n)
        etastart <- theta2eta(density.init, .link , earg = .earg )
      } else {
        myratio <- exp(lgamma( .ostatistic + 1/3) -
                       lgamma( .ostatistic ))
        density.init <- if (is.Numeric( .idensity ))
            rep( .idensity , len = n) else
            rep(3 * myratio^3 / (4 * pi * use.this^3), len = n)
        etastart <- theta2eta(density.init, .link , earg = .earg )
      }
    }
  }), list( .link = link, .earg = earg,
            .ostatistic = ostatistic,
            .dimension = dimension, .imethod = imethod,
            .idensity = idensity ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    density <- eta2theta(eta, .link, earg = .earg)
    if ( .dimension == 2) {
      myratio <- exp(lgamma( .ostatistic + 0.5) - lgamma( .ostatistic ))
      myratio / sqrt(density * pi)
    } else {
      myratio <- exp(lgamma( .ostatistic + 1/3) - lgamma( .ostatistic))
      myratio / (4 * density * pi/3)^(1/3)
    }
  }, list( .link = link, .earg = earg,
           .ostatistic = ostatistic,
           .dimension = dimension ))),
  last = eval(substitute(expression({
    misc$link <-    c("density" = .link)
    misc$earg <- list("density" = .earg)

    misc$expected <- TRUE
    misc$ostatistic <- .ostatistic
    misc$dimension <- .dimension
  }), list( .link = link, .earg = earg,
            .ostatistic = ostatistic,
            .dimension = dimension ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    density <- eta2theta(eta, .link, earg = .earg)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dpois.points(y, lambda = density,
                                     ostatistic = .ostatistic ,
                                     dimension = .dimension , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg,
           .ostatistic = ostatistic,
           .dimension = dimension ))),
  vfamily = c("poisson.points"),
  deriv = eval(substitute(expression({
    density <- eta2theta(eta, .link, earg = .earg)

    dl.ddensity <- if ( .dimension == 2) {
      .ostatistic / density - pi * y^2
    } else {
      .ostatistic / density - (4/3) * pi * y^3
    }

    ddensity.deta <- dtheta.deta(density, .link , earg = .earg )

    c(w) * dl.ddensity * ddensity.deta
  }), list( .link = link, .earg = earg,
            .ostatistic = ostatistic,
            .dimension = dimension ))),
  weight = eval(substitute(expression({
    ned2l.ddensity2 <- .ostatistic / density^2
    wz <- ddensity.deta^2 * ned2l.ddensity2
    c(w) * wz
  }), list( .link = link, .earg = earg,
            .ostatistic = ostatistic,
            .dimension = dimension ))))
}





