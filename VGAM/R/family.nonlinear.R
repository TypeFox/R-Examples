# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.











vnonlinear.control <- function(save.weights = TRUE, ...) {



    list(save.weights = as.logical(save.weights)[1])
}




subset.lohi <- function(xvec, yvec,
                        probs.x = c(0.15, 0.85),
                        type = c("median", "wtmean", "unwtmean"),
                        wtvec = rep(1, len = length(xvec))) {


  if (!is.Numeric(probs.x, length.arg = 2))
    stop("argument 'probs.x' must be numeric and of length two")

  min.q <- quantile(xvec, probs = probs.x[1] )
  max.q <- quantile(xvec, probs = probs.x[2] )

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("median", "wtmean", "unwtmean"))[1]


  if (type == "median") {
    y1bar <- median(yvec[xvec < min.q])
    x1bar <- median(xvec[xvec < min.q])
    y2bar <- median(yvec[xvec > max.q])
    x2bar <- median(xvec[xvec > max.q])
  }
  if (type == "wtmean") {
    y1bar <- weighted.mean(yvec[xvec < min.q], w = wtvec[xvec < min.q])
    x1bar <- weighted.mean(xvec[xvec < min.q], w = wtvec[xvec < min.q])
    y2bar <- weighted.mean(yvec[xvec > max.q], w = wtvec[xvec > max.q])
    x2bar <- weighted.mean(xvec[xvec > max.q], w = wtvec[xvec > max.q])
  }
  if (type == "unwtmean") {
    y1bar <- mean(yvec[xvec < min.q])
    x1bar <- mean(xvec[xvec < min.q])
    y2bar <- mean(yvec[xvec > max.q])
    x2bar <- mean(xvec[xvec > max.q])
  }

  if (x1bar >= x2bar)
    stop("cannot find two distinct x values; try decreasing the first ",
         "value of argument 'probs.x' and increasing the second value")

  list(x1bar =  x1bar,
       y1bar =  y1bar,
       x2bar =  x2bar,
       y2bar =  y2bar,
       slopeUp = (y2bar > y1bar))
}



micmen.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}





 micmen <- function(rpar = 0.001, divisor = 10,
                    init1 = NULL, init2 = NULL,
                    imethod = 1,
                    oim = TRUE,
                    link1 = "identitylink", link2 = "identitylink",
                    firstDeriv = c("nsimEIM", "rpar"),
                    probs.x = c(0.15, 0.85),
                    nsimEIM = 500,
                    dispersion = 0, zero = NULL) {



  firstDeriv <- match.arg(firstDeriv, c("nsimEIM", "rpar"))[1]

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("argument 'imethod' must be integer")
  if (!is.Numeric(probs.x, length.arg = 2))
    stop("argument 'probs.x' must be numeric and of length two")
  if (!is.logical(oim) || length(oim) != 1)
    stop("argument 'oim' must be single logical")

    stopifnot(nsimEIM > 10, length(nsimEIM) == 1,
              nsimEIM == round(nsimEIM))

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("'imethod' must be 1 or 2 or 3")


  estimated.dispersion <- (dispersion == 0)

  link1 <- as.list(substitute(link1))
  earg1 <- link2list(link1)
  link1 <- attr(earg1, "function.name")

  link2 <- as.list(substitute(link2))
  earg2 <- link2list(link2)
  link2 <- attr(earg2, "function.name")


  new("vglmff",
  blurb = c("Michaelis-Menton regression model\n",
            "Y_i = theta1 * u_i / (theta2 + u_i) + e_i\n\n",
            "Links:    ",
            namesof("theta1", link1, earg = earg1), ", ",
            namesof("theta2", link2, earg = earg2),
            "\n",
            "Variance: constant"),

  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero))),


  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M <- if (is.matrix(y)) ncol(y) else 1
    if (residuals) {
      if (M > 1) NULL else (y - mu) * sqrt(w)
    } else {
      ResSS.vgam(y - mu, w, M = M)
    }
  },


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("theta1", "theta2"),
         link1    = .link1 ,
         link2    = .link2 ,
         zero = .zero )
  }, list( .zero = zero,
           .link1 = link1, .link2 = link2
         ))),

  initialize = eval(substitute(expression({


    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w



    if (!length(Xm2))
      stop("regressor not found")
    if (ncol(as.matrix(Xm2)) != 1)
      stop("regressor not found or is not a vector. Use the ",
           "'form2' argument without an intercept")
    Xm2 <- as.vector(Xm2)  # Make sure
    extra$Xm2 <- Xm2          # Needed for @linkinv

    predictors.names <-
      c(namesof("theta1", .link1 , earg = .earg1, tag = FALSE),
        namesof("theta2", .link2 , earg = .earg2, tag = FALSE))

    if (length(mustart) || length(coefstart))
      stop("cannot handle 'mustart' or 'coefstart'")

    if (!length(etastart)) {
      if ( .imethod == 3 ) {
        index0 <- (1:n)[Xm2 <= quantile(Xm2, prob = .probs.x[2] )]
        init1 <- median(y[index0])
        init2 <- median(init1 * Xm2 / y - Xm2)
      }
      if ( .imethod == 1 || .imethod == 2) {
        mysubset <- subset.lohi(Xm2, y, probs.x = .probs.x,
                  type = ifelse( .imethod == 1, "median", "wtmean"),
                  wtvec = w)

        mat.x <- with(mysubset, cbind(c(x1bar, x2bar), -c(y1bar, y2bar)))
        theta.temp <- with(mysubset,
                           solve(mat.x, c(x1bar * y1bar, x2bar * y2bar)))
        init1 <- theta.temp[1]
        init2 <- theta.temp[2]



      }


      if (length( .init1 )) init1 <- .init1
      if (length( .init2 )) init2 <- .init2

      etastart <- cbind(
          rep(theta2eta(init1, .link1 , earg = .earg1 ), len = n),
          rep(theta2eta(init2, .link2 , earg = .earg2 ), len = n))
    } else {
      stop("cannot handle 'etastart' or 'mustart'")
    }
  }), list( .init1 = init1, .link1 = link1, .earg1 = earg1,
            .init2 = init2, .link2 = link2, .earg2 = earg2,
            .imethod = imethod,
            .probs.x = probs.x ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta1 <- eta2theta(eta[, 1], .link1 , earg = .earg1 )
    theta2 <- eta2theta(eta[, 2], .link2 , earg = .earg2 )
    theta1 * extra$Xm2 / (theta2 + extra$Xm2)
  }, list( .link1 = link1, .earg1 = earg1,
           .link2 = link2, .earg2 = earg2))),

  last = eval(substitute(expression({
    misc$link <-    c(theta1 = .link1 , theta2 = .link2 )

    misc$earg <- list(theta1 = .earg1 , theta2 = .earg2 )

    misc$rpar <- rpar
    fit$df.residual <- n - rank   # Not nrow.X.vlm - rank
    fit$df.total <- n             # Not nrow.X.vlm

    extra$Xm2 <- NULL             # Regressor is in control$regressor 
    dpar <- .dispersion
    if (!dpar) {
      dpar <- sum(c(w) * (y - mu)^2) / (n - ncol.X.vlm)
    }
    misc$dispersion <- dpar

    misc$default.dispersion <- 0
    misc$estimated.dispersion <- .estimated.dispersion

    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
    misc$firstDeriv <- .firstDeriv
    misc$oim <- .oim
    misc$rpar <- rpar
    misc$orig.rpar <- .rpar
    misc$multipleResponses <- FALSE
  }), list( .link1 = link1, .earg1 = earg1,
            .link2 = link2, .earg2 = earg2,
            .dispersion = dispersion,
            .imethod = imethod,
            .firstDeriv = firstDeriv,
            .oim = oim, .rpar = rpar,
            .nsimEIM = nsimEIM,
            .estimated.dispersion = estimated.dispersion ))),

  summary.dispersion = FALSE,

  vfamily = c("micmen", "vnonlinear"),

  deriv = eval(substitute(expression({
    theta1 <- eta2theta(eta[, 1], .link1 , earg = .earg1 )
    theta2 <- eta2theta(eta[, 2], .link2 , earg = .earg2 )
    dthetas.detas <- cbind(dtheta.deta(theta1, .link1 , earg = .earg1 ),
                           dtheta.deta(theta2, .link2 , earg = .earg2 ))

    rpar <- if ( .firstDeriv == "rpar") {
      if (iter > 1) {
        max(rpar / .divisor, 1000 * .Machine$double.eps)
      } else {
        d3 <- deriv3(~ theta1 * Xm2 / (theta2 + Xm2),
                     c("theta1", "theta2"), hessian = FALSE)
        .rpar
      }
    } else {
      .rpar
    }

    dmus.dthetas <- if (FALSE) {
      attr(eval(d3), "gradient")
    } else {
      dmu.dtheta1 <-           Xm2 / (theta2 + Xm2)
      dmu.dtheta2 <- -theta1 * Xm2 / (Xm2 + theta2)^2
      cbind(dmu.dtheta1, dmu.dtheta2)
    }

    myderiv <- if ( .firstDeriv == "rpar") {
      if (TRUE) {
        index <- iam(NA, NA, M = M, both = TRUE)
        temp200809 <- dmus.dthetas * dthetas.detas
        if (M > 1)
          temp200809[, 2:M] <- temp200809[, 2:M] + sqrt(rpar)
        c(w) * (y - mu) * temp200809
      } else {
        c(w) * (y - mu) *
          cbind(dmus.dthetas[, 1] * dthetas.detas[, 1],
                dmus.dthetas[, 2] * dthetas.detas[, 2] + sqrt(rpar))
      }
    } else {
      temp20101111 <- dmus.dthetas * dthetas.detas
      c(w) * (y - mu) * temp20101111
    }

    myderiv
  }), list( .link1 = link1, .earg1 = earg1,
            .link2 = link2, .earg2 = earg2,
            .firstDeriv = firstDeriv,
            .rpar = rpar, .divisor = divisor ))),

  weight = eval(substitute(expression({
    if ( .oim ) {
      wz <- matrix(0, n, dimm(M))
      wz[, iam(1, 1, M)] <- Xm2
      wz[, iam(1, 2, M)] <- y - 2 * mu
      wz[, iam(2, 2, M)] <- theta1 * (3 * mu - 2 * y) / (theta2 + Xm2)
      wz <- wz * Xm2 / (theta2 + Xm2)^2
    }


    if ( .firstDeriv == "rpar") {
      if (FALSE) {
        wz <-   dmus.dthetas[, index$row] *  dmus.dthetas[, index$col] *
               dthetas.detas[, index$row] * dthetas.detas[, index$col]
        if (M > 1)
          wz[, 2:M] <- wz[, 2:M] + rpar
      } else {
        wz <- cbind(( dmus.dthetas[, 1] * dthetas.detas[, 1])^2,
                    ( dmus.dthetas[, 2] * dthetas.detas[, 2])^2 + rpar,
                      dmus.dthetas[, 1] *  dmus.dthetas[, 2] * 
                     dthetas.detas[, 1] * dthetas.detas[, 2])
      }
    } else {
      run.varcov <- 0
      index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

      mysigma <- 1

      for (ii in 1:( .nsimEIM )) {
        ysim <- theta1 * Xm2 / (theta2 + Xm2) + rnorm(n, sd = mysigma)
        temp3 <- (ysim - mu) * dmus.dthetas * dthetas.detas
        run.varcov <- run.varcov +
                      temp3[, index0$row.index] *
                      temp3[, index0$col.index]
      }
      run.varcov <- run.varcov / .nsimEIM

      wz <- if (intercept.only)
              matrix(colMeans(run.varcov),
                     n, ncol(run.varcov), byrow = TRUE) else run.varcov

    }

    c(w) * wz
  }), list( .link1 = link1, .link2 = link2,
            .firstDeriv = firstDeriv,
            .nsimEIM = nsimEIM, .oim = oim ))))
}








skira.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 skira <- function(rpar = 0.1, divisor = 10,
           init1 = NULL, init2 = NULL,
           link1 = "identitylink", link2 = "identitylink",
           earg1 = list(),
           earg2 = list(),
           imethod = 1,
           oim = TRUE,
           probs.x = c(0.15, 0.85),
           smallno = 1.0e-3,
           nsimEIM = 500,
           firstDeriv = c("nsimEIM", "rpar"),
           dispersion = 0, zero = NULL) {

  firstDeriv <- match.arg(firstDeriv, c("nsimEIM", "rpar"))[1]

  if (!is.Numeric(probs.x, length.arg = 2))
    stop("argument 'probs.x' must be numeric and of length two")

  estimated.dispersion <- dispersion == 0
  if (mode(link1) != "character" && mode(link1) != "name")
    link1 <- as.character(substitute(link1))
  if (mode(link2) != "character" && mode(link2) != "name")
    link2 <- as.character(substitute(link2))

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE))
    stop("argument 'imethod' must be integer")

  if (imethod > 5)
    stop("argument 'imethod' must be 1, 2, 3, 4 or 5")

  if (!is.list(earg1))
    earg1 = list()
  if (!is.list(earg2))
    earg2 = list()

    stopifnot(nsimEIM > 10, length(nsimEIM) == 1,
              nsimEIM == round(nsimEIM))


  new("vglmff",
  blurb = c("Shinozaki-Kira regression model\n",
            "Y_i = 1 / (theta1 + theta2 * u_i) + e_i\n\n",
            "Links:  ",
            namesof("theta1", link1, earg = earg1), ", ",
            namesof("theta2", link2, earg = earg2)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    M <- if (is.matrix(y))
      ncol(y) else 1
    if (residuals) {
      if (M > 1) NULL else (y - mu) * sqrt(w)
    } else {
      ResSS.vgam(y - mu, w, M = M)
    }
  },

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("theta1", "theta2"),
         link1    = .link1 ,
         link2    = .link2 ,
         zero = .zero )
  }, list( .zero = zero,
           .link1 = link1, .link2 = link2
         ))),

  initialize = eval(substitute(expression({

 warning("20101105; need to fix a bug in the signs of initial vals")


    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    if (!length(Xm2)) stop("regressor not found")
    if (ncol(as.matrix(Xm2)) != 1)
      stop("regressor not found or is not a vector. ",
         "Use the 'form2' argument without an intercept")
    Xm2 <- as.vector(Xm2)
    extra$Xm2 <- Xm2

    predictors.names <-
       c(namesof("theta1", .link1 , earg = .earg1, tag = FALSE),
         namesof("theta2", .link2 , earg = .earg2, tag = FALSE))

    if (length(mustart) || length(coefstart))
      stop("cannot handle 'mustart' or 'coefstart'")

    if (!length(etastart)) {


        min.q <- quantile(Xm2, probs = .probs.x[1] )
        max.q <- quantile(Xm2, probs = .probs.x[2] )
      if ( .imethod == 3 || .imethod == 2 ) {

        mysubset <- subset.lohi(Xm2, y, probs.x = .probs.x,
                  type = ifelse( .imethod == 2, "median", "wtmean"),
                  wtvec = w)

        mat.x <- with(mysubset, cbind(c(1, 1),
                      c(x1bar, x2bar)) * c(y1bar, y2bar))
        theta.temp <- solve(mat.x, c(1, 1))
        init1 <- theta.temp[1]
        init2 <- theta.temp[2]
    } else if ( .imethod == 1 ) {
        yy <- as.vector(  y[(Xm2 > min.q) & (Xm2 < max.q)])
        xx <- as.vector(Xm2[(Xm2 > min.q) & (Xm2 < max.q)])
        ww <- as.vector(  w[(Xm2 > min.q) & (Xm2 < max.q)])
        yy[ abs(yy) < .smallno ] <- .smallno *
   sign(yy[ abs(yy) < .smallno ])

        wt.temp <- (yy^4) * ww
        wt.temp.max <- median(wt.temp) * 100
        wt.temp[wt.temp > wt.temp.max] <- wt.temp.max

        mylm.wfit <- lm.wfit(x = cbind(1, xx),
                             y = c(1 / yy),
                             w = c(wt.temp))
        init1 <- mylm.wfit$coef[1]
        init2 <- mylm.wfit$coef[2]
    } else if (( .imethod == 4) || ( .imethod == 5)) {

      tempfit <- if ( .imethod == 4 ) {
        fitted(loess(y ~ Xm2))
      } else {
        fitted(smooth.spline(Xm2, y, w = w, df = 2.0))
      }

      mysubset <- subset.lohi(Xm2, y, probs.x = .probs.x,
                              type = "wtmean", wtvec = w)


      mat.x <- with(mysubset, cbind(c(1, 1),
                    c(x1bar, x2bar)) * c(y1bar, y2bar))
      theta.temp <- solve(mat.x, c(1, 1))
      init1 <- theta.temp[1]
      init2 <- theta.temp[2]
    } else {
      stop("argument 'imethod' unmatched")
    }

    mu <- 1 / (init1 + init2 * Xm2)



 matplot(Xm2, cbind(y, mu), col = c("blue", "green"),
         main = "Initial values in green")
 if ( .imethod == 1 ) {
 points(Xm2, 1 / (init1 + init2 * Xm2), col = "green")
 } else {
 with(mysubset,
 points(c(x1bar, x2bar), c(y1bar, y2bar), col = "red", pch = "+", cex = 2))
 }


      if (length( .init1 )) init1 <- .init1
      if (length( .init2 )) init2 <- .init2
      etastart <- cbind(
          rep(theta2eta(init1, .link1 , earg = .earg1 ), len = n),
          rep(theta2eta(init2, .link2 , earg = .earg2 ), len = n))
    } else {
      stop("cannot handle 'etastart' or 'mustart'")
    }
  }), list( .init1 = init1, .link1 = link1, .earg1 = earg1,
            .init2 = init2, .link2 = link2, .earg2 = earg2,
            .smallno = smallno, .probs.x = probs.x,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    theta1 <- eta2theta(eta[, 1], .link1 , earg = .earg1 )
    theta2 <- eta2theta(eta[, 2], .link2 , earg = .earg2 )
    1 / (theta1 + theta2 * extra$Xm2)
  }, list( .link1 = link1, .earg1 = earg1,
           .link2 = link2, .earg2 = earg2 ))),
  last = eval(substitute(expression({
    misc$link <-    c(theta1 = .link1 , theta2 = .link2 )

    misc$earg <- list(theta1 = .earg1 , theta2 = .earg2 )

    misc$rpar <- rpar
    misc$orig.rpar <- .rpar
    fit$df.residual <- n - rank
    fit$df.total <- n
    dpar <- .dispersion
    if (!dpar) {
      dpar <- sum(c(w) * (y - mu)^2) / (n - ncol.X.vlm)
    }
    misc$dispersion <- dpar
    misc$default.dispersion <- 0
    misc$estimated.dispersion <- .estimated.dispersion

    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
    misc$firstDeriv <- .firstDeriv
    misc$oim <- .oim
    misc$multipleResponses <- FALSE
  }), list( .link1 = link1, .earg1 = earg1,
            .link2 = link2, .earg2 = earg2,
            .dispersion = dispersion, .rpar = rpar,
            .imethod = imethod, .nsimEIM = nsimEIM,
            .firstDeriv = firstDeriv, .oim = oim,
            .estimated.dispersion = estimated.dispersion ))),
  summary.dispersion = FALSE,
  vfamily = c("skira", "vnonlinear"),
  deriv = eval(substitute(expression({
    rpar <- if ( .firstDeriv == "rpar") {
      if (iter > 1) {
        max(rpar / .divisor, 1000 * .Machine$double.eps)
      } else {
        d3 <- deriv3( ~ 1 / (theta1 + theta2 * Xm2),
                     c("theta1", "theta2"), hessian = FALSE)
        .rpar
      }
    } else {
        .rpar
    }

    theta1 <- eta2theta(eta[, 1], .link1 , earg = .earg1 )
    theta2 <- eta2theta(eta[, 2], .link2 , earg = .earg2 )
    dthetas.detas <- cbind(dtheta.deta(theta1, .link1 , earg = .earg1 ),
                           dtheta.deta(theta2, .link2 , earg = .earg2 ))

    dmus.dthetas <- if (FALSE) {
      attr(eval(d3), "gradient")
    } else {
      dmu.dtheta1  <-   -1 / (theta1 + theta2 * Xm2)^2
      dmu.dtheta2  <- -Xm2 / (theta1 + theta2 * Xm2)^2
      cbind(dmu.dtheta1, dmu.dtheta2)
    }


      myderiv <- if ( .firstDeriv == "nsimEIM") {
        c(w) * (y - mu) * dmus.dthetas * dthetas.detas
      } else {
        c(w) * (y - mu) *
        cbind(dmus.dthetas[, 1] * dthetas.detas[, 1],
              dmus.dthetas[, 2] * dthetas.detas[, 2] + sqrt(rpar))
      }
      myderiv
  }), list( .link1 = link1, .earg1 = earg1,
            .link2 = link2, .earg2 = earg2,
            .firstDeriv = firstDeriv,
            .rpar = rpar, .divisor = divisor ))),
  weight = eval(substitute(expression({
    if ( .firstDeriv == "rpar") {
      if (FALSE) {
        index5 <- iam(NA, NA, M = M, both = TRUE)
        wz <- dmus.dthetas[, index5$row] *
              dmus.dthetas[, index5$col] *
             dthetas.detas[, index5$row] *
             dthetas.detas[, index5$col]

        if (M > 1) wz[, -(1:M)] <- wz[, -(1:M)] / 100
      } else {
        wz <- cbind((dmus.dthetas[, 1] * dthetas.detas[, 1])^2,
                    (dmus.dthetas[, 2] * dthetas.detas[, 2])^2 + rpar,
                     dmus.dthetas[, 1] *  dmus.dthetas[, 2] *
                    dthetas.detas[, 1] * dthetas.detas[, 2])
      }
    } else {
      run.varcov <- 0
      index0 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

      mysigma <- sqrt( median( (y - mu)^2 ) ) / 100
      mysigma <- 1

      for (ii in 1:( .nsimEIM )) {
        ysim <- 1 / (theta1 + theta2 * Xm2) + rnorm(n, sd = mysigma)
        temp3 <- (ysim - mu) * dmus.dthetas * dthetas.detas
        run.varcov <- run.varcov +
                      temp3[, index0$row.index] * temp3[, index0$col.index]
      }
      run.varcov <- run.varcov / .nsimEIM

      wz <- if (intercept.only)
              matrix(colMeans(run.varcov),
                     n, ncol(run.varcov), byrow = TRUE) else run.varcov
    }


      c(w) * wz
  }), list( .link1 = link1, .link2 = link2,
            .firstDeriv = firstDeriv,
            .nsimEIM = nsimEIM, .oim = oim ))))
}








