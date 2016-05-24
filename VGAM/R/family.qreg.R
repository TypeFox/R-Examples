# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.
























dlms.bcn <- function(x, lambda = 1, mu = 0, sigma = 1,
                     tol0 = 0.001, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  zedd   <- ((x/mu)^lambda - 1) / (lambda * sigma)
  log.dz.dy <- (lambda - 1) * log(x/mu) - log(mu * sigma)

  is.eff.0 <- abs(lambda) < tol0
  if (any(is.eff.0)) {
    zedd[is.eff.0] <- log(x[is.eff.0] / mu[is.eff.0]) / sigma[is.eff.0]
    log.dz.dy[is.eff.0] <- -log(x[is.eff.0] * sigma[is.eff.0])
  }
  logden <- dnorm(zedd, log = TRUE) + log.dz.dy
  if (log.arg) logden else exp(logden)
}



qlms.bcn <- function(p, lambda = 1, mu = 0, sigma = 1) {

  answer <- mu * (1 + lambda * sigma * qnorm(p))^(1/lambda)
  answer
}







lms.bcn.control <-
lms.bcg.control <-
lms.yjn.control <- function(trace = TRUE, ...)
   list(trace = trace) 



 lms.bcn <- function(percentiles = c(25, 50, 75),
                     zero = c("lambda", "sigma"),
                     llambda = "identitylink",
                     lmu = "identitylink",
                     lsigma = "loge",
                     idf.mu = 4,
                     idf.sigma = 2,
                     ilambda = 1,
                     isigma = NULL,
                     tol0 = 0.001) {
  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")


  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")


  if (!is.Numeric(tol0, positive = TRUE, length.arg = 1))
    stop("bad input for argument 'tol0'")

  if (!is.Numeric(ilambda))
    stop("bad input for argument 'ilambda'")
  if (length(isigma) &&
      !is.Numeric(isigma, positive = TRUE))
    stop("bad input for argument 'isigma'")



  new("vglmff",
  blurb = c("LMS ",
            "quantile",
            " regression (Box-Cox transformation to normality)\n",
            "Links:    ",
            namesof("lambda", link = llambda, earg = elambda), ", ",
            namesof("mu",     link = lmu,     earg = emu), ", ",
            namesof("sigma",  link = lsigma,  earg = esigma)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("lambda", "mu", "sigma"),
         llambda = .llambda ,
         lmu     = .lmu ,
         lsigma  = .lsigma ,
         zero = .zero )
  }, list( .zero = zero,
           .llambda = llambda, .lmu = lmu, .lsigma = lsigma ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
      c(namesof("lambda", .llambda, earg = .elambda, short= TRUE),
        namesof("mu",     .lmu,     earg = .emu,     short= TRUE),
        namesof("sigma",  .lsigma,  earg = .esigma,  short= TRUE))

    if (!length(etastart)) {

        Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                               y = y, w = w, df = .idf.mu )
        fv.init <- c(predict(Fit5, x = x[, min(ncol(x), 2)])$y)

        lambda.init <- if (is.Numeric( .ilambda )) .ilambda else 1.0
        sigma.init <- if (is.null(.isigma)) {
          myratio <- ((y/fv.init)^lambda.init - 1) / lambda.init
          if (is.Numeric( .idf.sigma )) {
            fit600 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                       y = myratio^2,
                                       w = w, df = .idf.sigma)
            sqrt(c(abs(predict(fit600, x = x[, min(ncol(x), 2)])$y)))
          } else {
            sqrt(var(myratio))
          }
        } else {
          .isigma
        }

        etastart <-
          cbind(theta2eta(lambda.init, .llambda, earg = .elambda),
                theta2eta(fv.init,     .lmu,     earg = .emu),
                theta2eta(sigma.init,  .lsigma,  earg = .esigma))
    }
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma, 
            .idf.mu = idf.mu,
            .idf.sigma = idf.sigma,
            .ilambda = ilambda, .isigma = isigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
      eta[, 1] <- eta2theta(eta[, 1], .llambda, earg = .elambda)
      eta[, 2] <- eta2theta(eta[, 2], .lmu,     earg = .emu)
      eta[, 3] <- eta2theta(eta[, 3], .lsigma,  earg = .esigma)
      qtplot.lms.bcn(percentiles = .percentiles, eta = eta)
  }, list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
           .elambda = elambda, .emu = emu, .esigma = esigma, 
           .percentiles = percentiles ))),
  last = eval(substitute(expression({
    misc$links <-    c(lambda = .llambda, mu = .lmu, sigma = .lsigma )

    misc$earg  <- list(lambda = .elambda, mu = .emu, sigma = .esigma )

    misc$tol0 <- .tol0
    misc$percentiles  <- .percentiles  # These are argument values
    misc$true.mu <- FALSE  # @fitted is not a true mu
    if (control$cdf) {
      post$cdf <-
        cdf.lms.bcn(y,
                    eta0 = matrix(c(lambda, mymu, sigma), ncol = 3,
                                  dimnames = list(dimnames(x)[[1]], NULL)))
    }
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma, 
            .percentiles = percentiles,
            .tol0 = tol0 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    muvec  <- eta2theta(eta[, 2], .lmu     , earg = .emu )
    sigma  <- eta2theta(eta[, 3], .lsigma  , earg = .esigma )


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- dlms.bcn(x = y, lambda = lambda, mu = mu, sigma = sigma,
                          tol0 = .tol0 , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
           .elambda = elambda, .emu = emu, .esigma = esigma,
           .tol0 = tol0 ))),
  vfamily = c("lms.bcn", "lmscreg"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    mymu   <- eta2theta(eta[, 2], .lmu, earg = .emu)
    sigma  <- eta2theta(eta[, 3], .lsigma, earg = .esigma)

    zedd <- ((y / mymu)^lambda - 1) / (lambda * sigma)
    z2m1 <- zedd * zedd - 1

    dl.dlambda <- zedd * (zedd - log(y/mymu) / sigma) / lambda -
                  z2m1 * log(y/mymu)
    dl.dmu <- zedd / (mymu * sigma) + z2m1 * lambda / mymu
    dl.dsigma <- z2m1 / sigma

    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    dmu.deta     <- dtheta.deta(mymu,   .lmu     , earg = .emu )
    dsigma.deta  <- dtheta.deta(sigma,  .lsigma  , earg = .esigma )

    c(w) * cbind(dl.dlambda  * dlambda.deta,
                 dl.dmu      * dmu.deta,
                 dl.dsigma   * dsigma.deta)
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, 6)
    wz[,iam(1, 1, M)] <- (7 * sigma^2 / 4) * dlambda.deta^2
    wz[,iam(2, 2, M)] <- (1 + 2*(lambda*sigma)^2)/(mymu*sigma)^2 *
                         dmu.deta^2
    wz[,iam(3, 3, M)] <- (2 / sigma^2) * dsigma.deta^2
    wz[,iam(1, 2, M)] <- (-1 / (2 * mymu)) * dlambda.deta * dmu.deta
    wz[,iam(1, 3, M)] <- (lambda * sigma) * dlambda.deta * dsigma.deta
    wz[,iam(2, 3, M)] <- (2*lambda/(mymu * sigma)) *
                         dmu.deta * dsigma.deta
    c(w) * wz
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma ))))
}  # End of lms.bcn






 lms.bcg <- function(percentiles = c(25, 50, 75),
                     zero = c("lambda", "sigma"),
                     llambda = "identitylink",
                     lmu = "identitylink",
                     lsigma = "loge",
                     idf.mu = 4,
                     idf.sigma = 2,
                     ilambda = 1,
                     isigma = NULL) {
  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")


    if (!is.Numeric(ilambda))
      stop("bad input for argument 'ilambda'")
    if (length(isigma) && !is.Numeric(isigma, positive = TRUE))
      stop("bad input for argument 'isigma'")

  new("vglmff",
  blurb = c("LMS Quantile Regression ",
            "(Box-Cox transformation to a Gamma distribution)\n",
            "Links:    ",
            namesof("lambda", link = llambda, earg = elambda), ", ",
            namesof("mu",     link = lmu,     earg = emu), ", ",
            namesof("sigma",  link = lsigma,  earg = esigma)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list(.zero = zero))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("lambda", "mu", "sigma"),
         llambda = .llambda ,
         lmu     = .lmu ,
         lsigma  = .lsigma ,
         zero = .zero )
  }, list( .zero = zero,
           .llambda = llambda, .lmu = lmu, .lsigma = lsigma ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1, ncol.y.max = 1)

        predictors.names <- c(
            namesof("lambda", .llambda, earg = .elambda, short = TRUE),
            namesof("mu",     .lmu,     earg = .emu,     short = TRUE),
            namesof("sigma",  .lsigma,  earg = .esigma,  short = TRUE))

        if (!length(etastart)) {

          Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                 y = y, w = w, df = .idf.mu )
          fv.init <- c(predict(Fit5, x = x[, min(ncol(x), 2)])$y)

          lambda.init <- if (is.Numeric( .ilambda )) .ilambda else 1.0

          sigma.init <- if (is.null( .isigma )) {
            myratio <- ((y/fv.init)^lambda.init-1) / lambda.init
            if (is.numeric( .idf.sigma ) &&
                is.finite( .idf.sigma )) {
              fit600 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                       y = (myratio)^2,
                                       w = w, df = .idf.sigma )
              sqrt(c(abs(predict(fit600, x = x[, min(ncol(x), 2)])$y)))
            } else {
              sqrt(var(myratio))
            }
            } else .isigma

            etastart <-
              cbind(theta2eta(lambda.init,  .llambda , earg = .elambda ),
                    theta2eta(fv.init,      .lmu ,     earg = .emu ),
                    theta2eta(sigma.init,   .lsigma ,  earg = .esigma ))
        }
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma, 
            .idf.mu = idf.mu,
            .idf.sigma = idf.sigma,
            .ilambda = ilambda, .isigma = isigma ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta[, 1] <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    eta[, 2] <- eta2theta(eta[, 2], .lmu ,     earg = .emu )
    eta[, 3] <- eta2theta(eta[, 3], .lsigma ,  earg = .esigma )
    qtplot.lms.bcg(percentiles = .percentiles, eta = eta)
  }, list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
           .elambda = elambda, .emu = emu, .esigma = esigma, 
           .percentiles = percentiles ))),
  last = eval(substitute(expression({
    misc$link <-    c(lambda = .llambda, mu = .lmu, sigma = .lsigma )

    misc$earg <- list(lambda = .elambda, mu = .emu, sigma = .esigma )

     misc$percentiles <- .percentiles  # These are argument values
    misc$true.mu <- FALSE    # $fitted is not a true mu
    if (control$cdf) {
      post$cdf <- cdf.lms.bcg(y, eta0 = matrix(c(lambda, mymu, sigma),
          ncol = 3, dimnames = list(dimnames(x)[[1]], NULL)))
    }
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma, 
            .percentiles = percentiles ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    mu     <- eta2theta(eta[, 2], .lmu, earg = .emu)
    sigma  <- eta2theta(eta[, 3], .lsigma, earg = .esigma)
    Gee <- (y / mu)^lambda
    theta <- 1 / (sigma * lambda)^2
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (log(abs(lambda)) + theta * (log(theta) +
                         log(Gee)-Gee) - lgamma(theta) - log(y))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
           .elambda = elambda, .emu = emu, .esigma = esigma ))),
  vfamily = c("lms.bcg", "lmscreg"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    mymu   <- eta2theta(eta[, 2], .lmu,     earg = .emu)
    sigma  <- eta2theta(eta[, 3], .lsigma,  earg = .esigma)

    Gee <- (y / mymu)^lambda
    theta <- 1 / (sigma * lambda)^2
    dd <- digamma(theta)

    dl.dlambda <- (1 + 2 * theta * (dd + Gee -1 -log(theta) -
                 0.5 * (Gee + 1) * log(Gee))) / lambda
    dl.dmu <- lambda * theta * (Gee-1) / mymu
    dl.dsigma <- 2*theta*(dd + Gee - log(theta * Gee)-1) / sigma

    dlambda.deta <- dtheta.deta(lambda, link = .llambda, earg = .elambda)
    dmu.deta     <- dtheta.deta(mymu, link = .lmu, earg = .emu)
    dsigma.deta  <- dtheta.deta(sigma, link = .lsigma, earg = .esigma)

    cbind(dl.dlambda * dlambda.deta,
          dl.dmu     * dmu.deta,
          dl.dsigma  * dsigma.deta) * w
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma ))),
  weight = eval(substitute(expression({
    tritheta <- trigamma(theta)
    wz <- matrix(0, n, 6)

    if (TRUE) {
        part2 <- dd + 2/theta - 2*log(theta)
        wz[,iam(1, 1, M)] <- ((1 + theta*(tritheta*(1+4*theta) -
                           4*(1+1/theta) - log(theta)*(2/theta -
                           log(theta)) + dd*part2)) / lambda^2) *
                           dlambda.deta^2
    } else {
        temp <- mean( Gee*(log(Gee))^2 )
        wz[,iam(1, 1, M)] <- ((4 * theta * (theta * tritheta-1) - 1 +
                          theta*temp) / lambda^2) * dlambda.deta^2
    }

    wz[,iam(2, 2, M)] <- dmu.deta^2 / (mymu * sigma)^2
    wz[,iam(3, 3, M)] <- (4 * theta * (theta * tritheta - 1) / sigma^2) *
                      dsigma.deta^2
    wz[,iam(1, 2, M)] <- (-theta * (dd + 1 / theta - log(theta)) / mymu) *
                      dlambda.deta * dmu.deta
    wz[,iam(1, 3, M)] <- 2 * theta^1.5 * (2 * theta * tritheta - 2 -
                      1 / theta) * dlambda.deta * dsigma.deta
    c(w) * wz
  }), list( .llambda = llambda, .lmu = lmu, .lsigma = lsigma,
            .elambda = elambda, .emu = emu, .esigma = esigma ))))
}




dy.dpsi.yeojohnson <- function(psi, lambda) {

    L <- max(length(psi), length(lambda))
    if (length(psi)    != L) psi    <- rep(psi,    length.out = L)
    if (length(lambda) != L) lambda <- rep(lambda, length.out = L)

    ifelse(psi > 0, (1 + psi * lambda)^(1/lambda - 1),
                    (1 - (2-lambda) * psi)^((lambda - 1) / (2-lambda)))
}


dyj.dy.yeojohnson <- function(y, lambda) {
    L <- max(length(y), length(lambda))
    if (length(y)      != L) y      <- rep(y,      length.out = L)
    if (length(lambda) != L) lambda <- rep(lambda, length.out = L)

    ifelse(y>0, (1 + y)^(lambda - 1), (1 - y)^(1 - lambda))
}


 yeo.johnson <- function(y, lambda, derivative = 0,
                        epsilon = sqrt(.Machine$double.eps),
                        inverse = FALSE) {

    if (!is.Numeric(derivative, length.arg = 1,
                    integer.valued = TRUE) ||
        derivative < 0)
      stop("argument 'derivative' must be a non-negative integer")

    ans <- y
    if (!is.Numeric(epsilon, length.arg = 1, positive = TRUE))
      stop("argument 'epsilon' must be a single positive number")
    L <- max(length(lambda), length(y))
    if (length(y) != L)
      y <- rep(y, length.out = L)
    if (length(lambda) != L)
      lambda <- rep(lambda, length.out = L)  # lambda may be of length 1

    if (inverse) {
        if (derivative != 0)
          stop("argument 'derivative' must 0 when inverse = TRUE")
        if (any(index <- y >= 0 & abs(lambda) > epsilon))
          ans[index] <- (y[index]*lambda[index] + 1)^(1/lambda[index]) - 1
        if (any(index <- y >= 0 & abs(lambda) <= epsilon))
          ans[index] <- expm1(y[index])
        if (any(index <- y <  0 & abs(lambda-2) > epsilon))
          ans[index] <- 1- (-(2-lambda[index]) *
                           y[index]+1)^(1/(2-lambda[index]))
        if (any(index <- y <  0 & abs(lambda-2) <= epsilon))
            ans[index] <- -expm1(-y[index])
        return(ans)
    }
    if (derivative == 0) {
        if (any(index <- y >= 0 & abs(lambda) > epsilon))
          ans[index] <- ((y[index]+1)^(lambda[index]) - 1) / lambda[index]
        if (any(index <- y >= 0 & abs(lambda) <= epsilon))
          ans[index] <- log1p(y[index])
        if (any(index <- y <  0 & abs(lambda-2) > epsilon))
          ans[index] <- -((-y[index]+1)^(2-lambda[index]) - 1)/(2 -
                         lambda[index])
        if (any(index <- y <  0 & abs(lambda-2) <= epsilon))
          ans[index] <- -log1p(-y[index])
    } else {
        psi <- Recall(y = y, lambda=lambda, derivative=derivative-1,
                      epsilon=epsilon, inverse=inverse)
        if (any(index <- y >= 0 & abs(lambda) > epsilon))
          ans[index] <- ( (y[index]+1)^(lambda[index]) *
                        (log1p(y[index]))^(derivative) - derivative *
                        psi[index] ) / lambda[index]
        if (any(index <- y >= 0 & abs(lambda) <= epsilon))
          ans[index] <- (log1p(y[index]))^(derivative + 1) / (derivative+1)
        if (any(index <- y <  0 & abs(lambda-2) > epsilon))
          ans[index] <- -( (-y[index]+1)^(2-lambda[index]) *
                        (-log1p(-y[index]))^(derivative) - derivative *
                        psi[index] ) / (2-lambda[index])
        if (any(index <- y <  0 & abs(lambda-2) <= epsilon))
          ans[index] <- (-log1p(-y[index]))^(derivative + 1) / (derivative+1)
    }
    ans
}


dpsi.dlambda.yjn <- function(psi, lambda, mymu, sigma,
                            derivative = 0, smallno = 1.0e-8) {

    if (!is.Numeric(derivative, length.arg = 1,
                    integer.valued = TRUE) ||
        derivative < 0)
      stop("argument 'derivative' must be a non-negative integer")
    if (!is.Numeric(smallno, length.arg = 1, positive = TRUE))
      stop("argument 'smallno' must be a single positive number")

    L <- max(length(psi), length(lambda), length(mymu), length(sigma))
    if (length(psi)    != L) psi    <- rep(psi,    length.out = L)
    if (length(lambda) != L) lambda <- rep(lambda, length.out = L)
    if (length(mymu)   != L) mymu   <- rep(mymu,   length.out = L)
    if (length(sigma)  != L) sigma  <- rep(sigma,  length.out = L)

    answer <- matrix(NA_real_, L, derivative+1)
    CC <- psi >= 0
    BB <- ifelse(CC, lambda, -2+lambda)
    AA <- psi * BB 
    temp8 <- if (derivative > 0) {
        answer[,1:derivative] <-
            Recall(psi = psi, lambda = lambda, mymu = mymu, sigma = sigma,
                   derivative = derivative-1, smallno = smallno) 
        answer[,derivative] * derivative
    } else { 
        0
    }
    answer[,1+derivative] <- ((AA+1) * (log1p(AA)/BB)^derivative - temp8) / BB

    pos <- (CC & abs(lambda) <= smallno) | (!CC & abs(lambda-2) <= smallno)
    if (any(pos)) 
      answer[pos,1+derivative] =
        (answer[pos, 1]^(1+derivative))/(derivative+1)
    answer
}


gh.weight.yjn.11 <- function(z, lambda, mymu, sigma, derivmat = NULL) {


    if (length(derivmat)) {
      ((derivmat[, 2]/sigma)^2 +
        sqrt(2) * z * derivmat[, 3] / sigma) / sqrt(pi)
    } else {
        # Long-winded way 
        psi <- mymu + sqrt(2) * sigma * z
        (1 / sqrt(pi)) *
        (dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                          derivative = 1)[, 2]^2 +
        (psi - mymu) * 
        dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                         derivative = 2)[, 3]) / sigma^2
    }
}


gh.weight.yjn.12 <- function(z, lambda, mymu, sigma, derivmat = NULL) {
    if (length(derivmat)) {
        (-derivmat[, 2]) / (sqrt(pi) * sigma^2)
    } else {
        psi <- mymu + sqrt(2) * sigma * z
        (1 / sqrt(pi)) * (- dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                                             derivative = 1)[, 2]) / sigma^2
    }
}


gh.weight.yjn.13 <- function(z, lambda, mymu, sigma, derivmat = NULL) {
    if (length(derivmat)) {
        sqrt(8 / pi) * (-derivmat[, 2]) * z / sigma^2
    } else {
        psi <- mymu + sqrt(2) * sigma * z
        (1 / sqrt(pi)) *
        (-2 * dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                               derivative = 1)[, 2]) *
        (psi - mymu) / sigma^3
    }
}


glag.weight.yjn.11 <- function(z, lambda, mymu, sigma, derivmat = NULL) {


  if (length(derivmat)) {
    derivmat[, 4] * (derivmat[, 2]^2 + sqrt(2) * sigma * z * derivmat[, 3])
  } else {
    psi <- mymu + sqrt(2) * sigma * z
    discontinuity <- -mymu / (sqrt(2) * sigma)
    (1 / (2 * sqrt((z-discontinuity^2)^2 + discontinuity^2))) *
    (1 / sqrt(pi)) *
    (dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative = 1)[, 2]^2 +
    (psi - mymu) * 
    dpsi.dlambda.yjn(psi, lambda, mymu,
                     sigma, derivative = 2)[, 3]) / sigma^2
  }
}

glag.weight.yjn.12 <- function(z, lambda, mymu, sigma, derivmat = NULL) {
  discontinuity <- -mymu / (sqrt(2) * sigma)
  if (length(derivmat)) {
    derivmat[, 4] * (-derivmat[, 2])
  } else {
    psi <- mymu + sqrt(2) * sigma * z
    (1 / (2 * sqrt((z-discontinuity^2)^2 + discontinuity^2))) *
    (1 / sqrt(pi)) *
    (- dpsi.dlambda.yjn(psi, lambda, mymu,
                        sigma, derivative = 1)[, 2]) / sigma^2
  }
}

glag.weight.yjn.13 <- function(z, lambda, mymu, sigma, derivmat = NULL) {
  if (length(derivmat)) {
    derivmat[, 4] * (-derivmat[, 2]) * sqrt(8) * z
  } else {
    psi <- mymu + sqrt(2) * sigma * z
    discontinuity <- -mymu / (sqrt(2) * sigma)
    (1 / (2 * sqrt((z-discontinuity^2)^2 + discontinuity^2))) *
    (1 / sqrt(pi)) *
    (-2 * dpsi.dlambda.yjn(psi, lambda, mymu,
                           sigma, derivative = 1)[, 2]) *
    (psi - mymu) / sigma^3
  }
}


gleg.weight.yjn.11 <- function(z, lambda, mymu, sigma, derivmat = NULL) {




  if (length(derivmat)) {
    derivmat[, 4] * (derivmat[, 2]^2 + sqrt(2) * sigma * z * derivmat[, 3])
  } else {
    psi <- mymu + sqrt(2) * sigma * z
    (exp(-z^2) / sqrt(pi)) *
    (dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative = 1)[, 2]^2 +
    (psi - mymu) * 
    dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                     derivative = 2)[, 3]) / sigma^2
  }
}


gleg.weight.yjn.12 <- function(z, lambda, mymu, sigma, derivmat = NULL) {
  if (length(derivmat)) {
    derivmat[, 4] * (- derivmat[, 2])
  } else {
    psi <- mymu + sqrt(2) * sigma * z
    (exp(-z^2) / sqrt(pi)) *
    (- dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                        derivative = 1)[, 2]) / sigma^2
  }
}


gleg.weight.yjn.13 <- function(z, lambda, mymu, sigma, derivmat = NULL) {
  if (length(derivmat)) {
    derivmat[, 4] * (-derivmat[, 2]) * sqrt(8) * z
  } else {
    psi <- mymu + sqrt(2) * sigma * z
    (exp(-z^2) / sqrt(pi)) *
    (-2 * dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative = 1)[, 2]) *
    (psi - mymu) / sigma^3
  }
}




lms.yjn2.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}

 lms.yjn2 <- function(percentiles = c(25, 50, 75),
                      zero = c("lambda", "sigma"),
                      llambda = "identitylink",
                      lmu = "identitylink",
                      lsigma = "loge",
                      idf.mu = 4,
                      idf.sigma = 2,
                      ilambda = 1.0,
                      isigma = NULL,
                      yoffset = NULL,
                      nsimEIM = 250) {

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")



  if (!is.Numeric(ilambda))
    stop("bad input for argument 'ilambda'")
  if (length(isigma) &&
      !is.Numeric(isigma, positive = TRUE))
    stop("bad input for argument 'isigma'")

  new("vglmff",
  blurb = c("LMS Quantile Regression (Yeo-Johnson transformation",
            " to normality)\n",
            "Links:    ",
            namesof("lambda", link = llambda, earg = elambda), ", ",
            namesof("mu",     link = lmu,     earg = emu    ), ", ",
            namesof("sigma",  link = lsigma,  earg = esigma )),
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
         parameters.names = c("lambda", "mu", "sigma"),
         llambda = .llambda ,
         lmu     = .lmu ,
         lsigma  = .lsigma ,
         zero = .zero )
  }, list( .zero = zero,
           .llambda = llambda, .lmu = lmu, .lsigma = lsigma ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)


    predictors.names <-
      c(namesof("lambda", .llambda, earg = .elambda, short= TRUE),
        namesof("mu",     .lmu,     earg = .emu,     short= TRUE),
        namesof("sigma",  .lsigma, earg = .esigma,  short= TRUE))

      y.save <- y
      yoff <- if (is.Numeric( .yoffset)) .yoffset else -median(y) 
      extra$yoffset <- yoff
      y <- y + yoff

      if (!length(etastart)) {
        lambda.init <- if (is.Numeric( .ilambda )) .ilambda else 1.

        y.tx <- yeo.johnson(y, lambda.init)
        fv.init =
        if (smoothok <-
         (length(unique(sort(x[, min(ncol(x), 2)]))) > 7)) {
          fit700 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                   y = y.tx, w = w, df = .idf.mu )
          c(predict(fit700, x = x[, min(ncol(x), 2)])$y)
        } else {
          rep(weighted.mean(y, w), length.out = n)
        }

        sigma.init <- if (!is.Numeric(.isigma)) {
                     if (is.Numeric( .idf.sigma) && smoothok) {
                     fit710 = vsmooth.spline(x = x[, min(ncol(x), 2)],
                                      y = (y.tx - fv.init)^2,
                                      w = w, df = .idf.sigma)
                          sqrt(c(abs(predict(fit710,
                               x = x[, min(ncol(x), 2)])$y)))
                   } else {
                    sqrt( sum( w * (y.tx - fv.init)^2 ) / sum(w) )
                   }
       } else
           .isigma

      etastart <- matrix(0, n, 3)
      etastart[, 1] <- theta2eta(lambda.init, .llambda, earg = .elambda)
      etastart[, 2] <- theta2eta(fv.init,     .lmu,     earg = .emu)
      etastart[, 3] <- theta2eta(sigma.init,  .lsigma,  earg = .esigma)

      }
  }), list(.llambda = llambda, .lmu = lmu, .lsigma = lsigma,
           .elambda = elambda, .emu = emu, .esigma = esigma, 
           .idf.mu = idf.mu,
           .idf.sigma = idf.sigma,
           .ilambda = ilambda,
           .yoffset=yoffset,
           .isigma = isigma))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta[, 1] <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    eta[, 3] <- eta2theta(eta[, 3], .lsigma, earg = .esigma)
    qtplot.lms.yjn(percentiles = .percentiles, eta = eta,
                   yoffset = extra$yoff)
  }, list(.percentiles = percentiles,
          .esigma = esigma, .elambda = elambda,
          .llambda = llambda,
          .lsigma = lsigma))),
  last = eval(substitute(expression({
    misc$link <-    c(lambda = .llambda, mu = .lmu, sigma = .lsigma)
    misc$earg <- list(lambda = .elambda, mu = .emu, sigma = .esigma)

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$percentiles  <- .percentiles  # These are argument values

    misc$true.mu <- FALSE # $fitted is not a true mu
    misc[["yoffset"]] <- extra$yoffset

    y <- y.save   # Restore back the value; to be attached to object

    if (control$cdf) {
        post$cdf <- cdf.lms.yjn(y + misc$yoffset,
            eta0=matrix(c(lambda,mymu,sigma), 
            ncol=3, dimnames = list(dimnames(x)[[1]], NULL)))
    }
  }), list(.percentiles = percentiles,
           .elambda = elambda, .emu = emu, .esigma = esigma, 
           .nsimEIM=nsimEIM,
           .llambda = llambda, .lmu = lmu, .lsigma = lsigma ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    mu     <- eta2theta(eta[, 2], .lmu     , earg = .emu )
    sigma  <- eta2theta(eta[, 3], .lsigma  , earg = .esigma )
    psi <- yeo.johnson(y, lambda)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (-log(sigma) - 0.5 * ((psi-mu)/sigma)^2 +
                        (lambda-1) * sign(y) * log1p(abs(y)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elambda = elambda, .emu = emu, .esigma = esigma, 
           .llambda = llambda, .lmu = lmu,
           .lsigma = lsigma ))),
  vfamily = c("lms.yjn2", "lmscreg"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    mymu <- eta2theta(eta[, 2], .lmu, earg = .emu)
    sigma <- eta2theta(eta[, 3], .lsigma, earg = .esigma)
    dlambda.deta <- dtheta.deta(lambda, link = .llambda, earg = .elambda)
    dmu.deta <- dtheta.deta(mymu, link = .lmu, earg = .emu)
    dsigma.deta <- dtheta.deta(sigma, link = .lsigma, earg = .esigma)

    psi <- yeo.johnson(y, lambda)
    d1 <- yeo.johnson(y, lambda, deriv = 1)
    AA <- (psi - mymu) / sigma 
    dl.dlambda <- -AA * d1 /sigma + sign(y) * log1p(abs(y))
    dl.dmu <- AA / sigma 
    dl.dsigma <- (AA^2 -1) / sigma
    dthetas.detas <- cbind(dlambda.deta, dmu.deta, dsigma.deta)
    c(w) * cbind(dl.dlambda, dl.dmu, dl.dsigma) * dthetas.detas
  }), list( .elambda = elambda, .emu = emu, .esigma = esigma, 
            .llambda = llambda, .lmu = lmu,
               .lsigma = lsigma ))),
  weight = eval(substitute(expression({


    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for (ii in 1:( .nsimEIM )) {
        psi <- rnorm(n, mymu, sigma)
        ysim <- yeo.johnson(y=psi, lam=lambda, inv = TRUE)
        d1 <- yeo.johnson(ysim, lambda, deriv = 1)
        AA <- (psi - mymu) / sigma 
        dl.dlambda <- -AA * d1 /sigma + sign(ysim) * log1p(abs(ysim))
        dl.dmu <- AA / sigma 
        dl.dsigma <- (AA^2 -1) / sigma
        rm(ysim)
        temp3 <- cbind(dl.dlambda, dl.dmu, dl.dsigma)
        run.varcov <- ((ii-1) * run.varcov +
                   temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
    }

        if (intercept.only)
            run.varcov <- matrix(colMeans(run.varcov),
                                nr=n, nc=ncol(run.varcov), byrow = TRUE)


    wz <- run.varcov * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
    dimnames(wz) <- list(rownames(wz), NULL)  # Remove the colnames
    c(w) * wz
  }), list(.lsigma = lsigma,
           .esigma = esigma, .elambda = elambda,
           .nsimEIM=nsimEIM,
           .llambda = llambda))))
}


 lms.yjn <- function(percentiles = c(25, 50, 75),
                    zero = c("lambda", "sigma"), 
                    llambda = "identitylink",
                    lsigma = "loge",
                    idf.mu = 4,
                    idf.sigma = 2,
                    ilambda = 1.0,
                    isigma = NULL,
                    rule = c(10, 5),
                    yoffset = NULL,
                    diagW = FALSE, iters.diagW = 6) {




  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lsigma <- as.list(substitute(lsigma))
  esigma <- link2list(lsigma)
  lsigma <- attr(esigma, "function.name")



  rule <- rule[1]  # Number of points (common) for all the quadrature schemes
  if (rule != 5 && rule != 10)
    stop("only rule=5 or 10 is supported")

  new("vglmff",
  blurb = c("LMS Quantile Regression ",
            "(Yeo-Johnson transformation to normality)\n",
            "Links:    ",
            namesof("lambda", link = llambda, earg = elambda), ", ",
            namesof("mu",     link = "identitylink", earg = list()), ", ",
            namesof("sigma",  link = lsigma,  earg = esigma)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list(.zero = zero))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("lambda", "mu", "sigma"),
         llambda = .llambda ,
         lmu     = "identitylink",
         lsigma  = .lsigma ,
         zero = .zero )
  }, list( .zero = zero,
           .llambda = llambda, .lsigma = lsigma ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1)



    predictors.names <-
      c(namesof("lambda", .llambda, earg = .elambda , short = TRUE),
                "mu",
        namesof("sigma",  .lsigma, earg = .esigma ,  short = TRUE))

    y.save <- y
    yoff <- if (is.Numeric( .yoffset )) .yoffset else -median(y) 
    extra$yoffset <- yoff
    y <- y + yoff

    if (!length(etastart)) {

          lambda.init <- if (is.Numeric( .ilambda )) .ilambda else 1.0

            y.tx <- yeo.johnson(y, lambda.init)
            if (smoothok <-
               (length(unique(sort(x[, min(ncol(x), 2)]))) > 7)) {
                fit700 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                        y = y.tx, w = w, df = .idf.mu )
                fv.init <- c(predict(fit700, x = x[, min(ncol(x), 2)])$y)
            } else {
                fv.init <- rep(weighted.mean(y, w), length.out = n)
            }

            sigma.init <- if (!is.Numeric( .isigma )) {
                           if (is.Numeric( .idf.sigma) &&
                               smoothok) {
                           fit710 = vsmooth.spline(x = x[, min(ncol(x), 2)],
                                      y = (y.tx - fv.init)^2,
                                      w = w, df = .idf.sigma)
                             sqrt(c(abs(predict(fit710,
                                        x = x[, min(ncol(x), 2)])$y)))
                           } else {
                             sqrt( sum( w * (y.tx - fv.init)^2 ) / sum(w) )
                           }
                         } else
                             .isigma

            etastart <-
              cbind(theta2eta(lambda.init,.llambda, earg = .elambda),
                    fv.init,
                    theta2eta(sigma.init, .lsigma, earg = .esigma))

        }
  }), list(.lsigma = lsigma,
           .llambda = llambda,
           .esigma = esigma, .elambda = elambda,
           .idf.mu = idf.mu,
           .idf.sigma = idf.sigma,
           .ilambda = ilambda,
           .yoffset=yoffset,
           .isigma = isigma))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta[, 1] <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    eta[, 3] <- eta2theta(eta[, 3], .lsigma, earg = .esigma)
    qtplot.lms.yjn(percentiles = .percentiles,
                   eta = eta, yoffset = extra$yoff)
  }, list(.percentiles = percentiles,
          .esigma = esigma,
          .elambda = elambda,
          .llambda = llambda,
          .lsigma = lsigma))),
  last = eval(substitute(expression({
    misc$link <-    c(lambda = .llambda, mu = "identitylink",
                      sigma = .lsigma)

    misc$earg <- list(lambda = .elambda, mu = list(theta = NULL),
                      sigma = .esigma)

    misc$percentiles  <- .percentiles  # These are argument values
    misc$true.mu <- FALSE    # $fitted is not a true mu
    misc[["yoffset"]] <- extra$yoff

    y <- y.save   # Restore back the value; to be attached to object

    if (control$cdf) {
        post$cdf =
          cdf.lms.yjn(y + misc$yoffset,
                      eta0 = matrix(c(lambda,mymu,sigma), 
                      ncol = 3,
                      dimnames = list(dimnames(x)[[1]], NULL)))
    }
  }), list(.percentiles = percentiles,
           .esigma = esigma, .elambda = elambda,
          .llambda = llambda,
          .lsigma = lsigma))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    lambda <- eta2theta(eta[, 1], .llambda , earg = .elambda )
    mu <- eta[, 2]
    sigma <- eta2theta(eta[, 3], .lsigma , earg = .esigma )
    psi <- yeo.johnson(y, lambda)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (-log(sigma) - 0.5 * ((psi-mu)/sigma)^2 +
                        (lambda-1) * sign(y) * log1p(abs(y)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .esigma = esigma, .elambda = elambda,
           .lsigma = lsigma, .llambda = llambda))),
  vfamily = c("lms.yjn", "lmscreg"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta[, 1], .llambda, earg = .elambda)
    mymu <- eta[, 2]
    sigma <- eta2theta(eta[, 3], .lsigma, earg = .esigma)

    psi <- yeo.johnson(y, lambda)
    d1 <- yeo.johnson(y, lambda, deriv = 1)
    AA <- (psi - mymu) / sigma 

    dl.dlambda <- -AA * d1 /sigma + sign(y) * log1p(abs(y))
    dl.dmu <- AA / sigma 
    dl.dsigma <- (AA^2 -1) / sigma
    dlambda.deta <- dtheta.deta(lambda, link = .llambda, earg = .elambda)
    dsigma.deta <- dtheta.deta(sigma, link = .lsigma, earg = .esigma)

    cbind(dl.dlambda * dlambda.deta,
          dl.dmu,
          dl.dsigma * dsigma.deta) * c(w)
  }), list( .esigma = esigma, .elambda = elambda,
            .lsigma = lsigma, .llambda = llambda ))),
  weight = eval(substitute(expression({
    wz <- matrix(0, n, 6)


        wz[,iam(2, 2, M)] <- 1 / sigma^2
        wz[,iam(3, 3, M)] <- 2 * wz[,iam(2, 2, M)]   # 2 / sigma^2


    if (.rule == 10) {
    glag.abs = c(0.13779347054,0.729454549503,
                 1.80834290174,3.40143369785,
                 5.55249614006,8.33015274676,
                 11.8437858379,16.2792578314,
                 21.996585812, 29.9206970123)
    glag.wts = c(0.308441115765, 0.401119929155, 0.218068287612,
                 0.0620874560987, 0.00950151697517, 0.000753008388588, 
                 2.82592334963e-5,
                 4.24931398502e-7, 1.83956482398e-9, 9.91182721958e-13)
    } else {
    glag.abs = c(0.2635603197180449, 1.4134030591060496,
                  3.5964257710396850,
                 7.0858100058570503, 12.6408008442729685)
    glag.wts = c(5.217556105826727e-01, 3.986668110832433e-01,
                 7.594244968176882e-02,
                 3.611758679927785e-03, 2.336997238583738e-05)
    }

    if (.rule == 10) {
    sgh.abs = c(0.03873852801690856, 0.19823332465268367,
                  0.46520116404433082,
                0.81686197962535023, 1.23454146277833154,
                  1.70679833036403172,
                2.22994030591819214, 2.80910399394755972,
                  3.46387269067033854,
                4.25536209637269280)
    sgh.wts = c(9.855210713854302e-02, 2.086780884700499e-01,
                 2.520517066468666e-01,
         1.986843323208932e-01,9.719839905023238e-02,
                 2.702440190640464e-02,
         3.804646170194185e-03, 2.288859354675587e-04,
                  4.345336765471935e-06,
         1.247734096219375e-08)
    } else {
  sgh.abs = c(0.1002421519682381, 0.4828139660462573,
                  1.0609498215257607,
              1.7797294185202606, 2.6697603560875995)
  sgh.wts = c(0.2484061520284881475,0.3923310666523834311,
                 0.2114181930760276606,
            0.0332466603513424663, 0.0008248533445158026)
    }

    if (.rule == 10) {
        gleg.abs = c(-0.973906528517, -0.865063366689, -0.679409568299,
                     -0.433395394129, -0.148874338982)
        gleg.abs = c(gleg.abs, rev(-gleg.abs))
        gleg.wts = c(0.0666713443087, 0.149451349151, 0.219086362516,
                     0.26926671931, 0.295524224715)
        gleg.wts = c(gleg.wts, rev(gleg.wts))
    } else {
        gleg.abs = c(-0.9061798459386643,-0.5384693101056820, 0,
                      0.5384693101056828, 0.9061798459386635)
        gleg.wts = c(0.2369268850561853,0.4786286704993680,
                 0.5688888888888889,
                   0.4786286704993661, 0.2369268850561916)
    }


    discontinuity = -mymu/(sqrt(2)*sigma)


    LL <- pmin(discontinuity, 0)
    UU <- pmax(discontinuity, 0)
    if (FALSE) {
      AA <- (UU-LL)/2
      for (kk in 1:length(gleg.wts)) {
        temp1 <- AA * gleg.wts[kk] 
        abscissae <- (UU+LL)/2 + AA * gleg.abs[kk]
        psi <- mymu + sqrt(2) * sigma * abscissae
        temp9 <- dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                                  derivative = 2)
        temp9 <- cbind(temp9, exp(-abscissae^2) / (sqrt(pi) * sigma^2))

        wz[,iam(1, 1, M)] <- wz[,iam(1, 1, M)] + temp1 *
              gleg.weight.yjn.11(abscissae, lambda, mymu, sigma, temp9)
        wz[,iam(1, 2, M)] <- wz[,iam(1, 2, M)] + temp1 *
              gleg.weight.yjn.12(abscissae, lambda, mymu, sigma, temp9)
        wz[,iam(1, 3, M)] <- wz[,iam(1, 3, M)] + temp1 *
              gleg.weight.yjn.13(abscissae, lambda, mymu, sigma, temp9)
      }
    } else {
      temp9 <- .Fortran("yjngintf", as.double(LL),
                 as.double(UU),
                 as.double(gleg.abs), as.double(gleg.wts), as.integer(n),
                 as.integer(length(gleg.abs)), as.double(lambda),
                 as.double(mymu), as.double(sigma), answer = double(3*n),
                     eps=as.double(1.0e-5))$ans
      dim(temp9) <- c(3,n)
      wz[,iam(1, 1, M)] <- temp9[1,]
      wz[,iam(1, 2, M)] <- temp9[2,]
      wz[,iam(1, 3, M)] <- temp9[3,]
    }



    for (kk in 1:length(sgh.wts)) {

      abscissae <- sign(-discontinuity) * sgh.abs[kk]
      psi <- mymu + sqrt(2) * sigma * abscissae   # abscissae = z
      temp9 <- dpsi.dlambda.yjn(psi, lambda, mymu, sigma,
                                 derivative = 2)
      wz[,iam(1, 1, M)] <- wz[,iam(1, 1, M)] + sgh.wts[kk] * 
            gh.weight.yjn.11(abscissae, lambda, mymu, sigma, temp9)
      wz[,iam(1, 2, M)] <- wz[,iam(1, 2, M)] + sgh.wts[kk] * 
            gh.weight.yjn.12(abscissae, lambda, mymu, sigma, temp9)
      wz[,iam(1, 3, M)] <- wz[,iam(1, 3, M)] + sgh.wts[kk] * 
            gh.weight.yjn.13(abscissae, lambda, mymu, sigma, temp9)
    }

    temp1 <- exp(-discontinuity^2)
    for (kk in 1:length(glag.wts)) {
      abscissae <- sign(discontinuity) * sqrt(glag.abs[kk]) + discontinuity^2
      psi <- mymu + sqrt(2) * sigma * abscissae
      temp9 <- dpsi.dlambda.yjn(psi, lambda, mymu, sigma, derivative = 2)
      temp9 <- cbind(temp9, 
                    1 / (2 * sqrt((abscissae-discontinuity^2)^2 +
                    discontinuity^2) *
                    sqrt(pi) * sigma^2))
      temp7 <- temp1 * glag.wts[kk]
      wz[,iam(1, 1, M)] <- wz[,iam(1, 1, M)] + temp7 * 
          glag.weight.yjn.11(abscissae, lambda, mymu, sigma, temp9)
      wz[,iam(1, 2, M)] <- wz[,iam(1, 2, M)] + temp7 * 
          glag.weight.yjn.12(abscissae, lambda, mymu, sigma, temp9)
      wz[,iam(1, 3, M)] <- wz[,iam(1, 3, M)] + temp7 * 
          glag.weight.yjn.13(abscissae, lambda, mymu, sigma, temp9)
    }

    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] * dlambda.deta^2
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)] * dlambda.deta
    wz[, iam(1, 3, M)] <- wz[, iam(1, 3, M)] * dsigma.deta * dlambda.deta
    if ( .diagW && iter <= .iters.diagW) {
      wz[,iam(1, 2, M)] <- wz[, iam(1, 3, M)] <- 0
    }
    wz[, iam(2, 3, M)] <- wz[, iam(2, 3, M)] * dsigma.deta
    wz[, iam(3, 3, M)] <- wz[, iam(3, 3, M)] * dsigma.deta^2

        c(w) * wz
  }), list(.lsigma = lsigma,
           .esigma = esigma, .elambda = elambda,
           .rule = rule,
           .diagW = diagW,
           .iters.diagW = iters.diagW,
           .llambda = llambda))))
}



lmscreg.control <- function(cdf = TRUE, at.arg = NULL, x0 = NULL, ...) {

  if (!is.logical(cdf)) {
    warning("'cdf' is not logical; using TRUE instead")
    cdf <- TRUE
  }
  list(cdf = cdf, at.arg = at.arg, x0 = x0)
}





Wr1 <- function(r, w) ifelse(r <= 0, 1, w)


Wr2 <- function(r, w) (r <= 0) * 1 + (r > 0) * w


amlnormal.deviance <- function(mu, y, w, residuals = FALSE,
                               eta, extra = NULL) {

  M <- length(extra$w.aml)

  if (M > 1) y <- matrix(y, extra$n, extra$M)

  devi <-  cbind((y - mu)^2)
  if (residuals) {
    stop("not sure here")
    wz <- VGAM.weights.function(w = w, M = extra$M, n = extra$n)
    return((y - mu) * sqrt(wz) * matrix(extra$w.aml,extra$n,extra$M))
  } else {
    all.deviances <- numeric(M)
    myresid <- matrix(y,extra$n,extra$M) - cbind(mu)
    for (ii in 1:M)
      all.deviances[ii] <- sum(c(w) * devi[, ii] *
                               Wr1(myresid[, ii],
                                   w = extra$w.aml[ii]))
  }
  if (is.logical(extra$individual) && extra$individual)
    all.deviances else sum(all.deviances)
}



 amlnormal <- function(w.aml = 1, parallel = FALSE,
                       lexpectile = "identitylink",
                       iexpectile = NULL,
                       imethod = 1, digw = 4) {




  if (!is.Numeric(w.aml, positive = TRUE))
    stop("argument 'w.aml' must be a vector of positive values")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1, 2 or 3")



  lexpectile <- as.list(substitute(lexpectile))
  eexpectile <- link2list(lexpectile)
  lexpectile <- attr(eexpectile, "function.name")


  if (length(iexpectile) && !is.Numeric(iexpectile))
    stop("bad input for argument 'iexpectile'")

  new("vglmff",
  blurb = c("Asymmetric least squares quantile regression\n\n",
            "Links:    ",
            namesof("expectile", link = lexpectile, earg = eexpectile)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
  }), list( .parallel = parallel ))),
  deviance = function(mu, y, w, residuals = FALSE, eta,
                      extra = NULL) {
    amlnormal.deviance(mu = mu, y = y, w = w, residuals = residuals,
                       eta = eta, extra = extra)
  },
  initialize = eval(substitute(expression({
    extra$w.aml <- .w.aml

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$M <- M <- length(extra$w.aml)  # Recycle if necessary
    extra$n <- n
    extra$y.names <- y.names <-
      paste("w.aml = ", round(extra$w.aml, digits = .digw ), sep = "")

    predictors.names <- c(namesof(
        paste("expectile(",y.names,")", sep = ""), .lexpectile ,
               earg = .eexpectile, tag = FALSE))

    if (!length(etastart)) {
      mean.init <-
        if ( .imethod == 1)
          rep(median(y), length = n) else
        if ( .imethod == 2 || .imethod == 3)
          rep(weighted.mean(y, w), length = n) else {
              junk <- lm.wfit(x = x, y = c(y), w = c(w))
              junk$fitted
        }


        if ( .imethod == 3)
          mean.init <- abs(mean.init) + 0.01


        if (length( .iexpectile))
          mean.init <- matrix( .iexpectile, n, M, byrow = TRUE)
        etastart <-
          matrix(theta2eta(mean.init, .lexpectile,
                           earg = .eexpectile), n, M)
    }
  }), list( .lexpectile = lexpectile, .eexpectile = eexpectile,
            .iexpectile = iexpectile,
            .imethod = imethod, .digw = digw, .w.aml = w.aml ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    ans <- eta <- as.matrix(eta)
    for (ii in 1:ncol(eta))
      ans[, ii] <- eta2theta(eta[, ii], .lexpectile, earg = .eexpectile)
    dimnames(ans) <- list(dimnames(eta)[[1]], extra$y.names)
    ans
  }, list( .lexpectile = lexpectile, .eexpectile = eexpectile ))),
  last = eval(substitute(expression({
    misc$link <- rep(.lexpectile , length = M)
    names(misc$link) <- extra$y.names

    misc$earg <- vector("list", M)
    for (ilocal in 1:M)
      misc$earg[[ilocal]] <- list(theta = NULL)
    names(misc$earg) <- names(misc$link)

    misc$parallel <- .parallel
    misc$expected <- TRUE
    extra$percentile <- numeric(M)  # These are estimates (empirical)
    misc$multipleResponses <- TRUE


    for (ii in 1:M) {
      use.w <- if (M > 1 && ncol(cbind(w)) == M) w[, ii] else w
      extra$percentile[ii] <- 100 *
        weighted.mean(myresid[, ii] <= 0, use.w)
    }
    names(extra$percentile) <- names(misc$link)

    extra$individual <- TRUE
    if (!(M > 1 && ncol(cbind(w)) == M)) {
      extra$deviance <-
        amlnormal.deviance(mu = mu, y = y, w = w,
                           residuals = FALSE, eta = eta, extra = extra)
      names(extra$deviance) <- extra$y.names
    }
  }), list( .lexpectile = lexpectile,
            .eexpectile = eexpectile, .parallel = parallel ))),
  vfamily = c("amlnormal"),

  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta, .lexpectile, earg = .eexpectile)
    dexpectile.deta <- dtheta.deta(mymu, .lexpectile, earg = .eexpectile)
    myresid <- matrix(y,extra$n,extra$M) - cbind(mu)
    wor1 <- Wr2(myresid, w = matrix(extra$w.aml, extra$n, extra$M,
                                   byrow = TRUE))
    c(w) * myresid * wor1 * dexpectile.deta
  }), list( .lexpectile = lexpectile,
            .eexpectile = eexpectile ))),

  weight = eval(substitute(expression({
    wz <- c(w) * wor1 * dexpectile.deta^2
    wz
  }), list( .lexpectile = lexpectile,
            .eexpectile = eexpectile ))))
}










amlpoisson.deviance <- function(mu, y, w, residuals = FALSE, eta,
                                extra = NULL) {

    M <- length(extra$w.aml)

    if (M > 1) y <- matrix(y,extra$n,extra$M)

    nz <- y > 0
    devi <-  cbind(-(y - mu))
    devi[nz] <- devi[nz] + y[nz] * log(y[nz]/mu[nz])
    if (residuals) {
        stop("not sure here")
        return(sign(y - mu) * sqrt(2 * abs(devi) * w) *
               matrix(extra$w,extra$n,extra$M))
    } else {
        all.deviances <- numeric(M)
        myresid <- matrix(y,extra$n,extra$M) - cbind(mu)
        for (ii in 1:M) all.deviances[ii] <- 2 * sum(c(w) * devi[, ii] *
                               Wr1(myresid[, ii], w=extra$w.aml[ii]))
    }
    if (is.logical(extra$individual) && extra$individual)
        all.deviances else sum(all.deviances)
}


 amlpoisson <- function(w.aml = 1, parallel = FALSE, imethod = 1,
                        digw = 4, link = "loge") {
  if (!is.Numeric(w.aml, positive = TRUE))
    stop("'w.aml' must be a vector of positive values")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Poisson expectile regression by",
            " asymmetric maximum likelihood estimation\n\n",
            "Link:     ", namesof("expectile", link, earg = earg)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
  }), list( .parallel = parallel ))),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    amlpoisson.deviance(mu = mu, y = y, w = w, residuals = residuals,
                        eta = eta, extra = extra)
  },
  initialize = eval(substitute(expression({
    extra$w.aml <- .w.aml

    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1, ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$M <- M <- length(extra$w.aml)  # Recycle if necessary
    extra$n <- n
    extra$y.names <- y.names <-
      paste("w.aml = ", round(extra$w.aml, digits = .digw ), sep = "")
    extra$individual <- FALSE
    predictors.names <-
      c(namesof(paste("expectile(",y.names,")", sep = ""),
                .link , earg = .earg , tag = FALSE))

    if (!length(etastart)) {
        mean.init <- if ( .imethod == 2)
              rep(median(y), length = n) else
            if ( .imethod == 1)
                rep(weighted.mean(y, w), length = n) else {
                    junk = lm.wfit(x = x, y = c(y), w = c(w))
                    abs(junk$fitted)
                }
        etastart <-
          matrix(theta2eta(mean.init, .link , earg = .earg ), n, M)
    }
  }), list( .link = link, .earg = earg, .imethod = imethod,
            .digw = digw, .w.aml = w.aml ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu.ans <- eta <- as.matrix(eta)
    for (ii in 1:ncol(eta))
      mu.ans[, ii] <- eta2theta(eta[, ii], .link , earg = .earg )
    dimnames(mu.ans) <- list(dimnames(eta)[[1]], extra$y.names)
    mu.ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$multipleResponses <- TRUE
    misc$expected <- TRUE
    misc$parallel <- .parallel


    misc$link <- rep(.link , length = M)
    names(misc$link) <- extra$y.names

    misc$earg <- vector("list", M)
    for (ilocal in 1:M)
      misc$earg[[ilocal]] <- list(theta = NULL)
    names(misc$earg) <- names(misc$link)

    extra$percentile <- numeric(M)  # These are estimates (empirical)
    for (ii in 1:M)
      extra$percentile[ii] <- 100 * weighted.mean(myresid[, ii] <= 0, w)
    names(extra$percentile) <- names(misc$link)

    extra$individual <- TRUE
    extra$deviance <- amlpoisson.deviance(mu = mu, y = y, w = w,
                                          residuals = FALSE,
                                          eta = eta, extra = extra)
    names(extra$deviance) <- extra$y.names
  }), list( .link = link, .earg = earg, .parallel = parallel ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link =  .link , earg = .earg )
  }, list( .link = link, .earg = earg ))),
  vfamily = c("amlpoisson"),
  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta, .link , earg = .earg )
    dexpectile.deta <- dtheta.deta(mymu, .link , earg = .earg )
    myresid <- matrix(y,extra$n,extra$M) - cbind(mu)
    wor1 <- Wr2(myresid, w = matrix(extra$w.aml, extra$n, extra$M,
                                    byrow = TRUE))
    c(w) * myresid * wor1 * (dexpectile.deta / mymu)
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    use.mu <- mymu
    use.mu[use.mu < .Machine$double.eps^(3/4)] <- .Machine$double.eps^(3/4)
    wz <- c(w) * wor1 * use.mu * (dexpectile.deta / mymu)^2
    wz
  }), list( .link = link, .earg = earg ))))
}





amlbinomial.deviance <- function(mu, y, w, residuals = FALSE,
                                 eta, extra = NULL) {

    M <- length(extra$w.aml)


    if (M > 1) y <- matrix(y,extra$n,extra$M)


    devy <- y
    nz <- y != 0
    devy[nz] <- y[nz] * log(y[nz])
    nz <- (1 - y) != 0
    devy[nz] <- devy[nz] + (1 - y[nz]) * log1p(-y[nz])
    devmu <- y * log(mu) + (1 - y) * log1p(-mu)
    if (any(small <- mu * (1 - mu) < .Machine$double.eps)) {
      warning("fitted values close to 0 or 1")
      smu <- mu[small]
      sy <- y[small]
      smu <- ifelse(smu < .Machine$double.eps,
                    .Machine$double.eps, smu)
      onemsmu <- ifelse((1 - smu) < .Machine$double.eps,
                        .Machine$double.eps, 1 - smu)
      devmu[small] <- sy * log(smu) + (1 - sy) * log(onemsmu)
    }
    devi <- 2 * (devy - devmu)
    if (residuals) {
      stop("not sure here")
      return(sign(y - mu) * sqrt(abs(devi) * w))
    } else {
      all.deviances <- numeric(M)
      myresid <- matrix(y,extra$n,extra$M) - matrix(mu,extra$n,extra$M)
      for (ii in 1:M) all.deviances[ii] <- sum(c(w) * devi[, ii] *
                             Wr1(myresid[, ii], w=extra$w.aml[ii]))
    }
    if (is.logical(extra$individual) && extra$individual)
      all.deviances else sum(all.deviances)
}


 amlbinomial <- function(w.aml = 1, parallel = FALSE, digw = 4,
                         link = "logit") {

  if (!is.Numeric(w.aml, positive = TRUE))
    stop("'w.aml' must be a vector of positive values")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Logistic expectile regression by ",
            "asymmetric maximum likelihood estimation\n\n",
            "Link:     ", namesof("expectile", link, earg = earg)),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
  }), list( .parallel = parallel ))),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    amlbinomial.deviance(mu = mu, y = y, w = w, residuals = residuals,
                         eta = eta, extra = extra)
  },
  initialize = eval(substitute(expression({


        {
            NCOL <- function (x)
                if (is.array(x) && length(dim(x)) > 1 ||
                is.data.frame(x)) ncol(x) else as.integer(1)

            if (NCOL(y) == 1) {
                if (is.factor(y)) y <- y != levels(y)[1]
                nn <- rep(1, n)
                if (!all(y >= 0 & y <= 1))
                    stop("response values must be in [0, 1]")
                if (!length(mustart) && !length(etastart))
                    mustart <- (0.5 + w * y) / (1 + w)
                no.successes <- w * y
                if (any(abs(no.successes - round(no.successes)) > 0.001))
                    stop("Number of successes must be integer-valued")
            } else if (NCOL(y) == 2) {
                if (any(abs(y - round(y)) > 0.001))
                    stop("Count data must be integer-valued")
                nn <- y[, 1] + y[, 2]
                y <- ifelse(nn > 0, y[, 1]/nn, 0)
                w <- w * nn
                if (!length(mustart) && !length(etastart))
                    mustart <- (0.5 + nn * y) / (1 + nn)
            } else
                 stop("Response not of the right form")
        }

        extra$w.aml <- .w.aml
        if (ncol(y <- cbind(y)) != 1)
            stop("response must be a vector or a one-column matrix")
        extra$M <- M <- length(extra$w.aml)  # Recycle if necessary
        extra$n <- n
        extra$y.names <- y.names <-
            paste("w.aml = ", round(extra$w.aml, digits = .digw ), sep = "")
        extra$individual <- FALSE
        predictors.names <-
            c(namesof(paste("expectile(", y.names, ")", sep = ""),
                      .link , earg = .earg , tag = FALSE))

        if (!length(etastart)) {
          etastart <- matrix(theta2eta(mustart, .link , earg = .earg ), n, M)
          mustart <- NULL
        }


  }), list( .link = link, .earg = earg,
            .digw = digw, .w.aml = w.aml ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu.ans <- eta <- as.matrix(eta)
    for (ii in 1:ncol(eta))
      mu.ans[, ii] <- eta2theta(eta[, ii], .link , earg = .earg )
    dimnames(mu.ans) <- list(dimnames(eta)[[1]], extra$y.names)
    mu.ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- rep(.link , length = M)
    names(misc$link) <- extra$y.names

    misc$earg <- vector("list", M)
    for (ilocal in 1:M)
      misc$earg[[ilocal]] <- list(theta = NULL)
    names(misc$earg) <- names(misc$link)

    misc$parallel <- .parallel
    misc$expected <- TRUE

    extra$percentile <- numeric(M)  # These are estimates (empirical)
    for (ii in 1:M)
      extra$percentile[ii] <- 100 * weighted.mean(myresid[, ii] <= 0, w)
    names(extra$percentile) <- names(misc$link)

    extra$individual <- TRUE
    extra$deviance <- amlbinomial.deviance(mu = mu, y = y, w = w,
                     residuals = FALSE, eta = eta, extra = extra)
    names(extra$deviance) <- extra$y.names
  }), list( .link = link, .earg = earg, .parallel = parallel ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link =  .link , earg = .earg )
  }, list( .link = link, .earg = earg ))),
  vfamily = c("amlbinomial"),
  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta, .link , earg = .earg )
    use.mu <- mymu
    use.mu[use.mu < .Machine$double.eps^(3/4)] = .Machine$double.eps^(3/4)
    dexpectile.deta <- dtheta.deta(use.mu, .link , earg = .earg )
    myresid <- matrix(y,extra$n,extra$M) - cbind(mu)
    wor1 <- Wr2(myresid, w = matrix(extra$w.aml, extra$n, extra$M,
                                   byrow = TRUE))
    c(w) * myresid * wor1 * (dexpectile.deta / (use.mu * (1-use.mu)))
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    wz <- c(w) * wor1 * (dexpectile.deta^2 / (use.mu * (1 - use.mu)))
    wz
  }), list( .link = link, .earg = earg))))
}










amlexponential.deviance <- function(mu, y, w, residuals = FALSE,
                                   eta, extra = NULL) {

  M <- length(extra$w.aml)

  if (M > 1) y <- matrix(y,extra$n,extra$M)

  devy <-  cbind(-log(y) - 1)
  devi <-  cbind(-log(mu) - y / mu)
  if (residuals) {
    stop("not sure here")
    return(sign(y - mu) * sqrt(2 * abs(devi) * w) *
           matrix(extra$w,extra$n,extra$M))
  } else {
    all.deviances <- numeric(M)
    myresid <- matrix(y,extra$n,extra$M) - cbind(mu)
    for (ii in 1:M) all.deviances[ii] = 2 * sum(c(w) *
                           (devy[, ii] - devi[, ii]) *
                           Wr1(myresid[, ii], w=extra$w.aml[ii]))
  }
  if (is.logical(extra$individual) && extra$individual)
    all.deviances else sum(all.deviances)
}




 amlexponential <- function(w.aml = 1, parallel = FALSE, imethod = 1,
                            digw = 4, link = "loge") {
  if (!is.Numeric(w.aml, positive = TRUE))
    stop("'w.aml' must be a vector of positive values")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1, 2 or 3")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  y.names <- paste("w.aml = ", round(w.aml, digits = digw), sep = "")
  predictors.names <- c(namesof(
      paste("expectile(", y.names,")", sep = ""), link, earg = earg))
  predictors.names <- paste(predictors.names, collapse = ", ")


  new("vglmff",
  blurb = c("Exponential expectile regression by",
            " asymmetric maximum likelihood estimation\n\n",
            "Link:     ", predictors.names),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel ,
                           constraints = constraints)
  }), list( .parallel = parallel ))),
  deviance = function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    amlexponential.deviance(mu = mu, y = y, w = w,
                            residuals = residuals,
                            eta = eta, extra = extra)
  },
  initialize = eval(substitute(expression({
    extra$w.aml <- .w.aml


    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              ncol.w.max = 1, ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$M <- M <- length(extra$w.aml)  # Recycle if necessary
    extra$n <- n
    extra$y.names <- y.names <-
        paste("w.aml = ", round(extra$w.aml, digits = .digw ), sep = "")
    extra$individual = FALSE


    predictors.names <- c(namesof(
        paste("expectile(", y.names, ")", sep = ""),
        .link , earg = .earg , tag = FALSE))

    if (!length(etastart)) {
      mean.init <- if ( .imethod == 1)
              rep(median(y), length = n) else
          if ( .imethod == 2)
              rep(weighted.mean(y, w), length = n) else {
                  1 / (y + 1)
              }
      etastart <- matrix(theta2eta(mean.init, .link , earg = .earg ),
                        n, M)
    }
  }), list( .link = link, .earg = earg, .imethod = imethod,
            .digw = digw, .w.aml = w.aml ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu.ans <- eta <- as.matrix(eta)
    for (ii in 1:ncol(eta))
      mu.ans[, ii] <- eta2theta(eta[, ii], .link , earg = .earg )
    dimnames(mu.ans) <- list(dimnames(eta)[[1]], extra$y.names)
    mu.ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$multipleResponses <- TRUE
    misc$expected <- TRUE
    misc$parallel <- .parallel

    misc$link <- rep(.link , length = M)
    names(misc$link) <- extra$y.names

    misc$earg <- vector("list", M)
    for (ilocal in 1:M)
      misc$earg[[ilocal]] <- list(theta = NULL)
    names(misc$earg) <- names(misc$link)


    extra$percentile <- numeric(M)  # These are estimates (empirical)
    for (ii in 1:M)
      extra$percentile[ii] <- 100 * weighted.mean(myresid[, ii] <= 0, w)
    names(extra$percentile) <- names(misc$link)

    extra$individual <- TRUE
    extra$deviance =
      amlexponential.deviance(mu = mu, y = y, w = w,
                              residuals = FALSE, eta = eta, extra = extra)
    names(extra$deviance) <- extra$y.names
  }), list( .link = link, .earg = earg, .parallel = parallel ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link =  .link , earg = .earg )
  }, list( .link = link, .earg = earg ))),
  vfamily = c("amlexponential"),
  deriv = eval(substitute(expression({
    mymu <- eta2theta(eta, .link , earg = .earg )
    bigy <- matrix(y,extra$n,extra$M)
    dl.dmu <- (bigy - mymu) / mymu^2

    dmu.deta <- dtheta.deta(mymu, .link , earg = .earg )
    myresid <- bigy - cbind(mymu)
    wor1 <- Wr2(myresid, w = matrix(extra$w.aml, extra$n, extra$M,
                                   byrow = TRUE))
    c(w) * wor1 * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    ned2l.dmu2 <- 1 / mymu^2
    wz <- c(w) * wor1 * ned2l.dmu2 * dmu.deta^2
    wz
  }), list( .link = link, .earg = earg ))))
}






rho1check <- function(u, tau = 0.5)
  u * (tau - (u <= 0))




dalap <- function(x, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau/(1-tau)), log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)




  NN <- max(length(x), length(location), length(scale), length(kappa),
            length(tau))
  if (length(x)        != NN) x        <- rep(x,        length.out = NN)
  if (length(location) != NN) location <- rep(location, length.out = NN)
  if (length(scale)    != NN) scale    <- rep(scale,    length.out = NN)
  if (length(kappa)    != NN) kappa    <- rep(kappa,    length.out = NN)
  if (length(tau)      != NN) tau      <- rep(tau,      length.out = NN)

  logconst <- 0.5 * log(2) - log(scale) + log(kappa) - log1p(kappa^2)
  exponent <- -(sqrt(2) / scale) * abs(x - location) *
             ifelse(x >= location, kappa, 1/kappa)

  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  logconst[!indexTF] <- NaN

  if (log.arg) logconst + exponent else exp(logconst + exponent)
}


ralap <- function(n, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau/(1-tau))) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  location <- rep(location, length.out = use.n);
  scale    <- rep(scale,    length.out = use.n)
  tau      <- rep(tau,      length.out = use.n);
  kappa    <- rep(kappa,    length.out = use.n);
  ans <- location + scale *
        log(runif(use.n)^kappa / runif(use.n)^(1/kappa)) / sqrt(2)
  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}



palap <- function(q, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau/(1-tau)),
                  lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  NN <- max(length(q), length(location), length(scale), length(kappa),
            length(tau))
  if (length(q)        != NN) q        <- rep(q,        length.out = NN)
  if (length(location) != NN) location <- rep(location, length.out = NN)
  if (length(scale)    != NN) scale    <- rep(scale,    length.out = NN)
  if (length(kappa)    != NN) kappa    <- rep(kappa,    length.out = NN)
  if (length(tau)      != NN) tau      <- rep(tau,      length.out = NN)

  exponent <- -(sqrt(2) / scale) * abs(q - location) *
              ifelse(q >= location, kappa, 1/kappa)
  temp5 <- exp(exponent) / (1 + kappa^2)
  index1 <- (q < location)


  if (lower.tail) {
    if (log.p) {
      ans <- log1p(-exp(exponent) / (1 + kappa^2))
      logtemp5 <- exponent - log1p(kappa^2)
      ans[index1] <- 2 * log(kappa[index1]) + logtemp5[index1]
    } else {
      ans <- (kappa^2 - expm1(exponent)) / (1 + kappa^2)
      ans[index1] <- (kappa[index1])^2 * temp5[index1]
    }
  } else {
    if (log.p) {
      ans <- exponent - log1p(kappa^2)  # logtemp5
      ans[index1] <- log1p(-(kappa[index1])^2 * temp5[index1])
    } else {
      ans <- temp5
      ans[index1] <- (1 + (kappa[index1])^2 *
                     (-expm1(exponent[index1]))) / (1+(kappa[index1])^2)
      }
  } 
  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}



qalap <- function(p, location = 0, scale = 1, tau = 0.5,
                  kappa = sqrt(tau / (1 - tau)),
                  lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  NN <- max(length(p), length(location), length(scale), length(kappa),
            length(tau))
  if (length(p)        != NN) p        <- rep(p,        length.out = NN)
  if (length(location) != NN) location <- rep(location, length.out = NN)
  if (length(scale)    != NN) scale    <- rep(scale,    length.out = NN)
  if (length(kappa)    != NN) kappa    <- rep(kappa,    length.out = NN)
  if (length(tau)      != NN) tau      <- rep(tau,      length.out = NN)



  temp5 <- kappa^2 / (1 + kappa^2)
  if (lower.tail) {
    if (log.p) {
      ans <- exp(p) 
      index1 <- (exp(p) <= temp5)
      exponent <- exp(p[index1]) / temp5[index1]
      ans[index1] <- location[index1] + (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] - (scale[!index1] / kappa[!index1]) *
        (log1p((kappa[!index1])^2) +
           log(-expm1(p[!index1]))) / sqrt(2)
    } else {
      ans <- p 
      index1 <- (p <= temp5)
      exponent <- p[index1] / temp5[index1]
      ans[index1] <- location[index1] + (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] - (scale[!index1] / kappa[!index1]) *
                      (log1p((kappa[!index1])^2) +
                      log1p(-p[!index1])) / sqrt(2)
    }
  } else {
    if (log.p) {
      ans <- -expm1(p) 
      index1 <- (-expm1(p)  <= temp5)
      exponent <- -expm1(p[index1]) / temp5[index1]
      ans[index1] <- location[index1] + (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] - (scale[!index1] / kappa[!index1]) *
        (log1p((kappa[!index1])^2) +
           p[!index1]) / sqrt(2)
    } else {
      ans <- exp(log1p(-p)) 
      index1 <- (p >= (1 / (1+kappa^2))) 
      exponent <- exp(log1p(-p[index1])) / temp5[index1]
      ans[index1] <- location[index1] + (scale[index1] * kappa[index1]) *
        log(exponent) / sqrt(2)
      ans[!index1] <- location[!index1] - (scale[!index1] / kappa[!index1]) *
        (log1p((kappa[!index1])^2) +
           log(p[!index1])) / sqrt(2)
    }
  }

  indexTF <- (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans  
}






rloglap <- function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                    kappa = sqrt(tau/(1-tau))) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
            stop("bad input for argument 'n'") else n
  location.ald <- rep(location.ald, length.out = use.n);
  scale.ald    <- rep(scale.ald,    length.out = use.n);
  tau          <- rep(tau,          length.out = use.n);
  kappa        <- rep(kappa,        length.out = use.n);
  ans <- exp(location.ald) *
     (runif(use.n)^kappa / runif(use.n)^(1/kappa))^(scale.ald / sqrt(2))
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}



dloglap <- function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
                    kappa = sqrt(tau/(1-tau)), log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  scale    <- scale.ald
  location <- location.ald
  NN <- max(length(x), length(location),
           length(scale), length(kappa), length(tau))

  if (length(x)        != NN) x        <- rep(x,        length.out = NN)
  if (length(location) != NN) location <- rep(location, length.out = NN)
  if (length(scale)    != NN) scale    <- rep(scale,    length.out = NN)
  if (length(kappa)    != NN) kappa    <- rep(kappa,    length.out = NN)
  if (length(tau)      != NN) tau      <- rep(tau,      length.out = NN)


  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- exp(location.ald)
  exponent <- ifelse(x >= Delta, -(Alpha+1), (Beta-1)) *
             (log(x) - location.ald)
  logdensity <- -location.ald + log(Alpha) + log(Beta) -
               log(Alpha + Beta) + exponent
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  if (log.arg) logdensity else exp(logdensity)
}



qloglap <- function(p, location.ald = 0, scale.ald = 1,
                    tau = 0.5, kappa = sqrt(tau/(1-tau)),
                    lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  NN <- max(length(p), length(location.ald), length(scale.ald),
            length(kappa))
  p        <- rep(p,            length.out = NN)
  location <- rep(location.ald, length.out = NN)
  scale    <- rep(scale.ald,    length.out = NN)
  kappa    <- rep(kappa,        length.out = NN)
  tau      <- rep(tau,          length.out = NN)


  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- exp(location.ald)
  temp9 <- Alpha + Beta


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- ifelse((exp(ln.p) > Alpha / temp9), 
                    Delta * (-expm1(ln.p) * temp9 / Beta)^(-1/Alpha),
                    Delta * (exp(ln.p) * temp9 / Alpha)^(1/Beta))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- ifelse((p > Alpha / temp9), 
                    Delta * exp((-1/Alpha) * (log1p(-p) + log(temp9/Beta))),
                    Delta * (p * temp9 / Alpha)^(1/Beta))
      ans[p <  0] <- NaN
      ans[p == 0] <- 0
      ans[p == 1] <- Inf
      ans[p >  1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- ifelse((-expm1(ln.p) > Alpha / temp9), 
                    Delta * (exp(ln.p) * temp9 / Beta)^(-1/Alpha),
                    Delta * (-expm1(ln.p) * temp9 / Alpha)^(1/Beta))
      ans[ln.p > 0] <- NaN
    } else { 
      ans <- ifelse((p < (temp9 - Alpha) / temp9), 
                    Delta * (p * temp9 / Beta)^(-1/Alpha),
                    Delta * exp((1/Beta)*(log1p(-p) + log(temp9/Alpha))))
      ans[p <  0] <- NaN
      ans[p == 0] <- Inf
      ans[p == 1] <- 0
      ans[p >  1] <- NaN
    }
  }
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)
  ans[!indexTF] <- NaN
  ans
}



ploglap <- function(q, location.ald = 0, scale.ald = 1,
                    tau = 0.5, kappa = sqrt(tau/(1-tau)),
                    lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  NN <- max(length(q), length(location.ald), length(scale.ald),
            length(kappa))
  location <- rep(location.ald, length.out = NN)
  scale    <- rep(scale.ald,    length.out = NN)
  kappa    <- rep(kappa,        length.out = NN)
  q        <- rep(q,            length.out = NN)
  tau      <- rep(tau,          length.out = NN)

  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- exp(location.ald)

  temp9 <- Alpha + Beta
  index1 <- (Delta <= q)


  if (lower.tail) {
    if (log.p) {
      ans <- log((Alpha / temp9) * (q / Delta)^(Beta))
      ans[index1] <- log1p((-(Beta/temp9) * (Delta/q)^(Alpha))[index1])
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- (Alpha / temp9) * (q / Delta)^(Beta)
      ans[index1] <- -expm1((log(Beta/temp9) + Alpha * log(Delta/q)))[index1]
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log1p(-(Alpha / temp9) * (q / Delta)^(Beta))
      ans[index1] <- log(((Beta/temp9) * (Delta/q)^(Alpha))[index1])
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- -expm1(log(Alpha/temp9) + Beta * log(q/Delta))
      ans[index1] <- ((Beta/temp9) * (Delta/q)^(Alpha))[index1]
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  } 

  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  ans[!indexTF] <- NaN
  ans
}




rlogitlap <- function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                      kappa = sqrt(tau/(1-tau))) {
  logit(ralap(n = n, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa),
        inverse = TRUE)  # earg = earg
}



dlogitlap <- function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
                      kappa = sqrt(tau/(1-tau)), log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  NN <- max(length(x), length(location.ald),
           length(scale.ald), length(kappa))
  location <- rep(location.ald, length.out = NN);
  scale <- rep(scale.ald, length.out = NN)
  kappa <- rep(kappa, length.out = NN);
  x <- rep(x, length.out = NN)
  tau <- rep(tau, length.out = NN)

  Alpha <- sqrt(2) * kappa / scale.ald
  Beta  <- sqrt(2) / (scale.ald * kappa)
  Delta <- logit(location.ald, inverse = TRUE)  # earg = earg

  exponent <- ifelse(x >= Delta, -Alpha, Beta) *
             (logit(x) - # earg = earg
              location.ald)
  logdensity <- log(Alpha) + log(Beta) - log(Alpha + Beta) -
               log(x) - log1p(-x) + exponent
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  logdensity[x >  1 & indexTF] <- -Inf
  if (log.arg) logdensity else exp(logdensity)
}



qlogitlap <- function(p, location.ald = 0, scale.ald = 1,
                      tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  qqq <- qalap(p = p, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa)
  ans <- logit(qqq, inverse = TRUE)  # earg = earg
  ans[(p < 0) | (p > 1)] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- 1
  ans
}



plogitlap <- function(q, location.ald = 0, scale.ald = 1,
                      tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  NN <- max(length(q), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep(location.ald, length.out = NN);
  scale.ald <- rep(scale.ald, length.out = NN)
  kappa <- rep(kappa, length.out = NN); q <- rep(q, length.out = NN)
  tau <- rep(tau, length.out = NN);

  indexTF <- (q > 0) & (q < 1)
  qqq <- logit(q[indexTF])  # earg = earg
  ans <- q
  ans[indexTF] <- palap(q = qqq, location = location.ald[indexTF],
                       scale = scale.ald[indexTF],
                       tau = tau[indexTF], kappa = kappa[indexTF])
  ans[q >= 1] <- 1
  ans[q <= 0] <- 0
  ans
}




rprobitlap <- function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                       kappa = sqrt(tau/(1-tau))) {



  probit(ralap(n = n, location = location.ald, scale = scale.ald,
               tau = tau, kappa = kappa),
               inverse = TRUE)
}



dprobitlap <-
  function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
           kappa = sqrt(tau/(1-tau)), log = FALSE,
           meth2 = TRUE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  NN <- max(length(x), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep(location.ald, length.out = NN);
  scale.ald <- rep(scale.ald, length.out = NN)
  kappa <- rep(kappa, length.out = NN); x = rep(x, length.out = NN)
  tau <- rep(tau, length.out = NN)

  logdensity <- x * NaN
  index1 <- (x > 0) & (x < 1)
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  if (meth2) {
    dx.dy <- x
    use.x <- probit(x[index1])  # earg = earg
    logdensity[index1] <-
      dalap(x = use.x, location = location.ald[index1],
            scale = scale.ald[index1], tau = tau[index1],
            kappa = kappa[index1], log = TRUE)
  } else {
    Alpha <- sqrt(2) * kappa / scale.ald
    Beta  <- sqrt(2) / (scale.ald * kappa)
    Delta <- pnorm(location.ald)
    use.x  <- qnorm(x)  # qnorm(x[index1])
    log.dy.dw <- dnorm(use.x, log = TRUE)

    exponent <- ifelse(x >= Delta, -Alpha, Beta) *
                     (use.x - location.ald) - log.dy.dw

    logdensity[index1] <- (log(Alpha) + log(Beta) -
                          log(Alpha + Beta) + exponent)[index1]
  }
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  logdensity[x >  1 & indexTF] <- -Inf

  if (meth2) {
    dx.dy[index1] <- probit(x[index1],  # earg = earg,
                            inverse = TRUE,
                            deriv = 1)
    dx.dy[!index1] <- 0
    dx.dy[!indexTF] <- NaN
    if (log.arg) logdensity - log(abs(dx.dy)) else
                 exp(logdensity) / abs(dx.dy)
  } else {
    if (log.arg) logdensity else exp(logdensity)
  }
}


qprobitlap <- function(p, location.ald = 0, scale.ald = 1,
                       tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  qqq <- qalap(p = p, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa)
  ans <- probit(qqq, inverse = TRUE)  # , earg = earg
  ans[(p < 0) | (p > 1)] = NaN
  ans[p == 0] <- 0
  ans[p == 1] <- 1
  ans
}



pprobitlap <- function(q, location.ald = 0, scale.ald = 1,
                       tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  NN <- max(length(q), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep(location.ald, length.out = NN);
  scale.ald <- rep(scale.ald, length.out = NN)
  kappa <- rep(kappa, length.out = NN);
  q <- rep(q, length.out = NN)
  tau <- rep(tau, length.out = NN);

  indexTF <- (q > 0) & (q < 1)
  qqq <- probit(q[indexTF])  # earg = earg
  ans <- q
  ans[indexTF] <- palap(q = qqq, location = location.ald[indexTF],
                       scale = scale.ald[indexTF],
                       tau = tau[indexTF], kappa = kappa[indexTF])
  ans[q >= 1] <- 1
  ans[q <= 0] <- 0
  ans
}





rclogloglap <- function(n, location.ald = 0, scale.ald = 1, tau = 0.5,
                        kappa = sqrt(tau/(1-tau))) {
  cloglog(ralap(n = n, location = location.ald, scale = scale.ald,
                tau = tau, kappa = kappa),  # earg = earg,
          inverse = TRUE)
}



dclogloglap <- function(x, location.ald = 0, scale.ald = 1, tau = 0.5,
                        kappa = sqrt(tau/(1-tau)), log = FALSE,
                        meth2 = TRUE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  NN <- max(length(x), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep(location.ald, length.out = NN)
  scale.ald <- rep(scale.ald, length.out = NN)
  kappa <- rep(kappa, length.out = NN)
  x <- rep(x, length.out = NN)
  tau <- rep(tau, length.out = NN)

  logdensity <- x * NaN
  index1 <- (x > 0) & (x < 1)
  indexTF <- (scale.ald > 0) & (tau > 0) & (tau < 1) & (kappa > 0)  # &
  if (meth2) {
    dx.dy <- x
    use.w <- cloglog(x[index1])  # earg = earg
    logdensity[index1] <-
      dalap(x = use.w, location = location.ald[index1],
            scale = scale.ald[index1],
            tau = tau[index1],
            kappa = kappa[index1], log = TRUE)

  } else {
    Alpha <- sqrt(2) * kappa / scale.ald
    Beta  <- sqrt(2) / (scale.ald * kappa)
    Delta <- cloglog(location.ald, inverse = TRUE)

    exponent <- ifelse(x >= Delta, -(Alpha+1), Beta-1) * log(-log1p(-x)) +
               ifelse(x >= Delta, Alpha, -Beta) * location.ald
    logdensity[index1] <- (log(Alpha) + log(Beta) -
                     log(Alpha + Beta) - log1p(-x) + exponent)[index1]
  }
  logdensity[!indexTF] <- NaN
  logdensity[x <  0 & indexTF] <- -Inf
  logdensity[x >  1 & indexTF] <- -Inf

  if (meth2) {
    dx.dy[index1] <- cloglog(x[index1],  # earg = earg,
                             inverse = TRUE, deriv = 1)
    dx.dy[!index1] <- 0
    dx.dy[!indexTF] <- NaN
    if (log.arg) logdensity - log(abs(dx.dy)) else
                 exp(logdensity) / abs(dx.dy)
  } else {
    if (log.arg) logdensity else exp(logdensity)
  }
}



qclogloglap <- function(p, location.ald = 0, scale.ald = 1,
                       tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  qqq <- qalap(p = p, location = location.ald, scale = scale.ald,
              tau = tau, kappa = kappa)
  ans <- cloglog(qqq, inverse = TRUE)  # , earg = earg
  ans[(p < 0) | (p > 1)] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- 1
  ans
}



pclogloglap <- function(q, location.ald = 0, scale.ald = 1,
                       tau = 0.5, kappa = sqrt(tau/(1-tau))) {
  NN <- max(length(q), length(location.ald), length(scale.ald),
           length(kappa))
  location.ald <- rep(location.ald, length.out = NN);
  scale.ald <- rep(scale.ald, length.out = NN)
  kappa <- rep(kappa, length.out = NN);
  q <- rep(q, length.out = NN)
  tau <- rep(tau, length.out = NN);

  indexTF <- (q > 0) & (q < 1)
  qqq <- cloglog(q[indexTF])  # earg = earg
  ans <- q
  ans[indexTF] <- palap(q = qqq, location = location.ald[indexTF],
                       scale = scale.ald[indexTF],
                       tau = tau[indexTF], kappa = kappa[indexTF])
  ans[q >= 1] <- 1
  ans[q <= 0] <- 0
  ans
}












alaplace2.control <- function(maxit = 100, ...) {
  list(maxit = maxit)
}


 alaplace2 <-
  function(tau = NULL,
           llocation = "identitylink", lscale = "loge",
           ilocation = NULL,           iscale = NULL,
           kappa = sqrt(tau / (1-tau)),
           ishrinkage = 0.95,

           parallel.locat = TRUE  ~ 0,
           parallel.scale = FALSE ~ 0,

           digt = 4,
           idf.mu = 3,
           imethod = 1,
           zero = "scale") {



  apply.parint.locat <- FALSE
  apply.parint.scale <- TRUE




  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  ilocat <- ilocation



  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")
  if (!is.Numeric(ishrinkage, length.arg = 1) ||
    ishrinkage < 0 ||
    ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  if (length(tau) &&
      max(abs(kappa - sqrt(tau / (1 - tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")





  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")





  new("vglmff",
  blurb = c("Two-parameter asymmetric Laplace distribution\n\n",
            "Links:      ",
            namesof("location1", llocat, earg = elocat), ", ",
            namesof("scale1",    lscale, earg = escale), ", ",
            namesof("location2", llocat, earg = elocat), ", ",
            namesof("scale2",    lscale, earg = escale),
            ", ..., ",
            "\n\n",
            "Mean:       ",
            "location + scale * (1/kappa - kappa) / sqrt(2)", "\n",
            "Quantiles:  location", "\n",
            "Variance:   scale^2 * (1 + kappa^4) / (2 * kappa^2)"),




  constraints = eval(substitute(expression({
 

    onemat <- matrix(1, Mdiv2, 1)
    constraints.orig <- constraints


    cm1.locat <- kronecker(diag(Mdiv2), rbind(1, 0))
    cmk.locat <- kronecker(onemat,      rbind(1, 0))
    con.locat <- cm.VGAM(cmk.locat,
                         x = x, bool = .parallel.locat ,
                         constraints = constraints.orig,
                         apply.int = .apply.parint.locat ,
                         cm.default           = cm1.locat,
                         cm.intercept.default = cm1.locat)
   
    

    cm1.scale <- kronecker(diag(Mdiv2), rbind(0, 1))
    cmk.scale <- kronecker(onemat,      rbind(0, 1))
    con.scale <- cm.VGAM(cmk.scale,
                         x = x, bool = .parallel.scale ,
                         constraints = constraints.orig,
                         apply.int = .apply.parint.scale ,
                         cm.default           = cm1.scale,
                         cm.intercept.default = cm1.scale)
   
    con.use <- con.scale
    for (klocal in 1:length(con.scale)) {
      con.use[[klocal]] <- cbind(con.locat[[klocal]],
                                 con.scale[[klocal]])
    }

    
    constraints <- con.use

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .parallel.locat = parallel.locat,
            .parallel.scale = parallel.scale,
            .zero = zero,
            .apply.parint.scale = apply.parint.scale,
            .apply.parint.locat = apply.parint.locat ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         summary.pvalues = FALSE,
         multipleResponses = FALSE,
         parameters.names = c("location1", "scale1", "location2", "scale2"),
         zero = .zero )
  }, list( .zero = zero ))),
  initialize = eval(substitute(expression({
    extra$M1 <- M1 <- 2


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = if (length( .kappa ) > 1) 1 else Inf,
              ncol.y.max = if (length( .kappa ) > 1) 1 else Inf,
              out.wy = TRUE,
              colsyperw = 1,  # Uncommented out 20140621
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$ncoly <- ncoly <- ncol(y)
    if ((ncoly > 1) && (length( .kappa ) > 1))
      stop("response must be a vector if 'kappa' or 'tau' ",
           "has a length greater than one")



    extra$kappa <- .kappa
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)

    extra$Mdiv2 <- Mdiv2 <- max(ncoly, length( .kappa ))
    extra$M <- M <- M1 * Mdiv2
    extra$n <- n



    extra$tau.names <- tau.names <-
      paste("(tau = ", round(extra$tau, digits = .digt), ")", sep = "")
    extra$Y.names <- Y.names <- if (ncoly > 1) dimnames(y)[[2]] else "y"
    if (is.null(Y.names) || any(Y.names == ""))
      extra$Y.names <- Y.names <- paste("y", 1:ncoly, sep = "")
    extra$y.names <- y.names <-
      if (ncoly > 1) paste(Y.names, tau.names, sep = "") else tau.names

    extra$individual <- FALSE


    mynames1 <- param.names("location", Mdiv2)
    mynames2 <- param.names("scale",    Mdiv2)
    predictors.names <-
        c(namesof(mynames1, .llocat , earg = .elocat, tag = FALSE),
          namesof(mynames2, .lscale , earg = .escale, tag = FALSE))
    predictors.names <-
    predictors.names[interleave.VGAM(M, M1 = M1)]




    locat.init <- scale.init <- matrix(0, n, Mdiv2)
    if (!length(etastart)) {
      for (jay in 1:Mdiv2) {
        y.use <- if (ncoly > 1) y[, jay] else y
        Jay   <- if (ncoly > 1) jay else 1
        if ( .imethod == 1) {
          locat.init[, jay] <- weighted.mean(y.use, w[, Jay])
          scale.init[, jay] <- sqrt(var(y.use) / 2)
        } else if ( .imethod == 2) {
          locat.init[, jay] <- median(y.use)
          scale.init[, jay] <- sqrt(sum(c(w[, Jay]) *
             abs(y - median(y.use))) / (sum(w[, Jay]) * 2))
        } else if ( .imethod == 3) {
          Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                 y = y.use, w = w[, Jay],
                                 df = .idf.mu )
          locat.init[, jay] <- predict(Fit5, x = x[, min(ncol(x), 2)])$y
          scale.init[, jay] <- sqrt(sum(c(w[, Jay]) *
                                    abs(y.use - median(y.use))) / (
                                        sum(w[, Jay]) * 2))
        } else {
          use.this <- weighted.mean(y.use, w[, Jay])
          locat.init[, jay] <- (1 - .ishrinkage ) * y.use + .ishrinkage * use.this
          scale.init[, jay] <-
            sqrt(sum(c(w[, Jay]) *
            abs(y.use - median(y.use ))) / (sum(w[, Jay]) * 2))
        }
      }



      if (length( .ilocat )) {
        locat.init <- matrix( .ilocat , n, Mdiv2, byrow = TRUE)
      }
      if (length( .iscale )) {
        scale.init <- matrix( .iscale , n, Mdiv2, byrow = TRUE)
      }

      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ),
                theta2eta(scale.init, .lscale , earg = .escale ))
      etastart <- etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
    }
  }), list( .imethod = imethod,
            .idf.mu = idf.mu,
            .ishrinkage = ishrinkage, .digt = digt,
            .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale, .kappa = kappa,
            .ilocat = ilocat, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Mdiv2 <- extra$Mdiv2
    locat <- eta2theta(eta[, 2 * (1:Mdiv2) - 1, drop = FALSE],
                      .llocat , earg = .elocat )
    dimnames(locat) <- list(dimnames(eta)[[1]], extra$y.names)
    myans <- if ( .fittedMean ) {
      kappamat <- matrix(extra$kappa, extra$n, extra$Mdiv2,
                         byrow = TRUE)
      Scale <- eta2theta(eta[, 2 * (1:Mdiv2)    , drop = FALSE],
                         .lscale , earg = .escale )
      locat + Scale * (1/kappamat - kappamat)
    } else {
      locat
    }
    dimnames(myans) <- list(dimnames(myans)[[1]], extra$y.names)
    myans
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .fittedMean = fittedMean,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    M1 <- extra$M1

    tmp34 <- c(rep( .llocat , length = Mdiv2),
               rep( .lscale , length = Mdiv2))
    names(tmp34) <- c(mynames1, mynames2) 
    tmp34 <- tmp34[interleave.VGAM(M, M1 = M1)]
    misc$link <- tmp34  # Already named

    misc$earg <- vector("list", M)
    misc$M1 <- M1
    for (ii in 1:Mdiv2) {
      misc$earg[[M1 * ii - 1]] <- .elocat
      misc$earg[[M1 * ii    ]] <- .escale
    }
    names(misc$earg) <- names(misc$link)


    misc$multipleResponses <- TRUE
    misc$expected <- TRUE
    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    misc$true.mu <- .fittedMean  # @fitted is not a true mu?

    extra$percentile <- numeric(Mdiv2)  # length(misc$kappa)
    locat <- as.matrix(locat)
    for (ii in 1:Mdiv2) {
      y.use <- if (ncoly > 1) y[, ii] else y
      Jay   <- if (ncoly > 1) ii else 1
      extra$percentile[ii] <- 100 * weighted.mean(y.use <= locat[, ii],
                                                  w[, Jay])
    }
    names(extra$percentile) <- y.names
  }), list( .elocat = elocat, .llocat = llocat,
            .escale = escale, .lscale = lscale,
            .fittedMean = fittedMean,
            .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 2
    Mdiv2 <- extra$Mdiv2
    ymat <- matrix(y, extra$n, extra$Mdiv2)
    kappamat <- matrix(extra$kappa, extra$n, extra$Mdiv2, byrow = TRUE)

    locat <- eta2theta(eta[, 2 * (1:Mdiv2) - 1, drop = FALSE],
                       .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2 * (1:Mdiv2)    , drop = FALSE],
                       .lscale , earg = .escale )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = c(ymat), location = c(locat),
                              scale = c(Scale), kappa = c(kappamat),
                              log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .kappa = kappa ))),
  vfamily = c("alaplace2"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    extra    <- object@extra
    locat    <- eta2theta(eta[, c(TRUE, FALSE)], .llocat , .elocat )
    Scale    <- eta2theta(eta[, c(FALSE, TRUE)], .lscale , .escale )
    kappamat <- matrix(extra$kappa, extra$n, extra$Mdiv2, byrow = TRUE)
    ralap(nsim * length(Scale), location = c(locat),
          scale = c(Scale), kappa = c(kappamat))
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .kappa = kappa ))),






  deriv = eval(substitute(expression({
    M1 <- 2
    Mdiv2 <- extra$Mdiv2
    ymat <- matrix(y, n, Mdiv2)

    locat <- eta2theta(eta[, M1 * (1:(Mdiv2)) - 1, drop = FALSE],
                      .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, M1 * (1:(Mdiv2))    , drop = FALSE],
                      .lscale , earg = .escale )


    kappamat <- matrix(extra$kappa, n, Mdiv2, byrow = TRUE)
    zedd <- abs(ymat - locat) / Scale
    dl.dlocat <- sqrt(2) * ifelse(ymat >= locat, kappamat, 1/kappamat) *
                 sign(ymat - locat) / Scale
    dl.dscale <- sqrt(2) * ifelse(ymat >= locat, kappamat, 1/kappamat) *
                 zedd / Scale - 1 / Scale
    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat,
            .kappa = kappa ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, M)

    d2l.dlocat2 <- 2 / Scale^2
    d2l.dscale2 <- 1 / Scale^2

    wz[, M1*(1:Mdiv2) - 1] <- d2l.dlocat2 * dlocat.deta^2
    wz[, M1*(1:Mdiv2)    ] <- d2l.dscale2 * dscale.deta^2

    c(w) * wz
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))))
}  # End of alaplace2().











alaplace1.control <- function(maxit = 100, ...) {
    list(maxit = maxit)
}









 alaplace1 <-
  function(tau = NULL,
           llocation = "identitylink",
           ilocation = NULL,
           kappa = sqrt(tau/(1-tau)),
           Scale.arg = 1,
           ishrinkage = 0.95,
           parallel.locat = TRUE  ~ 0,  # FALSE,
           digt = 4,
           idf.mu = 3,
           zero = NULL,
           imethod = 1) {



  apply.parint.locat <- FALSE



  
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")
  if (length(tau) &&
      max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")

  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")


  llocation <- llocation

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (!is.Numeric(Scale.arg, positive = TRUE))
    stop("bad input for argument 'Scale.arg'")




  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")





  new("vglmff",
  blurb = c("One-parameter asymmetric Laplace distribution\n\n",
            "Links:      ",
            namesof("location", llocat, earg = elocat),
            "\n", "\n",
            "Mean:       location + scale * (1/kappa - kappa) / ",
                         "sqrt(2)", "\n",
            "Quantiles:  location", "\n",
            "Variance:   scale^2 * (1 + kappa^4) / (2 * kappa^2)"),




  constraints = eval(substitute(expression({

    onemat <- matrix(1, M, 1)
    constraints.orig <- constraints


    cm1.locat <- diag(M)
    cmk.locat <- onemat
    con.locat <- cm.VGAM(cmk.locat,
                         x = x, bool = .parallel.locat ,
                         constraints = constraints.orig,
                         apply.int = .apply.parint.locat ,
                         cm.default           = cm1.locat,
                         cm.intercept.default = cm1.locat)
   
    
    constraints <- con.locat

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel.locat = parallel.locat,
            .zero = zero,
            .apply.parint.locat = apply.parint.locat ))),



  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         summary.pvalues = FALSE,
         tau   = .tau ,
         multipleResponses = FALSE,
         parameters.names = c("location"),
         kappa = .kappa)
  }, list( .kappa = kappa,
           .tau   = tau ))),
  initialize = eval(substitute(expression({
    extra$M1 <- M1 <- 1


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = if (length( .kappa ) > 1) 1 else Inf,
              ncol.y.max = if (length( .kappa ) > 1) 1 else Inf,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    extra$ncoly <- ncoly <- ncol(y)
    if ((ncoly > 1) && (length( .kappa ) > 1 ||
        length( .Scale.arg ) > 1))
      stop("response must be a vector if 'kappa' or 'Scale.arg' ",
           "has a length greater than one")

    extra$kappa <- .kappa
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)


    extra$M <- M <- max(length( .Scale.arg ),
                        ncoly,
                        length( .kappa ))  # Recycle
    extra$Scale <- rep( .Scale.arg , length = M)
    extra$kappa <- rep( .kappa , length = M)
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)
    extra$n <- n




    extra$tau.names <- tau.names <-
      paste("(tau = ", round(extra$tau, digits = .digt), ")", sep = "")
    extra$Y.names <- Y.names <- if (ncoly > 1) dimnames(y)[[2]] else "y"
    if (is.null(Y.names) || any(Y.names == ""))
      extra$Y.names <- Y.names <- paste("y", 1:ncoly, sep = "")
    extra$y.names <- y.names <-
      if (ncoly > 1) paste(Y.names, tau.names, sep = "") else tau.names

    extra$individual <- FALSE

    mynames1 <- param.names("location", M)
    predictors.names <-
        c(namesof(mynames1, .llocat , earg = .elocat , tag = FALSE))


    locat.init <- matrix(0, n, M)
    if (!length(etastart)) {

      for (jay in 1:M) {
        y.use <- if (ncoly > 1) y[, jay] else y
        if ( .imethod == 1) {
          locat.init[, jay] <- weighted.mean(y.use, w[, min(jay, ncol(w))])
        } else if ( .imethod == 2) {
          locat.init[, jay] <- median(y.use)
        } else if ( .imethod == 3) {
          Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                 y = y.use, w = w, df = .idf.mu )
          locat.init[, jay] <- c(predict(Fit5, x = x[, min(ncol(x), 2)])$y)
        } else {
          use.this <- weighted.mean(y.use, w[, min(jay, ncol(w))])
          locat.init[, jay] <- (1- .ishrinkage ) * y.use + .ishrinkage * use.this
        }


        if (length( .ilocat )) {
          locat.init <- matrix( .ilocat  , n, M, byrow = TRUE)
        }

        if ( .llocat == "loge") locat.init <- abs(locat.init)
        etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
      }
    }
    }), list( .imethod = imethod,
              .idf.mu = idf.mu,
              .ishrinkage = ishrinkage, .digt = digt,
              .elocat = elocat, .Scale.arg = Scale.arg,
              .llocat = llocat, .kappa = kappa,
              .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    if ( .fittedMean ) {
      kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
      locat <- eta2theta(eta, .llocat , earg = .elocat )
      Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
      locat + Scale * (1/kappamat - kappamat)
    } else {
      locat <- eta2theta(eta, .llocat , earg = .elocat )
      if (length(locat) > extra$n)
        dimnames(locat) <- list(dimnames(eta)[[1]], extra$y.names)
      locat
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .fittedMean = fittedMean, .Scale.arg = Scale.arg,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$M1 <- M1
    misc$multipleResponses <- TRUE

    tmp34 <- c(rep( .llocat , length = M))
    names(tmp34) <- mynames1 
    misc$link <- tmp34 # Already named

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M) {
      misc$earg[[ii]] <- .elocat
    }


    misc$expected <- TRUE
    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    misc$true.mu <- .fittedMean # @fitted is not a true mu?

    extra$percentile <- numeric(M)
    locat <- as.matrix(locat)
    for (ii in 1:M) {
      y.use <- if (ncoly > 1) y[, ii] else y
      extra$percentile[ii] <-
        100 * weighted.mean(y.use <= locat[, ii], w[, min(ii, ncol(w))])
    }
    names(extra$percentile) <- y.names

    extra$Scale.arg <- .Scale.arg
    }), list( .elocat = elocat,
              .llocat = llocat,
              .Scale.arg = Scale.arg, .fittedMean = fittedMean,
              .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    ymat <- matrix(y, extra$n, extra$M)
    locat <- eta2theta(eta, .llocat , earg = .elocat )
    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    Scale    <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = c(ymat), locat = c(locat),
                              scale = c(Scale), kappa = c(kappamat),
                              log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat,
           .llocat = llocat,
           .Scale.arg = Scale.arg, .kappa = kappa ))),
  vfamily = c("alaplace1"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    extra    <- object@extra
    locat    <- eta2theta(eta, .llocat , .elocat )
    Scale    <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    ralap(nsim * length(Scale), location = c(locat),
          scale = c(Scale), kappa = c(kappamat))
  }, list( .elocat = elocat, .llocat = llocat,
           .Scale.arg = Scale.arg, .kappa = kappa ))),



  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M)
    Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)

    locat <- eta2theta(eta, .llocat , earg = .elocat )

    kappamat <- matrix(extra$kappa, n, M, byrow = TRUE)
    zedd <- abs(ymat-locat) / Scale

    dl.dlocat <- ifelse(ymat >= locat, kappamat, 1/kappamat) *
                   sqrt(2) * sign(ymat - locat) / Scale
    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )

    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .Scale.arg = Scale.arg, .elocat = elocat,
            .llocat = llocat, .kappa = kappa ))),

  weight = eval(substitute(expression({
    d2l.dlocat2 <- 2 / Scale^2
    wz <- cbind(d2l.dlocat2 * dlocat.deta^2)

    c(w) * wz
  }), list( .Scale.arg = Scale.arg,
            .elocat = elocat, .llocat = llocat ))))
}









alaplace3.control <- function(maxit = 100, ...) {
  list(maxit = maxit)
}




 alaplace3 <-
  function(llocation = "identitylink", lscale = "loge", lkappa = "loge",
           ilocation = NULL,           iscale = NULL,   ikappa = 1.0,
           imethod = 1, zero = c("scale", "kappa")) {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lkappa <- as.list(substitute(lkappa))
  ekappa <- link2list(lkappa)
  lkappa <- attr(ekappa, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  new("vglmff",
  blurb = c("Three-parameter asymmetric Laplace distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale), ", ",
            namesof("kappa",    lkappa, earg = ekappa),
            "\n", "\n",
            "Mean:     location + scale * (1/kappa - kappa) / sqrt(2)",
            "\n",
            "Variance: Scale^2 * (1 + kappa^4) / (2 * kappa^2)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale", "kappa"),
         summary.pvalues = FALSE,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <- 
      c(namesof("location", .llocat , earg = .elocat, tag = FALSE),
        namesof("scale",    .lscale , earg = .escale, tag = FALSE),
        namesof("kappa",    .lkappa , earg = .ekappa, tag = FALSE))

    if (!length(etastart)) {
      kappa.init <- if (length( .ikappa ))
                   rep( .ikappa, length.out = n) else
                   rep( 1.0, length.out = n)
      if ( .imethod == 1) {
        locat.init <- median(y)
        scale.init <- sqrt(var(y) / 2)
      } else {
        locat.init <- y
        scale.init <- sqrt(sum(c(w)*abs(y-median(y ))) / (sum(w) *2))
      }
      locat.init <- if (length( .ilocat))
                       rep( .ilocat, length.out = n) else
                       rep(locat.init, length.out = n)
      scale.init <- if (length( .iscale))
                       rep( .iscale, length.out = n) else
                       rep(scale.init, length.out = n)
      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ),
                theta2eta(scale.init, .lscale , earg = .escale ),
                theta2eta(kappa.init, .lkappa, earg = .ekappa))
    }
  }), list( .imethod = imethod,
            .elocat = elocat, .escale = escale, .ekappa = ekappa,
            .llocat = llocat, .lscale = lscale, .lkappa = lkappa,
            .ilocat = ilocat, .iscale = iscale, .ikappa = ikappa ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kappa <- eta2theta(eta[, 3], .lkappa, earg = .ekappa)
    locat + Scale * (1/kappa - kappa) / sqrt(2)
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .ekappa = ekappa, .lkappa = lkappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat ,
                      scale    = .lscale ,
                      kappa    = .lkappa )

    misc$earg <- list(location = .elocat,
                      scale    = .escale,
                      kappa    = .ekappa )

    misc$expected = TRUE
  }), list( .elocat = elocat, .llocat = llocat,
            .escale = escale, .lscale = lscale,
            .ekappa = ekappa, .lkappa = lkappa ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kappa <- eta2theta(eta[, 3], .lkappa , earg = .ekappa )  # a matrix
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = y, locat = locat,
                              scale = Scale, kappa = kappa, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .ekappa = ekappa, .lkappa = lkappa ))),
  vfamily = c("alaplace3"),
  deriv = eval(substitute(expression({
    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    kappa <- eta2theta(eta[, 3], .lkappa, earg = .ekappa)

    zedd <- abs(y - locat) / Scale
    dl.dlocat <- sqrt(2) * ifelse(y >= locat, kappa, 1/kappa) *
                   sign(y-locat) / Scale
    dl.dscale <-  sqrt(2) * ifelse(y >= locat, kappa, 1/kappa) *
                 zedd / Scale - 1 / Scale
    dl.dkappa <-  1 / kappa - 2 * kappa / (1+kappa^2) -
                 (sqrt(2) / Scale) *
                 ifelse(y > locat, 1, -1/kappa^2) * abs(y-locat)  

    dlocat.deta <- dtheta.deta(locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dkappa.deta <- dtheta.deta(kappa, .lkappa, earg = .ekappa)

    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta,
                 dl.dkappa * dkappa.deta)
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat,
            .ekappa = ekappa, .lkappa = lkappa ))),
  weight = eval(substitute(expression({
    d2l.dlocat2 <- 2 / Scale^2
    d2l.dscale2 <- 1 / Scale^2
    d2l.dkappa2 <- 1 / kappa^2 + 4 / (1+kappa^2)^2
    d2l.dkappadloc <- -sqrt(8) / ((1+kappa^2) * Scale)
    d2l.dkappadscale <- -(1-kappa^2) / ((1+kappa^2) * kappa * Scale)
    wz <- matrix(0, nrow = n, dimm(M))
    wz[,iam(1, 1, M)] <- d2l.dlocat2 * dlocat.deta^2
    wz[,iam(2, 2, M)] <- d2l.dscale2 * dscale.deta^2
    wz[,iam(3, 3, M)] <- d2l.dkappa2 * dkappa.deta^2
    wz[,iam(1, 3, M)] <- d2l.dkappadloc * dkappa.deta * dlocat.deta
    wz[,iam(2, 3, M)] <- d2l.dkappadscale  * dkappa.deta * dscale.deta
    c(w) * wz
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))))
}









dlaplace <- function(x, location = 0, scale = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  logdensity <- (-abs(x-location)/scale) - log(2*scale)
  if (log.arg) logdensity else exp(logdensity)
}



plaplace <- function(q, location = 0, scale = 1,
                     lower.tail = TRUE, log.p =FALSE) {
  zedd <- (q - location) / scale

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  L <- max(length(q), length(location), length(scale))
  if (length(q)        != L) q        <- rep(q,        length.out = L)
  if (length(location) != L) location <- rep(location, length.out = L)
  if (length(scale)    != L) scale    <- rep(scale,    length.out = L)


  if (lower.tail) {
    if (log.p) {
      ans <- ifelse(q < location, log(0.5) + zedd, log1p(- 0.5 * exp(-zedd)))
    } else {
      ans <- ifelse(q < location, 0.5 * exp(zedd), 1 - 0.5 * exp(-zedd))
    }
  } else {
    if (log.p) {
      ans <- ifelse(q < location, log1p(- 0.5 * exp(zedd)), log(0.5) - zedd)
    } else {
      ans <- ifelse(q < location, 1 - 0.5 * exp(zedd), 0.5 * exp(-zedd))
    }
  }
  ans[scale <= 0] <- NaN
  ans
}



qlaplace <- function(p, location = 0, scale = 1,
                     lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")


  L <- max(length(p), length(location), length(scale))
  if (length(p)        != L) p        <- rep(p,        length.out = L)
  if (length(location) != L) location <- rep(location, length.out = L)
  if (length(scale)    != L) scale    <- rep(scale,    length.out = L)


  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- location - sign(exp(ln.p)-0.5) * scale *
             log(2 * ifelse(exp(ln.p) < 0.5, exp(ln.p), -expm1(ln.p)))
    } else {
      ans <- location - sign(p-0.5) * scale * log(2 * ifelse(p < 0.5, p, 1-p))
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- location - sign(0.5 - exp(ln.p)) * scale *
             log(2 * ifelse(-expm1(ln.p) < 0.5, -expm1(ln.p), exp(ln.p)))
     # ans[ln.p > 0] <- NaN
    } else { 
      ans <- location - sign(0.5 - p) * scale *
             log(2 * ifelse(p > 0.5, 1 - p, p))
    }
  }

  ans[scale <= 0] <- NaN
  ans  
}



rlaplace <- function(n, location = 0, scale = 1) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (!is.Numeric(scale, positive = TRUE))
    stop("'scale' must be positive")

  location <- rep(location, length.out = use.n)
  scale    <- rep(scale,    length.out = use.n)
  rrrr     <- runif(use.n)



  location - sign(rrrr - 0.5) * scale *
  (log(2) + ifelse(rrrr < 0.5, log(rrrr), log1p(-rrrr)))
}




  
 laplace <- function(llocation = "identitylink", lscale = "loge",
                     ilocation = NULL, iscale = NULL,
                     imethod = 1,
                     zero = "scale") {

  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  new("vglmff",
  blurb = c("Two-parameter Laplace distribution\n\n",
            "Links:    ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),
            "\n", "\n",
            "Mean:     location", "\n",
            "Variance: 2*scale^2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("location", "scale"),
         summary.pvalues = FALSE,
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    predictors.names <- 
      c(namesof("location", .llocat , earg = .elocat, tag = FALSE),
        namesof("scale",    .lscale , earg = .escale, tag = FALSE))


    if (!length(etastart)) {
      if ( .imethod == 1) {
        locat.init <- median(y)
        scale.init <- sqrt(var(y) / 2)
      } else if ( .imethod == 2) {
        locat.init <- weighted.mean(y, w)
        scale.init <- sqrt(var(y) / 2)
      } else {
        locat.init <- median(y)
        scale.init <- sqrt(sum(c(w)*abs(y-median(y ))) / (sum(w) *2))
      }
      locat.init <- if (length( .ilocat))
                       rep( .ilocat, length.out = n) else
                       rep(locat.init, length.out = n)
      scale.init <- if (length( .iscale))
                       rep( .iscale, length.out = n) else
                       rep(scale.init, length.out = n)
      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ),
                theta2eta(scale.init, .lscale , earg = .escale ))
    }
  }), list( .imethod = imethod,
            .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale,
            .ilocat = ilocat, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat )
  }, list( .elocat = elocat, .llocat = llocat ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale )

    misc$earg <- list(location = .elocat , scale = .escale )

    misc$expected <- TRUE
    misc$RegCondOK <- FALSE # Save this for later
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dlaplace(x = y, locat = locat,
                                 scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .escale = escale, .lscale = lscale,
           .elocat = elocat, .llocat = llocat ))),
  vfamily = c("laplace"),
  deriv = eval(substitute(expression({
    Locat <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    Scale <- eta2theta(eta[, 2], .lscale , earg = .escale )

    zedd <- abs(y-Locat) / Scale
    dl.dLocat <- sign(y - Locat) / Scale
    dl.dscale <-  zedd / Scale - 1 / Scale

    dLocat.deta <- dtheta.deta(Locat, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    c(w) * cbind(dl.dLocat * dLocat.deta,
                 dl.dscale    * dscale.deta)
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))),
  weight = eval(substitute(expression({
    d2l.dLocat2 <- d2l.dscale2 <- 1 / Scale^2
    wz <- matrix(0, nrow = n, ncol = M)  # diagonal
    wz[,iam(1, 1, M)] <- d2l.dLocat2 * dLocat.deta^2
    wz[,iam(2, 2, M)] <- d2l.dscale2 * dscale.deta^2
    c(w) * wz
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat ))))
}



fff.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 fff <- function(link = "loge",
                 idf1 = NULL, idf2 = NULL, nsimEIM = 100,  # ncp = 0,
                 imethod = 1, zero = NULL) {
  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10)
    stop("argument 'nsimEIM' should be an integer greater than 10")

  ncp <- 0
  if (any(ncp != 0))
    warning("not sure about ncp != 0 wrt dl/dtheta")



  new("vglmff",
  blurb = c("F-distribution\n\n",
            "Links:    ",
            namesof("df1", link, earg = earg), ", ",
            namesof("df2", link, earg = earg),
            "\n", "\n",
            "Mean:     df2/(df2-2) provided df2>2 and ncp = 0", "\n",
            "Variance: ",
            "2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4)) ",
            "provided df2>4 and ncp = 0"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("df1", "df2"),
         zero = .zero )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)



    predictors.names <-
      c(namesof("df1", .link , earg = .earg , tag = FALSE),
        namesof("df2", .link , earg = .earg , tag = FALSE))


    if (!length(etastart)) {
      if ( .imethod == 1) {
        df2.init <- b <- 2*mean(y) / (mean(y)-1)
        df1.init <- 2*b^2*(b-2)/(var(y)*(b-2)^2 * (b-4) - 2*b^2)
        if (df2.init < 4) df2.init <- 5
        if (df1.init < 2) df1.init <- 3
      } else {
            df2.init <- b <- 2*median(y) / (median(y)-1)
            summy <- summary(y)
            var.est <- summy[5] - summy[2]
            df1.init <- 2*b^2*(b-2)/(var.est*(b-2)^2 * (b-4) - 2*b^2)
        }
        df1.init <- if (length( .idf1 ))
                       rep( .idf1 , length.out = n) else
                       rep(df1.init, length.out = n)
        df2.init <- if (length( .idf2 ))
                       rep( .idf2 , length.out = n) else
                       rep(1, length.out = n)
        etastart <- cbind(theta2eta(df1.init, .link , earg = .earg ),
                          theta2eta(df2.init, .link , earg = .earg ))
    }
  }), list( .imethod = imethod, .idf1 = idf1, .earg = earg,
           .idf2 = idf2, .link = link ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    ans <- df2 * NA
    ans[df2>2] <- df2[df2>2] / (df2[df2>2]-2)
    ans
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(df1 = .link , df2 = .link )

    misc$earg <- list(df1 = .earg , df2 = .earg )

    misc$nsimEIM <- .nsimEIM
    misc$ncp <- .ncp
  }), list( .link = link, .earg = earg,
            .ncp = ncp,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    df1 <- eta2theta(eta[, 1], .link , earg = .earg )
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * df(x = y, df1 = df1, df2 = df2,
                           ncp = .ncp , log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg, .ncp=ncp ))),
  vfamily = c("fff"),
  deriv = eval(substitute(expression({
    df1 <- eta2theta(eta[, 1], .link , earg = .earg )
    df2 <- eta2theta(eta[, 2], .link , earg = .earg )
    dl.ddf1 <- 0.5*digamma(0.5*(df1+df2)) + 0.5 + 0.5*log(df1/df2) +
              0.5*log(y) - 0.5*digamma(0.5*df1) -
              0.5*(df1+df2)*(y/df2) / (1 + df1*y/df2) -
              0.5*log1p(df1*y/df2)
    dl.ddf2 <- 0.5*digamma(0.5*(df1+df2)) - 0.5*df1/df2 - 
              0.5*digamma(0.5*df2) -
              0.5*(df1+df2) * (-df1*y/df2^2) / (1 + df1*y/df2) -
              0.5*log1p(df1*y/df2)
    ddf1.deta <- dtheta.deta(df1, .link , earg = .earg )
    ddf2.deta <- dtheta.deta(df2, .link , earg = .earg )
    dthetas.detas <- cbind(ddf1.deta, ddf2.deta)
      c(w) * dthetas.detas * cbind(dl.ddf1, dl.ddf2)
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for (ii in 1:( .nsimEIM )) {
      ysim <- rf(n = n, df1=df1, df2=df2)
      dl.ddf1 <- 0.5*digamma(0.5*(df1+df2)) + 0.5 + 0.5*log(df1/df2) +
                0.5*log(ysim) - 0.5*digamma(0.5*df1) -
                0.5*(df1+df2)*(ysim/df2) / (1 + df1*ysim/df2) -
                0.5*log1p(df1*ysim/df2)
      dl.ddf2 <- 0.5*digamma(0.5*(df1+df2)) - 0.5*df1/df2 - 
                0.5*digamma(0.5*df2) -
                0.5*(df1+df2) * (-df1*ysim/df2^2)/(1 + df1*ysim/df2) -
                0.5*log1p(df1*ysim/df2)
      rm(ysim)
      temp3 <- cbind(dl.ddf1, dl.ddf2)
      run.varcov <- ((ii-1) * run.varcov +
                 temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
    }
    wz <- if (intercept.only)
        matrix(colMeans(run.varcov),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov

    wz <- c(w) * wz * dthetas.detas[, ind1$row] *
                     dthetas.detas[, ind1$col]
    wz
  }), list( .link = link, .earg = earg, .nsimEIM = nsimEIM,
            .ncp = ncp ))))
}




 hyperg <- function(N = NULL, D = NULL,
                    lprob = "logit",
                    iprob = NULL) {

  inputN <- is.Numeric(N, positive = TRUE)
  inputD <- is.Numeric(D, positive = TRUE)
  if (inputD && inputN)
    stop("only one of 'N' and 'D' is to be inputted")
  if (!inputD && !inputN)
    stop("one of 'N' and 'D' needs to be inputted")


  lprob <- as.list(substitute(lprob))
  earg <- link2list(lprob)
  lprob <- attr(earg, "function.name")



  new("vglmff",
  blurb = c("Hypergeometric distribution\n\n",
            "Link:     ",
            namesof("prob", lprob, earg = earg), "\n",
            "Mean:     D/N\n"),
  initialize = eval(substitute(expression({
    NCOL <- function (x)
        if (is.array(x) && length(dim(x)) > 1 ||
        is.data.frame(x)) ncol(x) else as.integer(1)
    if (NCOL(y) == 1) {
        if (is.factor(y)) y <- y != levels(y)[1]
        nn <- rep(1, length.out = n)
        if (!all(y >= 0 & y <= 1))
            stop("response values must be in [0, 1]")
        mustart <- (0.5 + w * y) / (1 + w)
        no.successes <- w * y
        if (any(abs(no.successes - round(no.successes)) > 0.001))
            stop("Number of successes must be integer-valued")
    } else if (NCOL(y) == 2) {
        if (any(abs(y - round(y)) > 0.001))
            stop("Count data must be integer-valued")
        nn <- y[, 1] + y[, 2]
        y <- ifelse(nn > 0, y[, 1]/nn, 0)
        w <- w * nn
        mustart <- (0.5 + nn * y) / (1 + nn)
        mustart[mustart >= 1] <- 0.95
    } else
         stop("Response not of the right form")

    predictors.names <-
      namesof("prob", .lprob , earg = .earg , tag = FALSE)
    extra$Nvector <- .N
    extra$Dvector <- .D
    extra$Nunknown <- length(extra$Nvector) == 0
    if (!length(etastart)) {
        init.prob <- if (length( .iprob))
                      rep( .iprob, length.out = n) else
                      mustart
            etastart <- matrix(init.prob, n, ncol(cbind(y )))

    }
  }), list( .lprob = lprob, .earg = earg, .N = N, .D = D,
            .iprob = iprob ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, .lprob, earg = .earg )
  }, list( .lprob = lprob, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c("prob" = .lprob) 

    misc$earg <- list("prob" = .earg ) 

    misc$Dvector <- .D
    misc$Nvector <- .N
  }), list( .N = N, .D = D, .lprob = lprob, .earg = earg ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, .lprob, earg = .earg )
  }, list( .lprob = lprob, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    N <- extra$Nvector
    Dvec <- extra$Dvector
    prob <- mu
    yvec <- w * y
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
      if (extra$Nunknown) {
        tmp12 <- Dvec * (1-prob) / prob


        (lgamma(1+tmp12) + lgamma(1+Dvec/prob-w) -
         lgamma(1+tmp12-w+yvec) - lgamma(1+Dvec/prob))
      } else {


        (lgamma(1+N*prob) + lgamma(1+N*(1-prob)) -
         lgamma(1+N*prob-yvec) -
         lgamma(1+N*(1-prob) -w + yvec))
      }
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .earg = earg ))), 
  vfamily = c("hyperg"),
  deriv = eval(substitute(expression({
    prob <- mu   # equivalently, eta2theta(eta, .lprob, earg = .earg )
    dprob.deta <- dtheta.deta(prob, .lprob, earg = .earg )
    Dvec <- extra$Dvector
    Nvec <- extra$Nvector
    yvec <- w * y
    if (extra$Nunknown) {
      tmp72 <- -Dvec / prob^2
      tmp12 <-  Dvec * (1-prob) / prob
      dl.dprob <- tmp72 * (digamma(1 + tmp12) +
                 digamma(1 + Dvec/prob -w) -
           digamma(1 + tmp12-w+yvec) - digamma(1 + Dvec/prob))
    } else {
      dl.dprob <- Nvec * (digamma(1+Nvec*prob) -
                 digamma(1+Nvec*(1-prob)) -
                 digamma(1+Nvec*prob-yvec) +
                 digamma(1+Nvec*(1-prob)-w+yvec))
    }
    c(w) * dl.dprob * dprob.deta
  }), list( .lprob = lprob, .earg = earg ))),
  weight = eval(substitute(expression({
    if (extra$Nunknown) {
      tmp722 <- tmp72^2
      tmp13 <- 2*Dvec / prob^3
      d2l.dprob2 <- tmp722 * (trigamma(1 + tmp12) + 
                   trigamma(1 + Dvec/prob - w) -
                   trigamma(1 + tmp12 - w + yvec) -
                   trigamma(1 + Dvec/prob)) +
                   tmp13 * (digamma(1 + tmp12) +
                   digamma(1 + Dvec/prob - w) -
                   digamma(1 + tmp12 - w + yvec) -
                   digamma(1 + Dvec/prob))
    } else {
      d2l.dprob2 <- Nvec^2 * (trigamma(1+Nvec*prob) +
                   trigamma(1+Nvec*(1-prob)) -
                   trigamma(1+Nvec*prob-yvec) -
                   trigamma(1+Nvec*(1-prob)-w+yvec))
    }
    d2prob.deta2 <- d2theta.deta2(prob, .lprob , earg = .earg )

    wz <- -(dprob.deta^2) * d2l.dprob2
    wz <- c(w) * wz
    wz[wz < .Machine$double.eps] <- .Machine$double.eps
    wz
    }), list( .lprob = lprob, .earg = earg ))))
}





dbenini <- function(x, y0, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(shape), length(y0))
  if (length(x)        != N) x        <- rep(x,        length.out = N)
  if (length(shape)    != N) shape    <- rep(shape,    length.out = N)
  if (length(y0)       != N) y0       <- rep(y0,       length.out = N)

  logdensity <- rep(log(0), length.out = N)
  xok <- (x > y0)
  tempxok <- log(x[xok]/y0[xok])
  logdensity[xok] <- log(2*shape[xok]) - shape[xok] * tempxok^2 +
                    log(tempxok) - log(x[xok])
  logdensity[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) logdensity else exp(logdensity)
}



pbenini <- function(q, y0, shape, lower.tail = TRUE, log.p = FALSE) {
  if (!is.Numeric(q))
    stop("bad input for argument 'q'")
  if (!is.Numeric(shape, positive = TRUE))
    stop("bad input for argument 'shape'")
  if (!is.Numeric(y0, positive = TRUE))
    stop("bad input for argument 'y0'")
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  N <- max(length(q), length(shape), length(y0))
  if (length(q)        != N) q        <- rep(q,        length.out = N)
  if (length(shape)    != N) shape    <- rep(shape,    length.out = N)
  if (length(y0)       != N) y0       <- rep(y0,       length.out = N)

  ans <- y0 * 0
  ok <- q > y0


  if (lower.tail) {
    if (log.p) {
      ans[ok] <- log(-expm1(-shape[ok] * (log(q[ok]/y0[ok]))^2))
      ans[q <= y0 ] <- -Inf
    } else {
      ans[ok] <- -expm1(-shape[ok] * (log(q[ok]/y0[ok]))^2)
    }
  } else {
    if (log.p) {
      ans[ok] <- -shape[ok] * (log(q[ok]/y0[ok]))^2
      ans[q <= y0] <- 0
    } else {
      ans[ok] <- exp(-shape[ok] * (log(q[ok]/y0[ok]))^2)
      ans[q <= y0] <- 1
    }
  } 

  ans
}



qbenini <- function(p, y0, shape, lower.tail = TRUE, log.p = FALSE) {


  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- y0 * exp(sqrt(-log(-expm1(ln.p)) / shape))
    } else {
      ans <- y0 * exp(sqrt(-log1p(-p) / shape))
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- y0 * exp(sqrt(-ln.p / shape))
    } else { 
      ans <-  y0 * exp(sqrt(-log(p) / shape))
    }
  }
  ans[y0 <= 0] <- NaN
  ans
}



rbenini <- function(n, y0, shape) {
  y0 * exp(sqrt(-log(runif(n)) / shape))
}





 benini1 <- function(y0 = stop("argument 'y0' must be specified"),
                     lshape = "loge",
                     ishape = NULL, imethod = 1, zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")
  if (!is.Numeric(y0, positive = TRUE))
   stop("bad input for argument 'y0'")




  new("vglmff",
  blurb = c("1-parameter Benini distribution\n\n",
            "Link:    ",
            namesof("shape", lshape, earg = eshape),
            "\n", "\n",
            "Median:     qbenini(p = 0.5, y0, shape)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("shape"),
         lshape = .lshape ,
         eshape = .eshape )
  }, list( .eshape = eshape,
           .lshape = lshape))),

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
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    mynames1 <- paste("shape", if (ncoly > 1) 1:ncoly else "", sep = "")
    predictors.names <-
      namesof(mynames1, .lshape , earg = .eshape , tag = FALSE)

    extra$y0 <- matrix( .y0 , n, ncoly, byrow = TRUE)
    if (any(y <= extra$y0))
      stop("some values of the response are > argument 'y0' values")


    if (!length(etastart)) {
      probs.y <- (1:3) / 4
      qofy <- quantile(rep(y, times = w), probs = probs.y)
      if ( .imethod == 1) {
        shape.init <- mean(-log1p(-probs.y) / (log(qofy))^2)
      } else {
        shape.init <- median(-log1p(-probs.y) / (log(qofy))^2)
      }
    shape.init <- matrix(if (length( .ishape )) .ishape else shape.init,
                        n, ncoly, byrow = TRUE)
    etastart <- cbind(theta2eta(shape.init, .lshape , earg = .eshape ))
  }
  }), list( .imethod = imethod,
            .ishape = ishape,
            .lshape = lshape, .eshape = eshape,
            .y0 = y0 ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )


    qbenini(p = 0.5, y0 = extra$y0, shape)
  }, list( .lshape = lshape, .eshape = eshape ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <- c(rep( .lshape , length = ncoly))
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:ncoly) {
      misc$earg[[ii]] <- .eshape
    }

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE


    extra$y0 <- .y0

  }), list( .lshape = lshape,
            .eshape = eshape, .y0 = y0 ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    y0 <- extra$y0
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dbenini(x = y, y0 = y0, shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape, .eshape = eshape ))),
  vfamily = c("benini1"),





  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    extra <- object@extra
    shape <- eta2theta(eta, .lshape , earg = .eshape )
    y0 <- extra$y0
    rbenini(nsim * length(shape), y0 = y0, shape = shape)
  }, list( .lshape = lshape, .eshape = eshape ))),





  deriv = eval(substitute(expression({
    shape <- eta2theta(eta, .lshape , earg = .eshape )

    y0 <- extra$y0
    dl.dshape <- 1/shape - (log(y/y0))^2

    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )

    c(w) * dl.dshape * dshape.deta
  }), list( .lshape = lshape, .eshape = eshape ))),
  weight = eval(substitute(expression({
    ned2l.dshape2 <- 1 / shape^2
    wz <- ned2l.dshape2 * dshape.deta^2
    c(w) * wz
  }), list( .lshape = lshape, .eshape = eshape ))))
}






 dpolono  <- function (x, meanlog = 0, sdlog = 1, bigx = 170, ...) {
  mapply(function(x, meanlog, sdlog, ...) {
    if (abs(x) > floor(x)) { # zero prob for -ve or non-integer
      0
    } else
    if (x == Inf) { # 20141215 KaiH  
      0
    } else
    if (x > bigx) {
      z <- (log(x) - meanlog) / sdlog
      (1 + (z^2 + log(x) - meanlog - 1) / (2 * x * sdlog^2)) *
      exp(-0.5 * z^2) / (sqrt(2 * pi) * sdlog * x)
    } else
       integrate( function(t) exp(t * x - exp(t) -
                              0.5 * ((t - meanlog) / sdlog)^2),
                 lower = -Inf, upper = Inf, ...)$value / (sqrt(2 * pi) *
                 sdlog * exp(lgamma(x + 1.0)))
  }, x, meanlog, sdlog, ...)
}




ppolono <- function(q, meanlog = 0, sdlog = 1,
                    isOne = 1 - sqrt( .Machine$double.eps ), ...) {


 .cumprob <- rep(0, length(q))
 .cumprob[q == Inf] <- 1  # special case


 q <- floor(q)
 ii <-  -1
 while (any(xActive <- ((.cumprob < isOne) & (q > ii))))
    .cumprob[xActive] <- .cumprob[xActive] +
      dpolono(ii <- (ii+1), meanlog, sdlog, ...)
 .cumprob
}









rpolono <- function(n, meanlog = 0, sdlog = 1) {
  lambda <- rlnorm(n = n, meanlog = meanlog, sdlog = sdlog)
  rpois(n = n, lambda = lambda)
}














dtriangle <- function(x, theta, lower = 0, upper = 1, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(theta), length(lower), length(upper))
  if (length(x)     != N) x     <- rep(x,     length.out = N)
  if (length(theta) != N) theta <- rep(theta, length.out = N)
  if (length(lower) != N) lower <- rep(lower, length.out = N)
  if (length(upper) != N) upper <- rep(upper, length.out = N)

  denom1 <- ((upper-lower)*(theta-lower))
  denom2 <- ((upper-lower)*(upper-theta))
  logdensity <- rep(log(0), length.out = N)
  xok.neg <- (lower <  x) & (x <= theta)
  xok.pos <- (theta <= x) & (x <  upper)
  logdensity[xok.neg] =
    log(2 * (x[xok.neg] - lower[xok.neg]) / denom1[xok.neg])
  logdensity[xok.pos] =
    log(2 * (upper[xok.pos] - x[xok.pos]) / denom2[xok.pos])
  logdensity[lower >= upper] <- NaN
  logdensity[lower >  theta] <- NaN
  logdensity[upper <  theta] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}


rtriangle <- function(n, theta, lower = 0, upper = 1) {


  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  if (!is.Numeric(theta))
    stop("bad input for argument 'theta'")
  if (!is.Numeric(lower))
    stop("bad input for argument 'lower'")
  if (!is.Numeric(upper))
    stop("bad input for argument 'upper'")
  if (!all(lower < theta & theta < upper))
    stop("lower < theta < upper values are required")

  N <- use.n
  lower <- rep(lower, length.out = N)
  upper <- rep(upper, length.out = N)
  theta <- rep(theta, length.out = N)
  t1 <- sqrt(runif(n))
  t2 <- sqrt(runif(n))
  ifelse(runif(n) < (theta - lower) / (upper - lower),
         lower + (theta - lower) * t1,
         upper - (upper - theta) * t2)
}



qtriangle <- function(p, theta, lower = 0, upper = 1,
                      lower.tail = TRUE, log.p = FALSE) {

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  N <- max(length(p), length(theta), length(lower), length(upper))
  if (length(p)     != N) p     <- rep(p,     length.out = N)
  if (length(theta) != N) theta <- rep(theta, length.out = N)
  if (length(lower) != N) lower <- rep(lower, length.out = N)
  if (length(upper) != N) upper <- rep(upper, length.out = N)

  ans <- NA_real_ * p
  if (lower.tail) {
    if (log.p) {
      Neg <- (exp(ln.p) <= (theta - lower) / (upper - lower))
      temp1 <- exp(ln.p) * (upper - lower) * (theta - lower)
      Pos <- (exp(ln.p) >= (theta - lower) / (upper - lower))
      pstar <- (exp(ln.p) - (theta - lower) / (upper - lower)) /
               ((upper - theta) / (upper - lower))
    } else {
      Neg <- (p <= (theta - lower) / (upper - lower))
      temp1 <- p * (upper - lower) * (theta - lower)
      Pos <- (p >= (theta - lower) / (upper - lower))
      pstar <- (p - (theta - lower) / (upper - lower)) /
               ((upper - theta) / (upper - lower))
    }
  } else {
    if (log.p) {
      ln.p <- p
      Neg <- (exp(ln.p) >= (upper- theta) / (upper - lower))
      temp1 <- -expm1(ln.p) * (upper - lower) * (theta - lower)
      Pos <- (exp(ln.p) <= (upper- theta) / (upper - lower))
      pstar <- (-expm1(ln.p) - (theta - lower) / (upper - lower)) /
               ((upper - theta) / (upper - lower))
    } else { 
      Neg <- (p >= (upper- theta) / (upper - lower))
      temp1 <- (1 - p) * (upper - lower) * (theta - lower)
      Pos <- (p <= (upper- theta) / (upper - lower))
      pstar <- ((upper- theta) / (upper - lower) - p) /
               ((upper - theta) / (upper - lower))
    }
  }
  ans[ Neg] <- lower[ Neg] + sqrt(temp1[ Neg])
  if (any(Pos)) {
    qstar <- cbind(1 - sqrt(1-pstar), 1 + sqrt(1-pstar))
    qstar <- qstar[Pos,, drop = FALSE]
    qstar <- ifelse(qstar[, 1] >= 0 & qstar[, 1] <= 1,
                    qstar[, 1],
                    qstar[, 2])
    ans[Pos] <- theta[Pos] + qstar * (upper - theta)[Pos]
  }
  
  ans[theta < lower | theta > upper] <- NaN
  ans
}



ptriangle <- function(q, theta, lower = 0, upper = 1, 
                      lower.tail = TRUE, log.p = FALSE) {

  N <- max(length(q), length(theta), length(lower), length(upper))
  if (length(q)     != N) q     <- rep(q,     length.out = N)
  if (length(theta) != N) theta <- rep(theta, length.out = N)
  if (length(lower) != N) lower <- rep(lower, length.out = N)
  if (length(upper) != N) upper <- rep(upper, length.out = N)

  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")

  ans <- q * 0
  qstar <- (q - lower)^2 / ((upper - lower) * (theta - lower))
  Neg <- (lower <= q & q <= theta)


  ans[Neg] <- if (lower.tail) {
    if (log.p) {
      (log(qstar))[Neg]
    } else {
      qstar[Neg]
    }
  } else {
    if (log.p) {
      (log1p(-qstar))[Neg]
    } else {
      1 - qstar[Neg]
    }
  } 

  Pos <- (theta <= q & q <= upper)
  qstar <- (q - theta) / (upper-theta)
  
  if (lower.tail) {
    if (log.p) {
      ans[Pos] <- log(((theta-lower)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) * (upper-theta) / (upper - lower))[Pos])
      ans[q <= lower] <- -Inf
      ans[q >= upper] <- 0
    } else {
      ans[Pos] <- ((theta-lower)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) * (upper-theta) / (upper - lower))[Pos]
      ans[q <= lower] <- 0
      ans[q >= upper] <- 1
    }
  } else {
    if (log.p) {
      ans[Pos] <- log(((upper - theta)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) * (upper-theta) / (upper - lower))[Pos])
      ans[q <= lower] <- 0
      ans[q >= upper] <- -Inf
    } else {
      ans[Pos] <- ((upper - theta)/(upper-lower))[Pos] +
                  (qstar * (2-qstar) * (upper-theta) / (upper - lower))[Pos]
      ans[q <= lower] <- 1
      ans[q >= upper] <- 0
    }
  } 

  ans[theta < lower | theta > upper] <- NaN
  ans
}







triangle.control <- function(stepsize = 0.33, maxit = 100, ...) {
  list(stepsize = stepsize, maxit = maxit)
}


 triangle <-
  function(lower = 0, upper = 1,
           link = extlogit(min = 0, max = 1),
           itheta = NULL) {






  if (!is.Numeric(lower))
    stop("bad input for argument 'lower'")
  if (!is.Numeric(upper))
    stop("bad input for argument 'upper'")
  if (!all(lower < upper))
    stop("lower < upper values are required")

  if (length(itheta) && !is.Numeric(itheta))
    stop("bad input for 'itheta'")




  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(earg$min) && any(earg$min != lower))
    stop("argument 'lower' does not match the 'link'")
  if (length(earg$max) && any(earg$max != upper))
    stop("argument 'upper' does not match the 'link'")



  new("vglmff",
  blurb = c("Triangle distribution\n\n",
            "Link:    ",
            namesof("theta", link, earg = earg)),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("theta"),
         link = .link )
  }, list( .link = link ))),

  initialize = eval(substitute(expression({

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)




    extra$lower <- rep( .lower , length.out = n)
    extra$upper <- rep( .upper , length.out = n)

    if (any(y <= extra$lower | y >= extra$upper))
      stop("some y values in [lower,upper] detected")

    predictors.names <-
      namesof("theta", .link , earg = .earg , tag = FALSE)


    if (!length(etastart)) {
      Theta.init <- if (length( .itheta )) .itheta else {
        weighted.mean(y, w)
      }
      Theta.init <- rep(Theta.init, length = n)
      etastart <- theta2eta(Theta.init, .link , earg = .earg )
    }
  }), list( .link = link, .earg = earg, .itheta=itheta,
            .upper = upper, .lower = lower ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Theta <- eta2theta(eta, .link , earg = .earg )
    lower <- extra$lower
    upper <- extra$upper

    mu1<-  (lower + upper + Theta) / 3

    mu1
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(theta = .link )

    misc$earg <- list(theta = .earg )

    misc$expected <- TRUE
  }), list( .link = link, .earg = earg ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    Theta <- eta2theta(eta, .link , earg = .earg )
    lower <- extra$lower
    upper <- extra$upper
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtriangle(x = y, theta = Theta, lower = lower,
                                  upper = upper, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("triangle"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    extra <- object@extra
    Theta <- eta2theta(eta, .link , earg = .earg )
    lower <- extra$lower
    upper <- extra$upper
    rtriangle(nsim * length(Theta),
              theta = Theta, lower = lower, upper = upper)
  }, list( .link = link, .earg = earg ))),




  deriv = eval(substitute(expression({
    Theta       <- eta2theta(eta,     .link , earg = .earg )
    dTheta.deta <- dtheta.deta(Theta, .link , earg = .earg )

    pos <- y > Theta
    neg <- y < Theta
    lower <- extra$lower
    upper <- extra$upper

    dl.dTheta <-  0 * y
    dl.dTheta[neg] <-  -1 / (Theta[neg]-lower[neg])
    dl.dTheta[pos] <-   1 / (upper[pos]-Theta[pos])

    c(w) * dl.dTheta * dTheta.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    var.dl.dTheta <-  1 / ((Theta - lower) * (upper - Theta))
    wz <- var.dl.dTheta * dTheta.deta^2
    c(w) * wz
  }), list( .link = link, .earg = earg ))))
}







adjust0.loglaplace1 <- function(ymat, y, w, rep0) {
  rangey0 <- range(y[y > 0])
  ymat[ymat <= 0] <- min(rangey0[1] / 2, rep0)
  ymat
}


loglaplace1.control <- function(maxit = 300, ...) {
  list(maxit = maxit)
}


 loglaplace1 <- function(tau = NULL,
                     llocation = "loge",
                     ilocation = NULL,
                     kappa = sqrt(tau/(1-tau)),
                     Scale.arg = 1,
                     ishrinkage = 0.95,
                     parallel.locat = FALSE, digt = 4,
                     idf.mu = 3,
                     rep0 = 0.5,  # 0.0001,
                     minquantile = 0, maxquantile = Inf,
                     imethod = 1, zero = NULL) {

  if (length(minquantile) != 1)
    stop("bad input for argument 'minquantile'")
  if (length(maxquantile) != 1)
    stop("bad input for argument 'maxquantile'")


  if (!is.Numeric(rep0, positive = TRUE, length.arg = 1) ||
      rep0 > 1)
    stop("bad input for argument 'rep0'")
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")

  if (length(tau) && max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
      stop("arguments 'kappa' and 'tau' do not match")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation


  llocat.identity <- as.list(substitute("identitylink"))
  elocat.identity <- link2list(llocat.identity)
  llocat.identity <- attr(elocat.identity, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  if (!is.Numeric(Scale.arg, positive = TRUE))
    stop("bad input for argument 'Scale.arg'")
  if (!is.logical(parallel.locat) ||
      length(parallel.locat) != 1)
    stop("bad input for argument 'parallel.locat'")

  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")


  mystring0 <- namesof("location", llocat, earg = elocat)
  mychars <- substring(mystring0, first = 1:nchar(mystring0),
                      last = 1:nchar(mystring0))
  mychars[nchar(mystring0)] <- ", inverse = TRUE)"
  mystring1 <- paste(mychars, collapse = "")




  new("vglmff",
  blurb = c("One-parameter ",
            if (llocat == "loge") "log-Laplace" else
              c(llocat, "-Laplace"),
            " distribution\n\n",
            "Links:      ", mystring0, "\n", "\n",
          "Quantiles:  ", mystring1),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel.locat ,
                           constraints = constraints, apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel.locat = parallel.locat,
            .Scale.arg = Scale.arg, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         parameters.names = c("location"),
         llocation = .llocat )
  }, list( .llocat = llocat,
           .zero   = zero ))),

  initialize = eval(substitute(expression({
    extra$M <- M <- max(length( .Scale.arg ), length( .kappa ))  # Recycle
    extra$Scale <- rep( .Scale.arg , length = M)
    extra$kappa <- rep( .kappa     , length = M)
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)


    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




        extra$n <- n
        extra$y.names <- y.names <-
          paste("tau = ", round(extra$tau, digits = .digt), sep = "")
        extra$individual <- FALSE


        predictors.names <-
          namesof(paste("quantile(", y.names, ")", sep = ""),
                  .llocat , earg = .elocat , tag = FALSE)


        if (FALSE) {
        if (min(y) < 0)
          stop("negative response values detected")
        if ((prop.0. <- weighted.mean(1*(y == 0), w)) >= min(extra$tau))
          stop("sample proportion of 0s == ", round(prop.0., digits = 4),
               " > minimum 'tau' value. Choose larger values for 'tau'.")
        if ( .rep0 == 0.5 &&
            (ave.tau <- (weighted.mean(1*(y <= 0), w) +
             weighted.mean(1*(y <= 1), w))/2) >= min(extra$tau))
            warning("the minimum 'tau' value should be greater than ",
                 round(ave.tau, digits = 4))
        }

        if (!length(etastart)) {
            if ( .imethod == 1) {
              locat.init <- quantile(rep(y, w), probs= extra$tau) + 1/16
            } else if ( .imethod == 2) {
              locat.init <- weighted.mean(y, w)
            } else if ( .imethod == 3) {
              locat.init <- median(y)
            } else if ( .imethod == 4) {
              Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                     y = y, w = w,
                                     df = .idf.mu )
              locat.init <- c(predict(Fit5, x = x[, min(ncol(x), 2)])$y)
            } else {
              use.this <- weighted.mean(y, w)
              locat.init <- (1- .ishrinkage )*y + .ishrinkage * use.this
            }
            locat.init <- if (length( .ilocat))
                             rep( .ilocat, length.out = M) else
                             rep(locat.init, length.out = M)
            locat.init <- matrix(locat.init, n, M, byrow = TRUE)
            if ( .llocat == "loge")
                locat.init <- abs(locat.init)
            etastart <-
                cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
        }
    }), list( .imethod = imethod,
              .idf.mu = idf.mu, .rep0 = rep0,
              .ishrinkage = ishrinkage, .digt = digt,
              .elocat = elocat, .Scale.arg = Scale.arg,
              .llocat = llocat, .kappa = kappa,
              .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat.y = eta2theta(eta, .llocat , earg = .elocat )
    if ( .fittedMean ) {
      stop("Yet to do: handle 'fittedMean = TRUE'")
      kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
      Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
      locat.y + Scale * (1/kappamat - kappamat)
    } else {
      if (length(locat.y) > extra$n)
        dimnames(locat.y) <- list(dimnames(eta)[[1]], extra$y.names)
      locat.y
    }
        locat.y[locat.y < .minquantile] = .minquantile
        locat.y[locat.y > .maxquantile] = .maxquantile
        locat.y
  }, list( .elocat = elocat, .llocat = llocat,
           .minquantile = minquantile, .maxquantile = maxquantile,
           .fittedMean = fittedMean, .Scale.arg = Scale.arg,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat)

    misc$earg <- list(location = .elocat )

    misc$expected <- TRUE

    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    extra$Scale.arg <- .Scale.arg

    misc$true.mu <- .fittedMean # @fitted is not a true mu?
    misc$rep0 <- .rep0
    misc$minquantile <- .minquantile
    misc$maxquantile <- .maxquantile

    extra$percentile <- numeric(length(misc$kappa))
    locat.y <- as.matrix(locat.y)
    for (ii in 1:length(misc$kappa))
      extra$percentile[ii] <- 100 * weighted.mean(y <= locat.y[, ii], w)
  }), list( .elocat = elocat, .llocat = llocat,
            .Scale.arg = Scale.arg, .fittedMean = fittedMean,
            .minquantile = minquantile, .maxquantile = maxquantile,
            .rep0 = rep0, .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    Scale.w <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    ymat <- matrix(y, extra$n, extra$M)

    if ( .llocat == "loge")
      ymat <- adjust0.loglaplace1(ymat = ymat, y = y, w = w, rep0 = .rep0)

        w.mat <- theta2eta(ymat, .llocat , earg = .elocat )  # e.g., logoff()



    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dalap(x = c(w.mat), locat = c(eta),
                              scale = c(Scale.w), kappa = c(kappamat),
                              log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .rep0 = rep0,
           .Scale.arg = Scale.arg, .kappa = kappa ))),
  vfamily = c("loglaplace1"),
  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M)
    Scale.w <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    locat.w <- eta
    locat.y <- eta2theta(locat.w, .llocat , earg = .elocat )
    kappamat <- matrix(extra$kappa, n, M, byrow = TRUE)

    ymat <- adjust0.loglaplace1(ymat = ymat, y = y, w = w, rep0= .rep0)
    w.mat <- theta2eta(ymat, .llocat , earg = .elocat )  # e.g., logit()
    zedd <- abs(w.mat-locat.w) / Scale.w
    dl.dlocat <- ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                   sqrt(2) * sign(w.mat-locat.w) / Scale.w


    dlocat.deta <- dtheta.deta(locat.w,
                              .llocat.identity ,
                              earg = .elocat.identity )
    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .Scale.arg = Scale.arg, .rep0 = rep0,
            .llocat = llocat, .elocat = elocat,
            .elocat.identity = elocat.identity,
            .llocat.identity = llocat.identity,

            .kappa = kappa ))),
  weight = eval(substitute(expression({
    ned2l.dlocat2 <- 2 / Scale.w^2
    wz <- cbind(ned2l.dlocat2 * dlocat.deta^2)
    c(w) * wz
  }), list( .Scale.arg = Scale.arg,
            .elocat = elocat, .llocat = llocat,
            .elocat.identity = elocat.identity,
            .llocat.identity = llocat.identity  ))))
}





loglaplace2.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}

 loglaplace2 <- function(tau = NULL,
                         llocation = "loge", lscale = "loge",
                         ilocation = NULL, iscale = NULL,
                         kappa = sqrt(tau/(1-tau)),
                         ishrinkage = 0.95,
                         parallel.locat = FALSE, digt = 4,
                         eq.scale = TRUE,
                         idf.mu = 3,
                         rep0 = 0.5, nsimEIM = NULL,
                         imethod = 1, zero = "(1 + M/2):M") {
 warning("it is best to use loglaplace1()")

  if (length(nsimEIM) &&
     (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 10))
    stop("argument 'nsimEIM' should be an integer greater than 10")
  if (!is.Numeric(rep0, positive = TRUE, length.arg = 1) ||
      rep0 > 1)
    stop("bad input for argument 'rep0'")
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")
  if (length(tau) && max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")




  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")
  if (length(iscale) && !is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")


  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")
  if (!is.logical(eq.scale) || length(eq.scale) != 1)
    stop("bad input for argument 'eq.scale'")
  if (!is.logical(parallel.locat) ||
      length(parallel.locat) != 1)
    stop("bad input for argument 'parallel.locat'")
  fittedMean <- FALSE
  if (!is.logical(fittedMean) || length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")

  if (llocat != "loge")
    stop("argument 'llocat' must be \"loge\"")


  new("vglmff",
  blurb = c("Two-parameter log-Laplace distribution\n\n",
            "Links:      ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale", lscale, earg = escale),
            "\n", "\n",
            "Mean:       zz location + scale * ",
                         "(1/kappa - kappa) / sqrt(2)", "\n",
            "Quantiles:  location", "\n",
            "Variance:   zz scale^2 * (1 + kappa^4) / (2 * kappa^2)"),
  constraints = eval(substitute(expression({
  .ZERO <- .zero
  if (is.character( .ZERO ))
    .ZERO <- eval(parse(text = .ZERO ))
  .PARALLEL <- .parallel.locat
      parelHmat <- if (is.logical( .PARALLEL ) && .PARALLEL )
                   matrix(1, M/2, 1) else diag(M/2)
      scaleHmat <- if (is.logical( .eq.scale ) && .eq.scale )
                   matrix(1, M/2, 1) else diag(M/2)
      mycmatrix <- cbind(rbind(  parelHmat, 0*parelHmat),
                         rbind(0*scaleHmat,   scaleHmat))
      constraints <- cm.VGAM(mycmatrix, x = x,
                             bool = .PARALLEL ,
                             constraints = constraints,
                             apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .ZERO , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)

      if ( .PARALLEL && names(constraints)[1] == "(Intercept)") {
          parelHmat <- diag(M/2)
          mycmatrix <- cbind(rbind(  parelHmat, 0*parelHmat),
                            rbind(0*scaleHmat,   scaleHmat))
          constraints[["(Intercept)"]] <- mycmatrix
      }
      if (is.logical( .eq.scale) && .eq.scale &&
       names(constraints)[1] == "(Intercept)") {
        temp3 <- constraints[["(Intercept)"]]
        temp3 <- cbind(temp3[,1:(M/2)], rbind(0*scaleHmat, scaleHmat))
        constraints[["(Intercept)"]] = temp3
      }
    }), list( .eq.scale = eq.scale, .parallel.locat = parallel.locat,
              .zero = zero ))),
  initialize = eval(substitute(expression({
    extra$kappa <- .kappa
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)



    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$M <- M <- 2 * length(extra$kappa)
    extra$n <- n
    extra$y.names <- y.names <-
      paste("tau = ", round(extra$tau, digits = .digt), sep = "")
    extra$individual = FALSE

    predictors.names <- 
        c(namesof(paste("quantile(", y.names, ")", sep = ""),
                  .llocat , earg = .elocat, tag = FALSE),
          namesof(if (M == 2) "scale" else
                  paste("scale", 1:(M/2), sep = ""),
                  .lscale ,    earg = .escale,    tag = FALSE))
        if (weighted.mean(1 * (y < 0.001), w) >= min(extra$tau))
          stop("sample proportion of 0s > minimum 'tau' value. ",
               "Choose larger values for 'tau'.")

        if (!length(etastart)) {
          if ( .imethod == 1) {
            locat.init.y <- weighted.mean(y, w)
            scale.init <- sqrt(var(y) / 2)
          } else if ( .imethod == 2) {
            locat.init.y <- median(y)
            scale.init <- sqrt(sum(c(w)*abs(y-median(y))) / (sum(w) *2))
          } else if ( .imethod == 3) {
            Fit5 <- vsmooth.spline(x = x[, min(ncol(x), 2)],
                                   y = y, w = w,
                                   df = .idf.mu )
            locat.init.y <- c(predict(Fit5, x = x[, min(ncol(x), 2)])$y)
            scale.init <- sqrt(sum(c(w)*abs(y-median(y))) / (sum(w) *2))
          } else {
            use.this <- weighted.mean(y, w)
            locat.init.y <- (1- .ishrinkage )*y + .ishrinkage * use.this
            scale.init <- sqrt(sum(c(w)*abs(y-median(y ))) / (sum(w) *2))
          }
          locat.init.y <- if (length( .ilocat ))
                           rep( .ilocat , length.out = n) else
                           rep(locat.init.y, length.out = n)
          locat.init.y <- matrix(locat.init.y, n, M/2)
          scale.init <- if (length( .iscale))
                           rep( .iscale, length.out = n) else
                           rep(scale.init, length.out = n)
          scale.init <- matrix(scale.init, n, M/2)
          etastart <-
            cbind(theta2eta(locat.init.y, .llocat , earg = .elocat ),
                  theta2eta(scale.init, .lscale , earg = .escale ))
        }
    }), list( .imethod = imethod,
              .idf.mu = idf.mu, .kappa = kappa,
              .ishrinkage = ishrinkage, .digt = digt,
              .llocat = llocat, .lscale = lscale,
              .elocat = elocat, .escale = escale,
              .ilocat = ilocat, .iscale = iscale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat.y <- eta2theta(eta[,1:(extra$M/2), drop = FALSE],
                               .llocat , earg = .elocat )
    if ( .fittedMean ) {
      kappamat <- matrix(extra$kappa, extra$n, extra$M/2,
                        byrow = TRUE)
      Scale.y <- eta2theta(eta[,(1+extra$M/2):extra$M],
                          .lscale , earg = .escale )
      locat.y + Scale.y * (1/kappamat - kappamat)
    } else {
      dimnames(locat.y) = list(dimnames(eta)[[1]], extra$y.names)
      locat.y
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .fittedMean = fittedMean,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat , scale = .lscale )

    misc$earg <- list(location = .elocat , scale = .escale )

    misc$expected <- TRUE
    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    misc$true.mu <- .fittedMean # @fitted is not a true mu?
    misc$nsimEIM <- .nsimEIM
    misc$rep0 <- .rep0
        extra$percentile <- numeric(length(misc$kappa))
        locat <- as.matrix(locat.y)
        for (ii in 1:length(misc$kappa))
          extra$percentile[ii] <- 100 *
                                 weighted.mean(y <= locat.y[, ii], w)
  }), list( .elocat = elocat, .llocat = llocat,
            .escale = escale, .lscale = lscale,
            .fittedMean = fittedMean,
            .nsimEIM = nsimEIM, .rep0 = rep0,
            .kappa = kappa ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {

    kappamat <- matrix(extra$kappa, extra$n, extra$M/2, byrow = TRUE)
    Scale.w <- eta2theta(eta[, (1+extra$M/2):extra$M],
                         .lscale , earg = .escale )
    ymat <- matrix(y, extra$n, extra$M/2)
    ymat[ymat <= 0] <- min(min(y[y > 0]), .rep0 )  # Adjust for 0s
    ell.mat <- matrix(c(dloglaplace(x = c(ymat),
                                    locat.ald = c(eta[, 1:(extra$M/2)]),
                                    scale.ald = c(Scale.w),
                                    kappa = c(kappamat), log = TRUE)),
                      extra$n, extra$M/2)
      if (residuals) {
        stop("loglikelihood residuals not implemented yet")
      } else {
      ll.elts <- c(w) * ell.mat
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .escale = escale, .lscale = lscale,
           .rep0 = rep0, .kappa = kappa ))),


  vfamily = c("loglaplace2"),
  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M/2)
    Scale.w <- eta2theta(eta[,(1+extra$M/2):extra$M],
                        .lscale , earg = .escale )
    locat.w <- eta[,1:(extra$M/2), drop = FALSE]
    locat.y <- eta2theta(locat.w, .llocat , earg = .elocat )
    kappamat <- matrix(extra$kappa, n, M/2, byrow = TRUE)
    w.mat <- ymat
    w.mat[w.mat <= 0] <- min(min(w.mat[w.mat > 0]), .rep0)  # Adjust for 0s
    w.mat <- theta2eta(w.mat, .llocat , earg = .elocat )  # w.mat=log(w.mat)
    zedd <- abs(w.mat-locat.w) / Scale.w
    dl.dlocat <- sqrt(2) *
                   ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                   sign(w.mat-locat.w) / Scale.w
    dl.dscale <-  sqrt(2) *
                 ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                 zedd / Scale.w - 1 / Scale.w
    dlocat.deta <- dtheta.deta(locat.w, .llocat , earg = .elocat )
    dscale.deta <- dtheta.deta(Scale.w, .lscale , earg = .escale )
    c(w) * cbind(dl.dlocat * dlocat.deta,
                 dl.dscale * dscale.deta)
  }), list( .escale = escale, .lscale = lscale,
            .elocat = elocat, .llocat = llocat,
            .rep0 = rep0, .kappa = kappa ))),
  weight = eval(substitute(expression({
    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    dthetas.detas <- cbind(dlocat.deta, dscale.deta)
    if (length( .nsimEIM )) {
        for (ii in 1:( .nsimEIM )) {
            wsim <- matrix(rloglap(n*M/2, loc = c(locat.w),
                                  sca = c(Scale.w),
                                  kappa = c(kappamat)), n, M/2)
            zedd <- abs(wsim-locat.w) / Scale.w
            dl.dlocat <- sqrt(2) *
                ifelse(wsim >= locat.w, kappamat, 1/kappamat) *
                sign(wsim-locat.w) / Scale.w
            dl.dscale <-  sqrt(2) *
                ifelse(wsim >= locat.w, kappamat, 1/kappamat) *
                zedd / Scale.w - 1 / Scale.w

            rm(wsim)
            temp3 <- cbind(dl.dlocat, dl.dscale)  # n x M matrix
            run.varcov <- ((ii-1) * run.varcov +
               temp3[,ind1$row.index]*temp3[,ind1$col.index]) / ii
        }
        wz <- if (intercept.only)
            matrix(colMeans(run.varcov),
                   n, ncol(run.varcov), byrow = TRUE) else run.varcov

        wz <- wz * dthetas.detas[,ind1$row] * dthetas.detas[,ind1$col]
        wz <- c(w) * matrix(wz, n, dimm(M))
        wz
    } else {
        d2l.dlocat2 <- 2 / (Scale.w * locat.w)^2
        d2l.dscale2 <- 1 / Scale.w^2
        wz <- cbind(d2l.dlocat2 * dlocat.deta^2,
                   d2l.dscale2 * dscale.deta^2)
        c(w) * wz
    }
  }), list( .elocat = elocat, .escale = escale,
            .llocat = llocat, .lscale = lscale,
            .nsimEIM = nsimEIM) )))
}






logitlaplace1.control <- function(maxit = 300, ...) {
    list(maxit = maxit)
}


adjust01.logitlaplace1 <- function(ymat, y, w, rep01) {
    rangey01 <- range(y[(y > 0) & (y < 1)])
    ymat[ymat <= 0] <- min(rangey01[1] / 2,           rep01 / w[y <= 0])
    ymat[ymat >= 1] <- max((1 + rangey01[2]) / 2, 1 - rep01 / w[y >= 1])
    ymat
}





 logitlaplace1 <-
  function(tau = NULL,
           llocation = "logit",
           ilocation = NULL,
           kappa = sqrt(tau/(1-tau)),
           Scale.arg = 1,
           ishrinkage = 0.95, parallel.locat = FALSE, digt = 4,
           idf.mu = 3,
           rep01 = 0.5,
           imethod = 1, zero = NULL) {

  if (!is.Numeric(rep01, positive = TRUE, length.arg = 1) ||
      rep01 > 0.5)
    stop("bad input for argument 'rep01'")
  if (!is.Numeric(kappa, positive = TRUE))
    stop("bad input for argument 'kappa'")

  if (length(tau) && max(abs(kappa - sqrt(tau/(1-tau)))) > 1.0e-6)
    stop("arguments 'kappa' and 'tau' do not match")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")
  ilocat <- ilocation


  llocat.identity <- as.list(substitute("identitylink"))
  elocat.identity <- link2list(llocat.identity)
  llocat.identity <- attr(elocat.identity, "function.name")




  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 4)
    stop("argument 'imethod' must be 1, 2 or ... 4")

  if (!is.Numeric(ishrinkage, length.arg = 1) ||
     ishrinkage < 0 ||
     ishrinkage > 1)
    stop("bad input for argument 'ishrinkage'")

  if (!is.Numeric(Scale.arg, positive = TRUE))
    stop("bad input for argument 'Scale.arg'")
  if (!is.logical(parallel.locat) ||
      length(parallel.locat) != 1)
    stop("bad input for argument 'parallel.locat'")
  fittedMean <- FALSE
  if (!is.logical(fittedMean) ||
      length(fittedMean) != 1)
    stop("bad input for argument 'fittedMean'")


  mystring0 <- namesof("location", llocat, earg = elocat)
  mychars <- substring(mystring0, first = 1:nchar(mystring0),
                      last = 1:nchar(mystring0))
  mychars[nchar(mystring0)] = ", inverse = TRUE)"
  mystring1 <- paste(mychars, collapse = "")




  new("vglmff",
  blurb = c("One-parameter ", llocat, "-Laplace distribution\n\n",
            "Links:      ", mystring0, "\n", "\n",
          "Quantiles:  ", mystring1),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel.locat ,
                           constraints = constraints, apply.int = FALSE)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel.locat = parallel.locat,
            .Scale.arg = Scale.arg, .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         multipleResponses = FALSE,
         parameters.names = c("location"),
         llocation = .llocat ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocat ))),

  initialize = eval(substitute(expression({
    extra$M <- M <- max(length( .Scale.arg ), length( .kappa ))  # Recycle
    extra$Scale <- rep( .Scale.arg, length = M)
    extra$kappa <- rep( .kappa, length = M)
    extra$tau <- extra$kappa^2 / (1 + extra$kappa^2)



    temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y





    extra$n <- n
    extra$y.names <- y.names <-
      paste("tau = ", round(extra$tau, digits = .digt), sep = "")
    extra$individual <- FALSE

    predictors.names <-
        namesof(paste("quantile(", y.names, ")", sep = ""),
                .llocat , earg = .elocat, tag = FALSE)

      if (all(y == 0 | y == 1))
        stop("response cannot be all 0s or 1s")
      if (min(y) < 0)
        stop("negative response values detected")
      if (max(y) > 1)
        stop("response values greater than 1 detected")
      if ((prop.0. <- weighted.mean(1*(y == 0), w)) >= min(extra$tau))
        stop("sample proportion of 0s == ", round(prop.0., digits = 4),
             " > minimum 'tau' value. Choose larger values for 'tau'.")
      if ((prop.1. <- weighted.mean(1*(y == 1), w)) >= max(extra$tau))
        stop("sample proportion of 1s == ", round(prop.1., digits = 4),
             " < maximum 'tau' value. Choose smaller values for 'tau'.")
      if (!length(etastart)) {
        if ( .imethod == 1) {
          locat.init <- quantile(rep(y, w), probs= extra$tau)
        } else if ( .imethod == 2) {
          locat.init <- weighted.mean(y, w)
          locat.init <- median(rep(y, w))
        } else if ( .imethod == 3) {
          use.this <- weighted.mean(y, w)
          locat.init <- (1- .ishrinkage )*y + use.this * .ishrinkage
        } else {
          stop("this option not implemented")
        }


      locat.init <- if (length( .ilocat ))
                       rep( .ilocat , length.out = M) else
                       rep(locat.init, length.out = M)
      locat.init <- matrix(locat.init, n, M, byrow = TRUE)
      locat.init <- abs(locat.init)
      etastart <-
          cbind(theta2eta(locat.init, .llocat , earg = .elocat ))
    }
  }), list( .imethod = imethod,
            .idf.mu = idf.mu,
            .ishrinkage = ishrinkage, .digt = digt,
            .elocat = elocat, .Scale.arg = Scale.arg,
            .llocat = llocat, .kappa = kappa,
            .ilocat = ilocat ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    locat.y <- eta2theta(eta, .llocat , earg = .elocat )
    if ( .fittedMean ) {
      stop("Yet to do: handle 'fittedMean = TRUE'")
      kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
      Scale <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
      locat.y + Scale * (1/kappamat - kappamat)
    } else {
      if (length(locat.y) > extra$n)
        dimnames(locat.y) <- list(dimnames(eta)[[1]], extra$y.names)
      locat.y
      }
  }, list( .elocat = elocat, .llocat = llocat,
           .fittedMean = fittedMean, .Scale.arg = Scale.arg,
           .kappa = kappa ))),
  last = eval(substitute(expression({
    misc$link <-    c(location = .llocat )
    misc$earg <- list(location = .elocat )

    misc$expected <- TRUE

    extra$kappa <- misc$kappa <- .kappa
    extra$tau <- misc$tau <- misc$kappa^2 / (1 + misc$kappa^2)
    extra$Scale.arg <- .Scale.arg

    misc$true.mu <- .fittedMean # @fitted is not a true mu?
    misc$rep01 <- .rep01

    extra$percentile <- numeric(length(misc$kappa))
    locat.y <- eta2theta(eta, .llocat , earg = .elocat )
    locat.y <- as.matrix(locat.y)
    for (ii in 1:length(misc$kappa))
      extra$percentile[ii] <- 100 *
                             weighted.mean(y <= locat.y[, ii], w)

  }), list( .elocat = elocat, .llocat = llocat,
            .Scale.arg = Scale.arg, .fittedMean = fittedMean,
            .rep01 = rep01,
            .kappa = kappa ))),


  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    kappamat <- matrix(extra$kappa, extra$n, extra$M, byrow = TRUE)
    Scale.w  <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    ymat     <- matrix(y,           extra$n, extra$M)
    ymat <- adjust01.logitlaplace1(ymat = ymat, y = y, w = w,
                                   rep01 = .rep01)
    w.mat <- theta2eta(ymat, .llocat , earg = .elocat )  # e.g., logit()
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dalap(x = c(w.mat), location = c(eta),
                     scale = c(Scale.w), kappa = c(kappamat),
                     log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .elocat = elocat, .llocat = llocat,
           .rep01 = rep01,
           .Scale.arg = Scale.arg, .kappa = kappa ))),


  vfamily = c("logitlaplace1"),
  deriv = eval(substitute(expression({
    ymat <- matrix(y, n, M)
    Scale.w <- matrix(extra$Scale, extra$n, extra$M, byrow = TRUE)
    locat.w <- eta
    kappamat <- matrix(extra$kappa, n, M, byrow = TRUE)
    ymat <- adjust01.logitlaplace1(ymat = ymat, y = y, w = w,
                                  rep01 = .rep01)
    w.mat <- theta2eta(ymat, .llocat , earg = .elocat )  # e.g., logit()
    zedd <- abs(w.mat-locat.w) / Scale.w
    dl.dlocat <- ifelse(w.mat >= locat.w, kappamat, 1/kappamat) *
                   sqrt(2) * sign(w.mat-locat.w) / Scale.w


    dlocat.deta <- dtheta.deta(locat.w,
                              "identitylink",
                              earg = .elocat.identity )


    c(w) * cbind(dl.dlocat * dlocat.deta)
  }), list( .Scale.arg = Scale.arg, .rep01 = rep01,
            .elocat = elocat,
            .llocat = llocat,

            .elocat.identity = elocat.identity,
            .llocat.identity = llocat.identity,

            .kappa = kappa ))),
  weight = eval(substitute(expression({
    d2l.dlocat2 <- 2 / Scale.w^2
    wz <- cbind(d2l.dlocat2 * dlocat.deta^2)
    c(w) * wz
  }), list( .Scale.arg = Scale.arg,
            .elocat = elocat, .llocat = llocat ))))
}





















