# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.















 cens.poisson <- function(link = "loge", imu = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Censored Poisson distribution\n\n",
            "Link:     ", namesof("mu", link, earg = earg), "\n",
            "Variance: mu"),
  initialize = eval(substitute(expression({
    if (any(is.na(y)))
      stop("NAs are not allowed in the response")


    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 3,
              Is.integer.y = TRUE)


    centype <- attr(y, "type")

    if (centype == "right") {
        temp <- y[, 2]
        extra$uncensored <- ifelse(temp == 1, TRUE, FALSE)
        extra$rightcensored <- ifelse(temp == 0, TRUE, FALSE)
        extra$leftcensored <- rep(FALSE, len = n)
        extra$interval <- rep(FALSE, len = n)
        init.mu <- pmax(y[, 1], 1/8)
    } else
    if (centype == "left") {
        temp <- y[, 2]
        extra$uncensored <- ifelse(temp == 1, TRUE, FALSE)
        extra$rightcensored <- rep(FALSE, len = n)
        extra$leftcensored <- ifelse(temp == 0, TRUE, FALSE)
        extra$interval <- rep(FALSE, len = n)
        init.mu <- pmax(y[, 1], 1/8)
    } else
    if (centype == "interval" ||
        centype == "interval2") {
        temp <- y[, 3]
        extra$uncensored <- ifelse(temp == 1, TRUE, FALSE)
        extra$rightcensored <- ifelse(temp == 0, TRUE, FALSE)
        extra$leftcensored <- ifelse(temp == 2, TRUE, FALSE)
        extra$intervalcensored <- ifelse(temp == 3, TRUE, FALSE)
        init.mu <- pmax((y[, 1] + y[, 2])/2, 1/8)  # for intervalcensored
        if (any(extra$uncensored))
        init.mu[extra$uncensored] <- pmax(y[extra$uncensored, 1], 1/8)
        if (any(extra$rightcensored))
     init.mu[extra$rightcensored] <- pmax(y[extra$rightcensored, 1], 1/8)
        if (any(extra$leftcensored))
       init.mu[extra$leftcensored] <- pmax(y[extra$leftcensored, 1], 1/8)
    } else
    if (centype == "counting") {
      stop("type == 'counting' not compatible with cens.poisson()")
      init.mu <- pmax(y[, 1], 1/8)
      stop("currently not working")
    } else
      stop("response have to be in a class of SurvS4")

      if (length( .imu )) init.mu <- 0 * y[, 1] + .imu
  
      predictors.names <-
        namesof("mu", .link, earg = .earg, short = TRUE)

      if (!length(etastart))
        etastart <- theta2eta(init.mu, link = .link, earg = .earg)
  }), list( .link = link, .earg = earg, .imu = imu))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mu <- eta2theta(eta, link = .link, earg = .earg)
    mu
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$expected <- FALSE

    misc$link <-    c("mu" = .link)

    misc$earg <- list("mu" = .earg)

    misc$multipleResponses <- FALSE
  }), list( .link = link, .earg = earg ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    theta2eta(mu, link = .link, earg = .earg)
  }, list( .link = link, .earg = earg ))),
  loglikelihood = function(mu, y, w, residuals = FALSE, eta,
                           extra = NULL) {
    cen0 <- extra$uncensored
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cenI <- extra$intervalcensored
    if (residuals){
      stop("loglikelihood residuals not implemented yet")
    } else {
      sum(w[cen0] * dpois(y[cen0, 1], mu[cen0], log = TRUE)) +
      sum(w[cenU] * log1p(-ppois(y[cenU, 1] - 1, mu[cenU]))) +
      sum(w[cenL] * ppois(y[cenL, 1] - 1, mu[cenL], log.p = TRUE)) +
      sum(w[cenI] * log(ppois(y[cenI, 2], mu[cenI]) -
                        ppois(y[cenI, 1], mu[cenI])))
    }
  },
  vfamily = "cens.poisson",
  deriv = eval(substitute(expression({
    cen0 <- extra$uncensored
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cenI <- extra$intervalcensored
      lambda <- eta2theta(eta, link = .link, earg = .earg)

      dl.dlambda <- (y[, 1] - lambda)/lambda   # uncensored

      yllim <- yulim <- y[, 1]   # uncensored

      if (any(cenU)) {
        yllim[cenU] <- y[cenU, 1]
        densm1 <- dpois(yllim-1, lambda)
        queue <- ppois(yllim-1, lambda, lower.tail = FALSE)
        dl.dlambda[cenU] <- densm1[cenU] / queue[cenU]
      }
      if (any(cenL)) {
          yulim[cenL] <- y[cenL, 1] - 1
          densm0 <- dpois(yulim, lambda)
          Queue <- ppois(yulim, lambda)    # Left tail probability
          dl.dlambda[cenL] <- -densm0[cenL] / Queue[cenL]
      }
      if (any(cenI)) {
          yllim[cenI] <- y[cenI, 1] + 1
          yulim[cenI] <- y[cenI, 2]
          Queue1 <- ppois(yllim-1, lambda)
          Queue2 <- ppois(yulim, lambda)
          densm02 <- dpois(yulim, lambda)
          densm12 <- dpois(yllim-1, lambda)
          dl.dlambda[cenI] <-
              (-densm02[cenI]+densm12[cenI]) / (Queue2[cenI]-Queue1[cenI])
      }

      dlambda.deta <- dtheta.deta(theta=lambda, link =  .link, earg = .earg)

    c(w) * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    d2lambda.deta2 <- d2theta.deta2(theta = lambda,
                                    link = .link, earg = .earg )
    d2l.dlambda2 <- 1 / lambda # uncensored; Fisher scoring

    if (any(cenU)) {
      densm2 <- dpois(yllim-2, lambda)
      d2l.dlambda2[cenU] <- (dl.dlambda[cenU])^2 -
          (densm2[cenU]-densm1[cenU])/queue[cenU]
    }
    if (any(cenL)) {
      densm1 <- dpois(yulim-1, lambda)
      d2l.dlambda2[cenL] <- (dl.dlambda[cenL])^2 -
          (densm0[cenL]-densm1[cenL])/Queue[cenL]
    }
    if (any(cenI)) {
      densm03 <- dpois(yulim-1, lambda)
      densm13 <- dpois(yllim-2, lambda)
      d2l.dlambda2[cenI] <- (dl.dlambda[cenI])^2 -
          (densm13[cenI]-densm12[cenI]-densm03[cenI] +
           densm02[cenI]) / (Queue2[cenI]-Queue1[cenI])
    }
    wz <-  c(w) * ((dlambda.deta^2) * d2l.dlambda2)
    wz
  }), list( .link = link, .earg = earg ))))
}




if (FALSE)
 cens.exponential <- 
 ecens.exponential <- function(link = "loge", location = 0) {
  if (!is.Numeric(location, length.arg = 1))
    stop("bad input for 'location'")

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Censored exponential distribution\n\n",
            "Link:     ", namesof("rate", link, tag = TRUE), "\n",
            "Mean:     ", "mu = ", location, " + 1 / ",
            namesof("rate", link, tag = FALSE), "\n",
            "Variance: ",
            if (location == 0) "Exponential: mu^2" else
            paste("(mu-",  location, ")^2", sep = "")),
  initialize = eval(substitute(expression({
    extra$location <- .location

    if (any(y[, 1] <= extra$location))
        stop("all responses must be greater than ", extra$location)

    predictors.names <- namesof("rate", .link , .earg , tag = FALSE)

    type <- attr(y, "type")
    if (type == "right" || type == "left"){
      mu <- y[, 1] + (abs(y[, 1] - extra$location) < 0.001) / 8
    } else if (type == "interval") {
      temp <- y[, 3]
      mu <- ifelse(temp == 3,
                   y[, 2] + (abs(y[, 2] - extra$location) < 0.001) / 8,
                   y[, 1] + (abs(y[, 1] - extra$location) < 0.001) / 8)
    }
    if (!length(etastart))
      etastart <- theta2eta(1/(mu-extra$location), .link , .earg )

    if (type == "right") {
      temp <- y[, 2]
      extra$uncensored <- ifelse(temp == 1, TRUE, FALSE)
      extra$rightcensored <- ifelse(temp == 0, TRUE, FALSE)
      extra$leftcensored <- rep(FALSE, len = n)
      extra$interval <- rep(FALSE, len = n)
    } else
    if (type == "left") {
      temp <- y[, 2]
      extra$uncensored <- ifelse(temp == 1, TRUE, FALSE)
      extra$rightcensored <- rep(FALSE, len = n)
      extra$leftcensored <- ifelse(temp == 0, TRUE, FALSE)
      extra$interval <- rep(FALSE, len = n)
    } else
    if (type == "counting") {
      stop("type == 'counting' not recognized")
      extra$uncensored <- rep(temp == 1, TRUE, FALSE)
      extra$interval <- rep(FALSE, len = n)
      extra$leftcensored <- rep(FALSE, len = n)
      extra$rightcensored <- rep(FALSE, len = n)
      extra$counting <- ifelse(temp == 0, TRUE, FALSE)
    } else
    if (type == "interval") {
      temp <- y[, 3]
      extra$uncensored <- ifelse(temp == 1, TRUE, FALSE)
      extra$rightcensored <- ifelse(temp == 0, TRUE, FALSE)
      extra$leftcensored <- ifelse(temp == 2, TRUE, FALSE)
      extra$interval <- ifelse(temp == 3, TRUE, FALSE)
    } else
      stop("'type' not recognized")
  }), list( .location = location, .link = link ))),
  linkinv = eval(substitute(function(eta, extra = NULL)
    extra$location + 1 / eta2theta(eta, .link , .earg ),
  list( .link = link ) )),
  last = eval(substitute(expression({
    misc$location <- extra$location
    misc$link   <- c("rate" = .link)
    misc$multipleResponses <- FALSE
  }), list( .link = link ))),
  link = eval(substitute(function(mu, extra = NULL)
    theta2eta(1 / (mu - extra$location), .link , .earg ),
  list( .link = link ) )),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    rate <- 1 / (mu - extra$location)
    cen0 <- extra$uncensored
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cenI <- extra$interval
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(w[cenL] * log1p(-exp(-rate[cenL] *
                  (y[cenL, 1] - extra$location)))) +
    sum(w[cenU] * (-rate[cenU]*(y[cenU, 1]-extra$location))) +
    sum(w[cen0] * (log(rate[cen0]) -
                   rate[cen0]*(y[cen0, 1]-extra$location))) +
    sum(w[cenI] * log(-exp(-rate[cenI]*(y[cenI, 2]-extra$location))+
    exp(-rate[cenI]*(y[cenI, 1]-extra$location))))
  }, list( .link = link ))),
  vfamily = c("ecens.exponential"),
  deriv = eval(substitute(expression({
    rate <- 1 / (mu - extra$location)
    cen0 <- extra$uncensored
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cenI <- extra$interval
    dl.drate <- 1/rate - (y[, 1]-extra$location)  # uncensored
    tmp200 <- exp(-rate*(y[, 1]-extra$location))
    tmp200b <- exp(-rate*(y[, 2]-extra$location))  # for interval censored
    if (any(cenL))
      dl.drate[cenL] <- (y[cenL, 1]-extra$location) *
                        tmp200[cenL] / (1 - tmp200[cenL])
    if (any(cenU))
      dl.drate[cenU] <- -(y[cenU, 1]-extra$location)
    if (any(cenI))
      dl.drate[cenI] <- ((y[cenI, 2] - extra$location) *
                        tmp200b[cenI] - (y[cenI, 1] - extra$location) *
                        tmp200[cenI]) / (-tmp200b[cenI] + tmp200[cenI])

    drate.deta <- dtheta.deta(rate, .link , .earg )

    c(w) * dl.drate * drate.deta
  }), list( .link = link ) )),
  weight = eval(substitute(expression({
    A123 <- ((mu-extra$location)^2)  # uncensored d2l.drate2
    Lowpt <- ifelse(cenL, y[, 1], extra$location)
    Lowpt <- ifelse(cenI, y[, 1], Lowpt)  #interval censored
    Upppt <- ifelse(cenU, y[, 1], Inf)
    Upppt <- ifelse(cenI, y[, 2], Upppt)  #interval censored
    tmp300 <- exp(-rate*(Lowpt - extra$location))

    d2l.drate2 <- 0 * y[, 1]
    ind50 <- Lowpt > extra$location

    d2l.drate2[ind50] <- (Lowpt[ind50]-extra$location)^2 *
                        tmp300[ind50] / (1-tmp300[ind50])
    d2l.drate2 <- d2l.drate2 + (exp(-rate*(Lowpt-extra$location)) -
                               exp(-rate*(Upppt-extra$location))) * A123

    wz <- c(w) * (drate.deta^2) * d2l.drate2
    wz
    }), list( .link = link ))))
}





 cennormal <-
 cens.normal <- function(lmu = "identitylink", lsd = "loge",
                         imethod = 1, zero = "sd") {


  lmu <- as.list(substitute(lmu))
  emu <- link2list(lmu)
  lmu <- attr(emu, "function.name")

  lsd <- as.list(substitute(lsd))
  esd <- link2list(lsd)
  lsd <- attr(esd, "function.name")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
    imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  new("vglmff",
  blurb = c("Censored univariate normal\n\n",
            "Links:    ", namesof("mu", lmu, tag = TRUE), "; ",
                          namesof("sd", lsd, tag = TRUE), "\n",
            "Conditional variance: sd^2"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         zero = .zero ,
         multiple.responses = FALSE,
         parameters.names = c("mu", "sd"),
         expected = TRUE )
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    if (!length(extra$leftcensored))
      extra$leftcensored <- rep(FALSE, len = n)
    if (!length(extra$rightcensored))
      extra$rightcensored <- rep(FALSE, len = n)
    if (any(extra$rightcensored & extra$leftcensored))
        stop("some observations are both right and left censored!")

    predictors.names <-
      c(namesof("mu", .lmu , earg = .emu , tag = FALSE),
        namesof("sd", .lsd , earg = .esd , tag = FALSE))

    if (!length(etastart)) {
      anyc <- extra$leftcensored | extra$rightcensored
      i11 <- if ( .imethod == 1) anyc else FALSE  # can be all data
      junk <- lm.wfit(x = cbind(x[!i11, ]),
                      y = y[!i11], w = w[!i11])
      sd.y.est <- sqrt(sum(w[!i11] * junk$resid^2) / junk$df.residual)
      etastart <- cbind(mu = y,
                        rep(theta2eta(sd.y.est, .lsd), length = n))
      if (any(anyc))
        etastart[anyc, 1] <- x[anyc, , drop = FALSE] %*% junk$coeff
    }
 }), list( .lmu = lmu, .lsd = lsd,
           .emu = emu, .esd = esd,
           .imethod = imethod ))),
  linkinv = eval(substitute( function(eta, extra = NULL) {
    eta2theta(eta[, 1], .lmu , earg = .emu )
  }, list( .lmu = lmu, .emu = emu ))),
  last = eval(substitute(expression({
    misc$link <-    c("mu" = .lmu , "sd" = .lsd )

    misc$earg <- list("mu" = .emu , "sd" = .esd )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lmu = lmu, .lsd = lsd,
            .emu = emu, .esd = esd ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cen0 <- !cenL & !cenU   # uncensored obsns

    mum <- eta2theta(eta[, 1], .lmu , earg = .emu )
    sdv <- eta2theta(eta[, 2], .lsd , earg = .esd )

    Lower <- ifelse(cenL, y, -Inf)
    Upper <- ifelse(cenU, y,  Inf)
    ell1 <- -log(sdv[cen0]) - 0.5 * ((y[cen0] - mum[cen0])/sdv[cen0])^2
    ell2 <- log1p(-pnorm((mum[cenL] - Lower[cenL]) / sdv[cenL]))
    ell3 <- log1p(-pnorm(( Upper[cenU] -  mum[cenU]) / sdv[cenU]))
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
    sum(w[cen0] * ell1) +
    sum(w[cenL] * ell2) +
    sum(w[cenU] * ell3)
  }, list( .lmu = lmu, .lsd = lsd,
           .emu = emu, .esd = esd ))),
  vfamily = c("cens.normal"),
  deriv = eval(substitute(expression({
    cenL <- extra$leftcensored
    cenU <- extra$rightcensored
    cen0 <- !cenL & !cenU   # uncensored obsns
    Lower <- ifelse(cenL, y, -Inf)
    Upper <- ifelse(cenU, y,  Inf)

    mum <- eta2theta(eta[, 1], .lmu)
    sdv <- eta2theta(eta[, 2], .lsd)

    dl.dmu <- (y-mum) / sdv^2
    dl.dsd <- (((y-mum)/sdv)^2 - 1) / sdv

    dmu.deta <- dtheta.deta(mum, .lmu , earg = .emu )
    dsd.deta <- dtheta.deta(sdv, .lsd , earg = .esd )

    if (any(cenL)) {
      mumL <- mum - Lower
      temp21L <- mumL[cenL] / sdv[cenL]
      PhiL <- pnorm(temp21L)
      phiL <- dnorm(temp21L)
      fred21 <- phiL / (1 - PhiL)
      dl.dmu[cenL] <- -fred21 / sdv[cenL]
      dl.dsd[cenL] <- mumL[cenL] * fred21 / sdv[cenL]^2
      rm(fred21)
    }
    if (any(cenU)) {
      mumU <- Upper - mum
      temp21U <- mumU[cenU] / sdv[cenU]
      PhiU <- pnorm(temp21U)
      phiU <- dnorm(temp21U)
      fred21 <- phiU / (1 - PhiU)
      dl.dmu[cenU] <- fred21 / sdv[cenU]  # Negated
      dl.dsd[cenU] <- mumU[cenU] * fred21 / sdv[cenU]^2
      rm(fred21)
    }
    c(w) * cbind(dl.dmu * dmu.deta,
                 dl.dsd * dsd.deta)
  }), list( .lmu = lmu, .lsd = lsd,
            .emu = emu, .esd = esd ))),
  weight = eval(substitute(expression({
    A1 <- 1 - pnorm((mum - Lower) / sdv)  # Lower
    A3 <- 1 - pnorm((Upper - mum) / sdv)  # Upper
    A2 <- 1 - A1 - A3                     # Middle; uncensored
    wz <- matrix(0, n, 3)
    wz[, iam(1, 1,M)] <- A2 * 1 / sdv^2  # ed2l.dmu2
    wz[, iam(2, 2,M)] <- A2 * 2 / sdv^2  # ed2l.dsd2
    mumL <- mum - Lower
    temp21L <- mumL / sdv
    PhiL <- pnorm(temp21L)
    phiL <- dnorm(temp21L)
    temp31L <- ((1-PhiL) * sdv)^2 
    wz.cenL11 <- phiL * (phiL - (1-PhiL)*temp21L) / temp31L
    wz.cenL22 <- mumL * phiL * ((1-PhiL) * (2 - temp21L^2) +
                 mumL * phiL / sdv) / (sdv * temp31L)
    wz.cenL12 <- phiL * ((1-PhiL)*(temp21L^2 - 1) -
                 temp21L*phiL) / temp31L
    wz.cenL11[!is.finite(wz.cenL11)] <- 0
    wz.cenL22[!is.finite(wz.cenL22)] <- 0
    wz.cenL12[!is.finite(wz.cenL12)] <- 0
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] + A1 * wz.cenL11
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] + A1 * wz.cenL22
    wz[, iam(1, 2, M)] <- A1 * wz.cenL12
    mumU <- Upper - mum    # often Inf
    temp21U <- mumU / sdv    # often Inf
    PhiU <- pnorm(temp21U)  # often 1
    phiU <- dnorm(temp21U)  # often 0
    temp31U <- ((1-PhiU) * sdv)^2  # often 0
    tmp8 <- (1-PhiU)*temp21U
    wzcenU11 <- phiU * (phiU - tmp8) / temp31U
    tmp9 <- (1-PhiU) * (2 - temp21U^2)
    wzcenU22 <- mumU * phiU * (tmp9 + mumU * phiU / sdv) / (sdv * temp31U)
    wzcenU12 <- -phiU * ((1-PhiU)*(temp21U^2 - 1) -
                 temp21U*phiU) / temp31U
    wzcenU11[!is.finite(wzcenU11)] <- 0  # Needed when Upper==Inf
    wzcenU22[!is.finite(wzcenU22)] <- 0  # Needed when Upper==Inf
    wzcenU12[!is.finite(wzcenU12)] <- 0  # Needed when Upper==Inf
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] + A3 * wzcenU11
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] + A3 * wzcenU22
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)] + A3 * wzcenU12
    wz[, iam(1, 1, M)] <- wz[, iam(1, 1, M)] * dmu.deta^2
    wz[, iam(2, 2, M)] <- wz[, iam(2, 2, M)] * dsd.deta^2
    wz[, iam(1, 2, M)] <- wz[, iam(1, 2, M)] * dmu.deta * dsd.deta
    c(w) * wz
  }), list( .lmu = lmu, .lsd = lsd ))))
}




 cens.rayleigh <- function(lscale = "loge",
                           oim  = TRUE) {

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (!is.logical(oim) || length(oim) != 1)
    stop("bad input for argument 'oim'")

  new("vglmff",
  blurb = c("Censored Rayleigh distribution\n\n",
            "f(y) = y*exp(-0.5*(y/scale)^2)/scale^2, y>0, scale>0\n",
            "Link:    ",
            namesof("scale", lscale, earg = escale ), "\n", "\n",
            "Mean:    scale * sqrt(pi / 2)"),
  initialize = eval(substitute(expression({
    if (ncol(cbind(y)) != 1)
      stop("response must be a vector or a one-column matrix")

    if (length(extra$leftcensored))
      stop("cannot handle left-censored data")

    if (!length(extra$rightcensored))
      extra$rightcensored <- rep(FALSE, len = n)

    predictors.names <-
      namesof("scale", .lscale , earg = .escale , tag = FALSE)

    if (!length(etastart)) {
      a.init <- (y+1/8) / sqrt(pi/2)
      etastart <- theta2eta(a.init, .lscale , earg = .escale )
    }
  }), list( .lscale = lscale, .escale = escale ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta, .lscale , earg = .escale )
    Scale * sqrt(pi/2)
  }, list( .lscale = lscale, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link <-    c("scale" = .lscale)

    misc$earg <- list("scale" = .escale)

    misc$oim <- .oim
  }), list( .lscale = lscale, .escale = escale,
            .oim = oim ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Scale <- eta2theta(eta, .lscale , earg = .escale )

    cen0 <- !extra$rightcensored   # uncensored obsns
    cenU <- extra$rightcensored

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else
      sum(w[cen0] * (log(y[cen0]) - 2*log(Scale[cen0]) -
                     0.5*(y[cen0]/Scale[cen0])^2)) -
      sum(w[cenU] * (y[cenU]/Scale[cenU])^2) * 0.5
  }, list( .lscale = lscale,
           .escale = escale ))),
  vfamily = c("cens.rayleigh"),
  deriv = eval(substitute(expression({
    cen0 <- !extra$rightcensored   # uncensored obsns
    cenU <- extra$rightcensored

    Scale <- eta2theta(eta, .lscale , earg = .escale )

    dl.dScale <- ((y/Scale)^2 - 2) / Scale

    dScale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dl.dScale[cenU] <- y[cenU]^2 / Scale[cenU]^3

    c(w) * dl.dScale * dScale.deta
  }), list( .lscale = lscale,
            .escale = escale ))),
  weight = eval(substitute(expression({
    ned2l.dScale2 <- 4 / Scale^2
    wz <- dScale.deta^2 * ned2l.dScale2

    if ( .oim ) {
      d2l.dScale2 <- 3 * (y[cenU])^2 / (Scale[cenU])^4
      d2Scale.deta2 <- d2theta.deta2(Scale[cenU], .lscale , earg = .escale )
      wz[cenU] <- (dScale.deta[cenU])^2 * d2l.dScale2 -
                   dl.dScale[cenU] * d2Scale.deta2
    } else {
      ned2l.dScale2[cenU] <- 6 / (Scale[cenU])^2
      wz[cenU] <- (dScale.deta[cenU])^2 * ned2l.dScale2[cenU]
    }

    c(w) * wz
  }), list( .lscale = lscale, .escale = escale,
            .oim = oim ))))
}










 weibull.mean <-
  function(lmean = "loge", lshape = "loge",
           imean = NULL,   ishape = NULL,
           probs.y = c(0.2, 0.5, 0.8),
           imethod = 1, zero = "shape") {




  imeann <- imean


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lmeann <- as.list(substitute(lmean))
  emeann <- link2list(lmeann)
  lmeann <- attr(emeann, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(probs.y, positive  = TRUE) ||
      length(probs.y) < 2 ||
      max(probs.y) >= 1)
    stop("bad input for argument 'probs.y'")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")

  if (length(imeann))
    if (!is.Numeric(imeann, positive = TRUE))
      stop("argument 'imean' values must be positive")

  blurb.vec <- c(namesof("mean",  lmeann, earg = emeann),
                 namesof("shape", lshape, earg = eshape))

  new("vglmff",
  blurb = c("Weibull distribution (parameterized by the mean)\n\n",
            "Links:    ",
            blurb.vec[1], ", ",
            blurb.vec[2], "\n",
            "Mean:     mean\n",
            "Variance: mean^2 * (gamma(1 + 2/shape) / ",
                      "gamma(1 + 1/shape)^2 - 1)"),
 constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero,
            .lmeann = lmeann ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("mean", "shape"),
         lmean  = .lmeann ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero,
           .lmeann = lmeann, .lshape = lshape ))),

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


    if (is.SurvS4(y))
      stop("only uncensored observations are allowed; ",
           "don't use SurvS4()")


    mynames1 <- param.names("mean" , ncoly)
    mynames2 <- param.names("shape", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lmeann , earg = .emeann , tag = FALSE),
          namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]


    Meann.init <- matrix(if (length( .imeann )) .imeann else 0.5 * colMeans(y),
                         n, ncoly, byrow = TRUE) + 0.5 * y
    Shape.init <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                         n, ncoly, byrow = TRUE)

    if (!length(etastart)) {
      if (!length( .ishape ) ||
          !length( .imeann )) {
        for (ilocal in 1:ncoly) {

          anyc <- FALSE  # extra$leftcensored | extra$rightcensored
          i11 <- if ( .imethod == 1) anyc else FALSE  # Can be all data
          probs.y <- .probs.y
          xvec <- log(-log1p(-probs.y))
          fit0 <- lsfit(x  = xvec,
                        y  = log(quantile(y[!i11, ilocal],
                                 probs = probs.y )))


          if (!is.Numeric(Shape.init[, ilocal]))
            Shape.init[, ilocal] <- 1 / fit0$coef["X"]
        }  # ilocal

        etastart <- 
          cbind(theta2eta(Meann.init, .lmeann , earg = .emeann ),
                theta2eta(Shape.init, .lshape , earg = .eshape ))[,
                interleave.VGAM(M, M1 = M1)]
      }
    }
  }), list( .lmeann = lmeann, .lshape = lshape,
            .emeann = emeann, .eshape = eshape,
            .imeann = imeann, .ishape = ishape,
            .probs.y = probs.y,
            .imethod = imethod ) )),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Meann <- eta2theta(eta[, c(TRUE, FALSE)], .lmeann , earg = .emeann )
    Meann
  }, list( .lmeann = lmeann, .lshape = lshape,
           .emeann = emeann, .eshape = eshape ) )),
  last = eval(substitute(expression({
    regnotok <- any(Shape <= 2)
    if (any(Shape <= 1)) {
      warning("MLE regularity conditions are violated",
              "(shape <= 1) at the final iteration: ",
              "MLEs are not consistent")
    } else if (any(1 < Shape & Shape < 2)) {
      warning("MLE regularity conditions are violated",
              "(1 < shape < 2) at the final iteration: ",
              "MLEs exist but are not asymptotically normal")
    } else if (any(2 == Shape)) {
      warning("MLE regularity conditions are violated",
              "(shape == 2) at the final iteration: ",
              "MLEs exist and are normal and asymptotically ",
              "efficient but with a slower convergence rate than when ",
              "shape > 2")
    }



    M1 <- extra$M1
    avector <- c(rep( .lmeann , length = ncoly),
                 rep( .lshape , length = ncoly))
    misc$link <- avector[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .emeann
      misc$earg[[M1*ii  ]] <- .eshape
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE


    misc$RegCondOK <- !regnotok # Save this for later
    misc$expected <- TRUE   # all(cen0)
  }), list( .lmeann = lmeann, .lshape = lshape,
            .emeann = emeann, .eshape = eshape,
            .imethod = imethod ) )),
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Meann <- eta2theta(eta[, c(TRUE, FALSE)], .lmeann , earg = .emeann )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dweibull(x = y, shape = Shape,
                          scale = Meann / gamma(1 + 1/Shape),
                          log = TRUE))
    }
  }, list( .lmeann = lmeann, .lshape = lshape,
           .emeann = emeann, .eshape = eshape ) )),
  vfamily = c("weibull.mean"),
  deriv = eval(substitute(expression({
    M1 <- 2
    Meann <- eta2theta(eta[, c(TRUE, FALSE)], .lmeann , earg = .emeann )
    Shape <- eta2theta(eta[, c(FALSE, TRUE)], .lshape , earg = .eshape )

    if (FALSE) {
    } else {
      EulerM <- -digamma(1.0)
      AA <- (EulerM - 1)^2 + (pi^2) / 6
      BB <- digamma(1 + 1/Shape)
      CC <- y * gamma(1 + 1/Shape) / Meann
      dl.dmeann <- (CC^Shape - 1) * Shape / Meann  # Agrees
      dl.dshape <- 1/Shape -
                   (log(y/Meann) + lgamma(1 + 1/Shape)) * (CC^Shape - 1) + 
                   (BB / Shape) * (CC^Shape - 1)
    }


    dmeann.deta <- dtheta.deta(Meann, .lmeann , earg = .emeann )
    dshape.deta <- dtheta.deta(Shape, .lshape , earg = .eshape )

    myderiv <- c(w) * cbind(dl.dmeann * dmeann.deta,
                            dl.dshape * dshape.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lmeann = lmeann, .lshape = lshape,
            .emeann = emeann, .eshape = eshape ) )),
  weight = eval(substitute(expression({

    if (FALSE) {
    } else {
      ned2l.dmeann <- (Shape / Meann)^2  #
      ned2l.dshape <- AA / Shape^2  # Unchanged
      ned2l.dshapemeann <- (EulerM - 1 + BB) / Meann
    }


    wz <- array(c(c(w) * ned2l.dmeann * dmeann.deta^2,
                  c(w) * ned2l.dshape * dshape.deta^2,
                  c(w) * ned2l.dshapemeann * dmeann.deta * dshape.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .eshape = eshape ))))
}






 weibullR <-
  function(lscale = "loge", lshape = "loge",
           iscale = NULL,   ishape = NULL,
           lss = TRUE,
           nrfs = 1,
           probs.y = c(0.2, 0.5, 0.8),
           imethod = 1, zero = "shape") {



  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(probs.y, positive  = TRUE) ||
      length(probs.y) < 2 ||
      max(probs.y) >= 1)
    stop("bad input for argument 'probs.y'")


  if (!is.Numeric(nrfs, length.arg = 1) ||
      nrfs < 0 ||
      nrfs > 1)
    stop("bad input for argument 'nrfs'")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE))
      stop("argument 'ishape' values must be positive")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
      stop("argument 'iscale' values must be positive")

  scale.TF <- if (lss) c(TRUE, FALSE) else c(FALSE, TRUE)
  scale.12 <- if (lss) 1:2 else 2:1
  blurb.vec <- c(namesof("scale", lscale, earg = escale),
                 namesof("shape", lshape, earg = eshape))
  blurb.vec <- blurb.vec[scale.12]

  new("vglmff",
  blurb = c("Weibull distribution\n\n",
            "Links:    ",
            blurb.vec[1], ", ",
            blurb.vec[2], "\n",
            "Mean:     scale * gamma(1 + 1/shape)\n",
            "Variance: scale^2 * (gamma(1 + 2/shape) - ",
                      "gamma(1 + 1/shape)^2)"),
 constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = if ( .lss )
           c("scale", "shape") else
           c("shape", "scale"),
         lss = .lss ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .scale.12 = scale.12, .scale.TF = scale.TF,
           .lscale = lscale ,
           .lshape = lshape ,
           .lss = lss 
         ))),

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


    if (is.SurvS4(y))
      stop("only uncensored observations are allowed; ",
           "don't use SurvS4()")


    if ( .lss ) {
      mynames1 <- param.names("scale", ncoly)
      mynames2 <- param.names("shape", ncoly)
      predictors.names <-
          c(namesof(mynames1, .lscale , earg = .escale , tag = FALSE),
            namesof(mynames2, .lshape , earg = .eshape , tag = FALSE))
            
    } else {
      mynames1 <- param.names("shape", ncoly)
      mynames2 <- param.names("scale", ncoly)
      predictors.names <-
          c(namesof(mynames1, .lshape , earg = .eshape , tag = FALSE),
            namesof(mynames2, .lscale , earg = .escale , tag = FALSE))
    }
    predictors.names <- predictors.names[
          interleave.VGAM(M, M1 = M1)]


    Shape.init <- matrix(if (length( .ishape )) .ishape else 0 + NA,
                         n, ncoly, byrow = TRUE)
    Scale.init <- matrix(if (length( .iscale )) .iscale else 0 + NA,
                         n, ncoly, byrow = TRUE)

    if (!length(etastart)) {
      if (!length( .ishape ) ||
          !length( .iscale )) {
        for (ilocal in 1:ncoly) {

          anyc <- FALSE  # extra$leftcensored | extra$rightcensored
          i11 <- if ( .imethod == 1) anyc else FALSE  # Can be all data
          probs.y <- .probs.y
          xvec <- log(-log1p(-probs.y))
          fit0 <- lsfit(x  = xvec,
                        y  = log(quantile(y[!i11, ilocal],
                                 probs = probs.y )))


          if (!is.Numeric(Shape.init[, ilocal]))
            Shape.init[, ilocal] <- 1 / fit0$coef["X"]
          if (!is.Numeric(Scale.init[, ilocal]))
            Scale.init[, ilocal] <- exp(fit0$coef["Intercept"])
        }  # ilocal

        etastart <- if ( .lss )
          cbind(theta2eta(Scale.init, .lscale , earg = .escale ),
                theta2eta(Shape.init, .lshape , earg = .eshape ))[,
                interleave.VGAM(M, M1 = M1)] else
          cbind(theta2eta(Shape.init, .lshape , earg = .eshape ),
                theta2eta(Scale.init, .lscale , earg = .escale ))[,
                interleave.VGAM(M, M1 = M1)]
      }
    }
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .iscale = iscale, .ishape = ishape,
            .probs.y = probs.y,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss,
            .imethod = imethod ) )),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale )
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )
    Scale * gamma(1 + 1 / Shape)
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape,
           .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
  last = eval(substitute(expression({
    regnotok <- any(Shape <= 2)
    if (any(Shape <= 1)) {
      warning("MLE regularity conditions are violated",
              "(shape <= 1) at the final iteration: ",
              "MLEs are not consistent")
    } else if (any(1 < Shape & Shape < 2)) {
      warning("MLE regularity conditions are violated",
              "(1 < shape < 2) at the final iteration: ",
              "MLEs exist but are not asymptotically normal")
    } else if (any(2 == Shape)) {
      warning("MLE regularity conditions are violated",
              "(shape == 2) at the final iteration: ",
              "MLEs exist and are normal and asymptotically ",
              "efficient but with a slower convergence rate than when ",
              "shape > 2")
    }



    M1 <- extra$M1
    avector <- if ( .lss ) c(rep( .lscale , length = ncoly),
                             rep( .lshape , length = ncoly)) else
                           c(rep( .lshape , length = ncoly),
                             rep( .lscale , length = ncoly))
    misc$link <- avector[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- if ( .lss ) .escale else .eshape
      misc$earg[[M1*ii  ]] <- if ( .lss ) .eshape else .escale
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE


    misc$nrfs <- .nrfs
    misc$RegCondOK <- !regnotok # Save this for later
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .imethod = imethod,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss,
            .nrfs = nrfs ) )),
  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale )
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * dweibull(y, shape = Shape, scale = Scale, log = TRUE))
    }
  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape,
           .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
  vfamily = c("weibullR"),
  deriv = eval(substitute(expression({
    M1 <- 2
    Scale <- eta2theta(eta[,    .scale.TF  ], .lscale , earg = .escale )
    Shape <- eta2theta(eta[, !( .scale.TF )], .lshape , earg = .eshape )

    dl.dshape <- 1 / Shape + log(y / Scale) -
                 log(y / Scale) * (y / Scale)^Shape
    dl.dscale <- (Shape / Scale) * (-1.0 + (y / Scale)^Shape)

    dshape.deta <- dtheta.deta(Shape, .lshape, earg = .eshape )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )

    myderiv <- if ( .lss )
                 c(w) * cbind(dl.dscale * dscale.deta,
                              dl.dshape * dshape.deta) else
                 c(w) * cbind(dl.dshape * dshape.deta,
                              dl.dscale * dscale.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ) )),
  weight = eval(substitute(expression({
    EulerM <- -digamma(1.0)





    ned2l.dscale <- (Shape / Scale)^2
    ned2l.dshape <- (6*(EulerM - 1)^2 + pi^2)/(6*Shape^2)  # KK (2003)
    ned2l.dshapescale <- (EulerM-1) / Scale




    wz <- if ( .lss )
            array(c(c(w) * ned2l.dscale * dscale.deta^2,
                    c(w) * ned2l.dshape * dshape.deta^2,
                    c(w) * ned2l.dshapescale * dscale.deta * dshape.deta),
                  dim = c(n, M / M1, 3)) else
            array(c(c(w) * ned2l.dshape * dshape.deta^2,
                    c(w) * ned2l.dscale * dscale.deta^2,
                    c(w) * ned2l.dshapescale * dscale.deta * dshape.deta),
                  dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .eshape = eshape, .nrfs = nrfs,
            .scale.12 = scale.12, .scale.TF = scale.TF, .lss = lss ))))
}








setOldClass(c("SurvS4", "Surv"))


 SurvS4 <-
function (time, time2, event, type = c("right", "left", "interval",
    "counting", "interval2"), origin = 0) {
  nn <- length(time)
  ng <- nargs()
  if (missing(type)) {
    if (ng == 1 || ng == 2)
      type <- "right" else if (ng == 3)
      type <- "counting" else stop("Invalid number of arguments")
  } else {
    type <- match.arg(type)
    ng <- ng - 1
    if (ng != 3 && (type == "interval" || type == "counting"))
      stop("Wrong number of args for this type of survival data")
    if (ng != 2 && (type == "right" || type == "left" ||
        type == "interval2"))
      stop("Wrong number of args for this type of survival data")
  }
  who <- !is.na(time)
  if (ng == 1) {
    if (!is.numeric(time))
      stop("Time variable is not numeric")
    ss <- cbind(time, 1)
    dimnames(ss) <- list(NULL, c("time", "status"))
  } else if (type == "right" || type == "left") {
    if (!is.numeric(time))
      stop("Time variable is not numeric")
    if (length(time2) != nn)
      stop("Time and status are different lengths")
    if (is.logical(time2))
      status <- 1 * time2 else if (is.numeric(time2)) {
      who2 <- !is.na(time2)
      if (max(time2[who2]) == 2)
          status <- time2 - 1 else status <- time2
      if (any(status[who2] != 0 & status[who2] != 1))
        stop("Invalid status value")
    } else stop("Invalid status value")
    ss <- cbind(time, status)
    dimnames(ss) <- list(NULL, c("time", "status"))
  } else if (type == "counting") {
    if (length(time2) != nn)
      stop("Start and stop are different lengths")
    if (length(event) != nn)
      stop("Start and event are different lengths")
    if (!is.numeric(time))
      stop("Start time is not numeric")
    if (!is.numeric(time2))
      stop("Stop time is not numeric")
    who3 <- who & !is.na(time2)
    if (any(time[who3] >= time2[who3]))
      stop("Stop time must be > start time")
    if (is.logical(event))
      status <- 1 * event else if (is.numeric(event)) {
      who2 <- !is.na(event)
      if (max(event[who2]) == 2)
          status <- event - 1 else status <- event
      if (any(status[who2] != 0 & status[who2] != 1))
          stop("Invalid status value")
    } else stop("Invalid status value")
    ss <- cbind(time - origin, time2 - origin, status)
  } else {
    if (type == "interval2") {
      event <- ifelse(is.na(time), 2, ifelse(is.na(time2),
          0, ifelse(time == time2, 1, 3)))
      if (any(time[event == 3] > time2[event == 3]))
        stop("Invalid interval: start > stop")
      time <- ifelse(event != 2, time, time2)
      type <- "interval"
    } else {
      temp <- event[!is.na(event)]
      if (!is.numeric(temp))
        stop("Status indicator must be numeric")
      if (length(temp) > 0 && any(temp != floor(temp) |
          temp < 0 | temp > 3))
        stop("Status indicator must be 0, 1, 2 or 3")
    }
    status <- event
    ss <- cbind(time, ifelse(!is.na(event) & event == 3,
        time2, 1), status)
  }
  attr(ss, "type") <- type
  class(ss) <- "SurvS4"
  ss
}





is.SurvS4 <- function(x) inherits(x, "SurvS4")




 setIs(class1 = "SurvS4",
       class2 = "matrix")  # Forces vglm()@y to be a matrix



as.character.SurvS4 <- function (x, ...) {
  class(x) <- NULL
  type <- attr(x, "type")

  if (type == "right") {
    temp <- x[, 2]
    temp <- ifelse(is.na(temp), "?", ifelse(temp == 0, "+", " "))
    paste(format(x[, 1]), temp, sep = "")
  } else if (type == "counting") {
    temp <- x[, 3]
    temp <- ifelse(is.na(temp), "?", ifelse(temp == 0, "+", " "))
    paste("(", format(x[, 1]), ",", format(x[, 2]), temp, "]", sep = "")
  } else if (type == "left") {
    temp <- x[, 2]
    temp <- ifelse(is.na(temp), "?", ifelse(temp == 0, "<", " "))
    paste(temp, format(x[, 1]), sep = "")
  } else {
    stat <- x[, 3]
    temp <- c("+", "", "-", "]")[stat + 1]
    temp2 <- ifelse(stat == 3, paste("(", format(x[, 1]),
        ", ", format(x[, 2]), sep = ""), format(x[, 1]))
    ifelse(is.na(stat), as.character(NA), paste(temp2, temp, sep = ""))
  }
}



"[.SurvS4" <- function(x, i, j, drop = FALSE) {
    if (missing(j)) {
        temp <- class(x)
        type <- attr(x, "type")
        class(x) <- NULL
        x <- x[i, , drop = FALSE]
        class(x) <- temp
        attr(x, "type") <- type
        x
    } else {

        class(x) <- NULL
        NextMethod("[")
    }
}


is.na.SurvS4 <- function(x) {
  as.vector( (1* is.na(unclass(x)))%*% rep(1, ncol(x)) >0)
}







show.SurvS4 <- function (object)
  print.default(as.character.SurvS4(object), quote = FALSE)




setMethod("show", "SurvS4",
         function(object)
         show.SurvS4(object))








pgamma.deriv.unscaled <- function(q, shape) {




  gam0 <- exp(lgamma(shape) + pgamma(q = q, shape = shape, log.p = TRUE))

  I.sq <- pgamma(q = q, shape = shape)

  alld <- pgamma.deriv(q = q, shape = shape)  # 6-coln matrix
  tmp3 <- alld[, 3] / I.sq  # RHS of eqn (4.5) of \cite{wing:1989}
    
  G1s <- digamma(shape) + tmp3  # eqn (4.9)
  gam1 <- gam0 * G1s


  dG1s <- trigamma(shape) + alld[, 4] / I.sq - tmp3^2  # eqn (4.13)
  G2s <- dG1s + G1s^2  # eqn (4.12)

  gam2 <- gam0 * G2s

  cbind("0" = gam0,
        "1" = gam1,
        "2" = gam2)
}






 truncweibull <-
  function(lower.limit = 1e-5,
           lAlpha = "loge", lBetaa = "loge",
           iAlpha = NULL,   iBetaa = NULL,
           nrfs = 1,
           probs.y = c(0.2, 0.5, 0.8),
           imethod = 1,
           zero = "Betaa") {









  lAlpha <- as.list(substitute(lAlpha))
  eAlpha <- link2list(lAlpha)
  lAlpha <- attr(eAlpha, "function.name")

  lBetaa <- as.list(substitute(lBetaa))
  eBetaa <- link2list(lBetaa)
  lBetaa <- attr(eBetaa, "function.name")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")

  if (!is.Numeric(probs.y, positive  = TRUE) ||
      length(probs.y) < 2 ||
      max(probs.y) >= 1)
    stop("bad input for argument 'probs.y'")


  if (!is.Numeric(nrfs, length.arg = 1) ||
      nrfs < 0 ||
      nrfs > 1)
    stop("bad input for argument 'nrfs'")


  if (length(iAlpha))
    if (!is.Numeric(iAlpha, positive = TRUE))
      stop("argument 'iAlpha' values must be positive")

  if (length(iBetaa))
    if (!is.Numeric(iBetaa, positive = TRUE))
      stop("argument 'iBetaa' values must be positive")


  new("vglmff",
  blurb = c("Truncated weibull distribution\n\n",
            "Links:    ",
            namesof("Alpha", lAlpha, earg = eAlpha), ", ", 
            namesof("Betaa", lBetaa, earg = eBetaa), "\n", 
            if (length( lower.limit ) < 5)
              paste("Truncation point(s):     ",
                    lower.limit, sep = ", ") else
              ""),
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
         parameters.names = c("Alpha", "Betaa"),
         lower.limit = .lower.limit ,
         lAlpha = .lAlpha ,
         lBetaa = .lBetaa ,
         zero = .zero )
  }, list( .zero = zero,
           .lAlpha = lAlpha ,
           .lBetaa = lBetaa ,
           .lower.limit = lower.limit
         ))),

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

    extra$lower.limit <- matrix( .lower.limit , n, ncoly, byrow = TRUE)

    if (any(y < extra$lower.limit)) {
      stop("some response values less than argument 'lower.limit'")
    }



    if (is.SurvS4(y))
      stop("only uncensored observations are allowed; ",
           "don't use SurvS4()")


    mynames1 <- param.names("Alpha", ncoly)
    mynames2 <- param.names("Betaa", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lAlpha , earg = .eAlpha , tag = FALSE),
          namesof(mynames2, .lBetaa , earg = .eBetaa , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


    Alpha.init <- matrix(if (length( .iAlpha )) .iAlpha else 0 + NA,
                         n, ncoly, byrow = TRUE)
    Betaa.init <- matrix(if (length( .iBetaa )) .iBetaa else 0 + NA,
                         n, ncoly, byrow = TRUE)

    if (!length(etastart)) {
      if (!length( .iAlpha ) ||
          !length( .iBetaa )) {
        for (ilocal in 1:ncoly) {

          anyc <- FALSE  # extra$leftcensored | extra$rightcensored
          i11 <- if ( .imethod == 1) anyc else FALSE  # Can be all data
          probs.y <- .probs.y
          xvec <- log(-log1p(-probs.y))
          fit0 <- lsfit(x  = xvec,
                        y  = log(quantile(y[!i11, ilocal],
                                 probs = probs.y )))
          aaa.init <- 1 / fit0$coef["X"]
          bbb.init <- exp(fit0$coef["Intercept"])

          if (!is.Numeric(Betaa.init[, ilocal]))
            Betaa.init[, ilocal] <- aaa.init
          if (!is.Numeric(Alpha.init[, ilocal]))
            Alpha.init[, ilocal] <- (1 / bbb.init)^aaa.init
        }  # ilocal
      } else {
        Alpha.init <- rep( .iAlpha , length = n)
        Betaa.init <- rep( .iBetaa , length = n)
      }

      etastart <-
        cbind(theta2eta(Alpha.init, .lAlpha , earg = .eAlpha ),
              theta2eta(Betaa.init, .lBetaa , earg = .eBetaa ))[,
              interleave.VGAM(M, M1 = M1)]
    }
  }), list( .lBetaa = lBetaa, .lAlpha = lAlpha,
            .eBetaa = eBetaa, .eAlpha = eAlpha,
            .iBetaa = iBetaa, .iAlpha = iAlpha,
            .lower.limit = lower.limit,
            .probs.y = probs.y,
            .imethod = imethod ) )),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    Alpha <- eta2theta(eta[, c(TRUE, FALSE)], .lAlpha , earg = .eAlpha )
    Betaa <- eta2theta(eta[, c(FALSE, TRUE)], .lBetaa , earg = .eBetaa )


    aTb <- Alpha * extra$lower.limit^Betaa
    wingo3 <- pgamma.deriv.unscaled(q = aTb, shape = 1 + 1 / Betaa)
    exp.aTb <- exp(aTb)

    (gamma(1 + 1 / Betaa) - wingo3[, 1]) *
    exp.aTb / Alpha^(1 / Betaa)
  }, list( .lBetaa = lBetaa, .lAlpha = lAlpha,
           .eBetaa = eBetaa, .eAlpha = eAlpha,
           .lower.limit = lower.limit) )),
  last = eval(substitute(expression({

    aaa.hat <- Betaa
    regnotok <- any(aaa.hat <= 2)
    if (any(aaa.hat <= 1)) {
      warning("MLE regularity conditions are violated",
              "(Betaa <= 1) at the final iteration: ",
              "MLEs are not consistent")
    } else if (any(1 < aaa.hat & aaa.hat < 2)) {
      warning("MLE regularity conditions are violated",
              "(1 < Betaa < 2) at the final iteration: ",
              "MLEs exist but are not asymptotically normal")
    } else if (any(2 == aaa.hat)) {
      warning("MLE regularity conditions are violated",
              "(Betaa == 2) at the final iteration: ",
              "MLEs exist and are normal and asymptotically ",
              "efficient but with a slower convergence rate than when ",
              "Betaa > 2")
    }



    M1 <- extra$M1
    misc$link <-
      c(rep( .lAlpha , length = ncoly),
        rep( .lBetaa , length = ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .eAlpha
      misc$earg[[M1*ii  ]] <- .eBetaa
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE


    misc$nrfs <- .nrfs
    misc$RegCondOK <- !regnotok # Save this for later



  }), list( .lBetaa = lBetaa, .lAlpha = lAlpha,
            .eBetaa = eBetaa, .eAlpha = eAlpha,
            .imethod = imethod,
            .lower.limit = lower.limit,
            .nrfs = nrfs ) )),

  loglikelihood = eval(substitute(
          function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    Alpha <- eta2theta(eta[, c(TRUE, FALSE)], .lAlpha , earg = .eAlpha )
    Betaa <- eta2theta(eta[, c(FALSE, TRUE)], .lBetaa , earg = .eBetaa )
    Shape <- Betaa
    Scale <- 1 / Alpha^(1/Betaa)

    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(c(w) * (dweibull(x = y, shape = Shape,
                           scale = Scale, log = TRUE) -
                  pweibull(q = extra$lower.limit, shape = Shape,
                           scale = Scale, log.p = TRUE,
                           lower.tail = FALSE)))
    }
  }, list( .lBetaa = lBetaa, .lAlpha = lAlpha,
           .eBetaa = eBetaa, .eAlpha = eAlpha,
           .lower.limit = lower.limit ) )),

  vfamily = c("truncweibull"),

  deriv = eval(substitute(expression({
    M1 <- 2
    Alpha <- eta2theta(eta[, c(TRUE, FALSE)], .lAlpha , earg = .eAlpha )
    Betaa <- eta2theta(eta[, c(FALSE, TRUE)], .lBetaa , earg = .eBetaa )

    Shape <- Betaa
    Scale <- 1 / Alpha^(1/Betaa)
    TTT <- extra$lower.limit
    dl.dAlpha <- 1 / Alpha - y^Betaa + TTT^Betaa
    dl.dBetaa <- (1 / Betaa) + log(y) -
                 Alpha * (y^Betaa * log(y) -
                          TTT^Betaa * log(TTT))

    dAlpha.deta <- dtheta.deta(Alpha, .lAlpha, earg = .eAlpha )
    dBetaa.deta <- dtheta.deta(Betaa, .lBetaa, earg = .eBetaa )

    myderiv <- c(w) * cbind(dl.dAlpha * dAlpha.deta,
                            dl.dBetaa * dBetaa.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lBetaa = lBetaa, .lAlpha = lAlpha,
            .eBetaa = eBetaa, .eAlpha = eAlpha,
            .lower.limit = lower.limit ) )),


  weight = eval(substitute(expression({
    aTb <- Alpha * TTT^Betaa
    exp.aTb <- exp(aTb)
    TblogT <- (TTT^Betaa) * log(TTT)
    wingo3 <- pgamma.deriv.unscaled(q = aTb,
                                    shape = 2)  # 3-cols

 
    Eyblogy <- (exp.aTb * (digamma(2) - wingo3[, 2]) -
               (aTb + 1) * log(Alpha)) / (Alpha * Betaa)


 
    Eyblog2y <- (exp.aTb * (digamma(2)^2 + trigamma(2) -
                 wingo3[, 3]) - 2 * log(Alpha) *
                (digamma(2) - wingo3[, 2])) / (Alpha * Betaa^2) +
                (log(Alpha)^2) * (aTb + 1) / (Alpha * Betaa^2)

    ned2l.daa <- 1 / Alpha^2
    ned2l.dab <- Eyblogy - TblogT
    ned2l.dbb <- (1 / Betaa)^2 + Alpha * Eyblog2y -
                 aTb * (log(TTT))^2





    wz <- array(c(c(w) * ned2l.daa * dAlpha.deta^2,
                  c(w) * ned2l.dbb * dBetaa.deta^2,
                  c(w) * ned2l.dab * dBetaa.deta * dAlpha.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)
    wz
  }), list( .nrfs = nrfs ))))
}





