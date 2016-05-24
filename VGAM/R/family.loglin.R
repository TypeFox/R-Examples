# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.





 loglinb2 <- function(exchangeable = FALSE, zero = "u12") {


  if (!is.logical(exchangeable))
    warning("argument 'exchangeable' should be a single logical")


  new("vglmff",
  blurb = c("Log-linear model for binary data\n\n",
            "Links:    ",
            "Identity: u1, u2, u12",
            "\n"),
  constraints = eval(substitute(expression({
    cm.intercept.default <- diag(3)

    constraints <- cm.VGAM(matrix(c(1,1,0, 0,0,1), 3, 2), x = x,
                           bool = .exchangeable ,
                           constraints = constraints,
                           apply.int = TRUE,
                           cm.default           = cm.intercept.default,
                           cm.intercept.default = cm.intercept.default)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .exchangeable = exchangeable, .zero = zero ))),



  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 4,  # ncol(fitted(object))
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("u1", "u2", "u12"),
         zero = .zero )
  }, list( .zero = zero
         ))),


  initialize = expression({


    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 2,
              out.wy = TRUE,
              colsyperw = 2,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (ncol(y) != 2)
      stop("ncol(y) must be = 2")

    predictors.names <- c("u1", "u2", "u12")

    if (length(mustart) + length(etastart) == 0) {
      mustart <- matrix(NA_real_, nrow(y), 4)
      mustart[,1] <- weighted.mean((1-y[,1])*(1-y[,2]), w)
      mustart[,2] <- weighted.mean((1-y[,1])*y[,2], w)
      mustart[,3] <- weighted.mean(y[,1]*(1-y[,2]), w)
      mustart[,4] <- weighted.mean(y[,1]*y[,2], w)
      if (any(mustart == 0))
        stop("some combinations of the response not realized") 
    }
  }),
  linkinv = function(eta, extra = NULL) {
    u1 <-  eta[,1]
    u2 <-  eta[,2]
    u12 <- eta[,3]
    denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
    cbind("00" = 1/denom,
          "01" = exp(u2) / denom,
          "10" = exp(u1) / denom,
          "11" = exp(u1+u2+u12) / denom)
  },
  last = expression({
    misc$link <-    c("u1" = "identitylink", "u2" = "identitylink",
                      "u12" = "identitylink")
    misc$earg <- list("u1"  = list(),    "u2"  = list(),
                      "u12"  = list())

    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }),
  linkfun = function(mu, extra = NULL)  {
    u0 <-  log(mu[,1]) 
    u2 <-  log(mu[,2]) - u0
    u1 <-  log(mu[,3]) - u0
    u12 <- log(mu[,4]) - u0 - u1 - u2 
    cbind(u1, u2, u12)
  },
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    u1 <-  eta[,1]
    u2 <-  eta[,2]
    u12 <- eta[,3]
    denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
    u0 <- -log(denom)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * (u0 + u1*y[,1] + u2*y[,2] + u12*y[,1]*y[,2])
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("loglinb2"),
  deriv = expression({
    u1 <-  eta[,1]
    u2 <-  eta[,2]
    u12 <- eta[,3]
    denom <- 1 + exp(u1) + exp(u2) + exp(u1 + u2 + u12)
    du0.du1 <- -(exp(u1) + exp(u1 + u2 + u12)) / denom 
    du0.du2 <- -(exp(u2) + exp(u1 + u2 + u12)) / denom 
    du0.du12 <- -exp(u1 + u2 + u12) / denom 
    c(w) * cbind(du0.du1  + y[,1], 
                 du0.du2  + y[,2],
                 du0.du12 + y[,1] * y[,2]) 
  }),
  weight = expression({
    d2u0.du1.2 <- -(exp(u1) + exp(u1 + u2 + u12)) * (1+exp(u2)) / denom^2 
    d2u0.du22 <-  -(exp(u2) + exp(u1 + u2 + u12)) * (1+exp(u1)) / denom^2 
    d2u0.du122 <- -exp(u1 + u2 + u12) * (1+exp(u1)+exp(u2)) / denom^2 
    d2u0.du1u2 <- -(exp(u1 + u2 + u12) - exp(u1 + u2)) / denom^2 
    d2u0.du1u3 <- -(1 + exp(u2)) * exp(u1 + u2 + u12) / denom^2 
    d2u0.du2u3 <- -(1 + exp(u1)) * exp(u1 + u2 + u12) / denom^2 

    wz <- matrix(NA_real_, n, dimm(M)) 
    wz[,iam(1,1,M)] <- -d2u0.du1.2 
    wz[,iam(2,2,M)] <- -d2u0.du22
    wz[,iam(3,3,M)] <- -d2u0.du122 
    wz[,iam(1,2,M)] <- -d2u0.du1u2
    wz[,iam(1,3,M)] <- -d2u0.du1u3
    wz[,iam(2,3,M)] <- -d2u0.du2u3
    c(w) * wz
  }))
}




 loglinb3 <- function(exchangeable = FALSE,
                      zero = c("u12", "u13", "u23")) {


  if (!is.logical(exchangeable))
    warning("argument 'exchangeable' should be a single logical")


  new("vglmff",
  blurb = c("Log-linear model for trivariate binary data\n\n",
            "Links:    ",
            "Identity: u1, u2, u3, u12, u13, u23",
            "\n"),
  constraints = eval(substitute(expression({
    cm.intercept.default <- diag(6)

    constraints <- cm.VGAM(matrix(c(1,1,1,0,0,0, 0,0,0,1,1,1), 6, 2), x = x,
                           bool = .exchangeable ,
                           constraints = constraints,
                           apply.int = TRUE,
                           cm.default           = cm.intercept.default,
                           cm.intercept.default = cm.intercept.default)
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 6)
  }), list( .exchangeable = exchangeable, .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 6,
         Q1 = 8,  # ncol(fitted(object))
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("u1", "u2", "u3", "u12", "u13", "u23"),
         zero = .zero )
  }, list( .zero = zero
         ))),


  initialize = expression({
    predictors.names <- c("u1", "u2", "u3", "u12", "u13", "u23")


    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.max = 3,
              out.wy = TRUE,
              colsyperw = 3,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    if (ncol(y) != 3)
      stop("ncol(y) must be = 3")


    if (FALSE)
    extra$my.expression <- expression({
      u1  <- eta[, 1]
      u2  <- eta[, 2]
      u3  <- eta[, 3]
      u12 <- eta[, 4]
      u13 <- eta[, 5]
      u23 <- eta[, 6]
      denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
               exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
               exp(u1 + u2 + u3 + u12 + u13 + u23)
    })



    if (length(mustart) + length(etastart) == 0) {
      mustart <- matrix(NA_real_, nrow(y), 2^3)
      mustart[,1] <- weighted.mean((1-y[,1])*(1-y[,2])*(1-y[,3]), w)
      mustart[,2] <- weighted.mean((1-y[,1])*(1-y[,2])*y[,3], w)
      mustart[,3] <- weighted.mean((1-y[,1])*y[,2]*(1-y[,3]), w)
      mustart[,4] <- weighted.mean((1-y[,1])*y[,2]*y[,3], w)
      mustart[,5] <- weighted.mean(y[,1]*(1-y[,2])*(1-y[,3]), w)
      mustart[,6] <- weighted.mean(y[,1]*(1-y[,2])*y[,3], w)
      mustart[,7] <- weighted.mean(y[,1]*y[,2]*(1-y[,3]), w)
      mustart[,8] <- weighted.mean(y[,1]*y[,2]*y[,3], w)
      if (any(mustart == 0))
        stop("some combinations of the response not realized") 
    }
  }),
  linkinv = function(eta, extra = NULL) {
      u1  <- eta[, 1]
      u2  <- eta[, 2]
      u3  <- eta[, 3]
      u12 <- eta[, 4]
      u13 <- eta[, 5]
      u23 <- eta[, 6]
      denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
               exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
               exp(u1 + u2 + u3 + u12 + u13 + u23)


    cbind("000" = 1,
          "001" = exp(u3),
          "010" = exp(u2),
          "011" = exp(u2+u3+u23),
          "100" = exp(u1),
          "101" = exp(u1+u3+u13),
          "110" = exp(u1+u2+u12),
          "111" = exp(u1+u2+u3+u12+u13+u23)) / denom
  },
  last = expression({
    misc$link <- rep("identitylink", length = M)
    names(misc$link) <- predictors.names

    misc$earg <- list(u1  = list(), u2  = list(), u3  = list(),
                      u12 = list(), u13 = list(), u23 = list())

    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

  }),
  linkfun = function(mu, extra = NULL)  {
    u0  <- log(mu[,1])
    u3  <- log(mu[,2]) - u0
    u2  <- log(mu[,3]) - u0
    u23 <- log(mu[,4]) - u0 - u2 - u3
    u1  <- log(mu[,5]) - u0
    u13 <- log(mu[,6]) - u0 - u1 - u3
    u12 <- log(mu[,7]) - u0 - u1 - u2
    cbind(u1, u2, u3, u12, u13, u23)
  },
  loglikelihood =
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    u1  <- eta[, 1]
    u2  <- eta[, 2]
    u3  <- eta[, 3]
    u12 <- eta[, 4]
    u13 <- eta[, 5]
    u23 <- eta[, 6]
    denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
             exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
             exp(u1 + u2 + u3 + u12 + u13 + u23)

    u0 <- -log(denom)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * (u0 + u1*y[,1] + u2*y[,2] + u3*y[,3] +u12*y[,1]*y[,2] +
                u13*y[,1]*y[,3] + u23*y[,2]*y[,3])
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  },
  vfamily = c("loglinb3"),
  deriv = expression({
    u1  <- eta[, 1]
    u2  <- eta[, 2]
    u3  <- eta[, 3]
    u12 <- eta[, 4]
    u13 <- eta[, 5]
    u23 <- eta[, 6]
    denom <- 1 + exp(u1) + exp(u2) + exp(u3) + exp(u1 + u2 + u12) +
             exp(u1 + u3 + u13) + exp(u2 + u3 + u23) +
             exp(u1 + u2 + u3 + u12 + u13 + u23)



    allterms <- exp(u1+u2+u3+u12+u13+u23)
    A1 <- exp(u1) + exp(u1 + u2 + u12) + exp(u1 + u3 + u13) +
          allterms
    A2 <- exp(u2) + exp(u1 + u2 + u12) + exp(u2 + u3 + u23) +
          allterms
    A3 <- exp(u3) + exp(u3 + u2 + u23) + exp(u1 + u3 + u13) +
          allterms
    A12 <- exp(u1 + u2 + u12) + allterms
    A13 <- exp(u1 + u3 + u13) + allterms
    A23 <- exp(u2 + u3 + u23) + allterms


    c(w) * cbind(-A1/denom + y[,1], 
                 -A2/denom + y[,2],
                 -A3/denom + y[,3],
                 -A12/denom + y[,1]*y[,2],
                 -A13/denom + y[,1]*y[,3],
                 -A23/denom + y[,2]*y[,3])
  }),
  weight = expression({
    u0 <- -log(denom)
    dA2.du1 <- exp(u1 + u2 + u12) + allterms
    dA3.du1 <- exp(u1 + u3 + u13) + allterms
    dA3.du2 <- exp(u2 + u3 + u23) + allterms

    wz <- matrix(NA_real_, n, dimm(6)) 
    expu0 <- exp(u0)

    wz[,iam(1,1,M)] <- A1 * (1 - expu0 * A1)
    wz[,iam(2,2,M)] <- A2 * (1 - expu0 * A2)
    wz[,iam(3,3,M)] <- A3 * (1 - expu0 * A3)
    wz[,iam(1,2,M)] <- (dA2.du1 - expu0 * A1 * A2)
    wz[,iam(1,3,M)] <- (dA3.du1 - expu0 * A1 * A3)
    wz[,iam(2,3,M)] <- (dA3.du2 - expu0 * A2 * A3)
    wz[,iam(4,4,M)] <- A12 * (1 - expu0 * A12)
    wz[,iam(5,5,M)] <- A13 * (1 - expu0 * A13)
    wz[,iam(6,6,M)] <- A23 * (1 - expu0 * A23)
    wz[,iam(4,6,M)] <- (allterms - expu0 * A12 * A23)
    wz[,iam(5,6,M)] <- (allterms - expu0 * A12 * A23)
    wz[,iam(4,5,M)] <- (allterms - expu0 * A12 * A13)
    wz[,iam(1,4,M)] <- A12 * (1 - expu0 * A1)
    wz[,iam(1,5,M)] <- A13 * (1 - expu0 * A1)
    wz[,iam(1,6,M)] <- (allterms - expu0 * A1 * A23)
    wz[,iam(2,4,M)] <- A12 * (1 - expu0 * A2)
    wz[,iam(2,5,M)] <- (allterms - expu0 * A2 * A13)
    wz[,iam(2,6,M)] <- A23 * (1 - expu0 * A2)
    wz[,iam(3,4,M)] <- (allterms - expu0 * A3 * A12)
    wz[,iam(3,5,M)] <- A13 * (1 - expu0 * A3)
    wz[,iam(3,6,M)] <- A23 * (1 - expu0 * A3)
    wz <- expu0 * wz 
    c(w) * wz
  }))
}


