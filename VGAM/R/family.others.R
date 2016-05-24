# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.


















dexppois <- function(x, rate = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(shape), length(rate))
  x     <- rep(x,     len = N)
  shape <- rep(shape, len = N)
  rate  <- rep(rate,  len = N)

  logdensity <- rep(log(0), len = N)
  xok <- (0 < x)
 
  logdensity[xok] <- log(shape[xok]) + log(rate[xok]) -
                     log1p(-exp(-shape[xok])) - shape[xok] - 
                     rate[xok] * x[xok] + shape[xok] * 
                     exp(-rate[xok] * x[xok])
   
  logdensity[shape <= 0] <- NaN
  logdensity[rate <= 0] <- NaN
  if (log.arg) logdensity else exp(logdensity)
}






qexppois<- function(p, rate = 1, shape, 
                    lower.tail = TRUE, log.p = FALSE) { 
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- -log(log(exp(ln.p) * (-expm1(shape)) + exp(shape)) / shape) / rate
      ans[ln.p > 0] <- NaN
    } else {
      ans <- -log(log(p * (-expm1(shape)) + exp(shape)) / shape) / rate
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- -log(log(expm1(ln.p) * expm1(shape) + exp(shape)) / shape) / rate
      ans[ln.p > 0] <- NaN
    } else { 
      ans <- -log(log(p * expm1(shape) + 1) / shape) / rate
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  }
  ans[(shape <= 0) | (rate <= 0)] <- NaN
  ans
}














pexppois<- function(q, rate = 1, shape, 
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ans <- log((exp(shape * exp(-rate * q)) -
                    exp(shape)) / -expm1(shape))
      ans[q <= 0 ] <- -Inf
      ans[q == Inf] <- 0
    } else {
      ans <- (exp(shape * exp(-rate * q)) - exp(shape)) / (-expm1(shape))
      ans[q <= 0] <- 0
      ans[q == Inf] <- 1
    }
  } else {
    if (log.p) {
      ans <- log(expm1(shape * exp(-rate * q)) / expm1(shape))
      ans[q <= 0] <- 0
      ans[q == Inf] <- -Inf
    } else {
      ans <- expm1(shape * exp(-rate * q)) / expm1(shape)
      ans[q <= 0] <- 1
      ans[q == Inf] <- 0
    }
  } 
  ans[(shape <= 0) | (rate <= 0)] <- NaN
  ans
}




rexppois <- function(n, rate = 1, shape) {
  ans <- -log(log(runif(n) * (-expm1(shape)) +
         exp(shape)) / shape) / rate
  ans[(shape <= 0) | (rate <= 0)] <- NaN
  ans
}









 exppoisson <- function(lrate = "loge", lshape = "loge",
                        irate = 2.0, ishape = 1.1,   
                        zero = NULL) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lratee <- as.list(substitute(lrate))
  eratee <- link2list(lratee)
  lratee <- attr(eratee, "function.name")


  iratee <- irate




  if (length(ishape) &&
      !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (length(iratee) &&
      !is.Numeric(iratee, positive = TRUE))
    stop("bad input for argument 'irate'")

  ishape[abs(ishape - 1) < 0.01] = 1.1


  new("vglmff",
  blurb = c("Exponential Poisson distribution \n \n",
            "Links:    ",
            namesof("rate",  lratee, earg = eratee), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     shape/(expm1(shape) * rate)) * ",
                      "genhypergeo(c(1, 1), c(2, 2), shape)"),


  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("rate", "shape"),
         lrate  = .lratee ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lratee = lratee, .lshape = lshape ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <- c(
      namesof("rate",  .lratee , earg = .eratee , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {
      ratee.init <- if (length( .iratee ))
              rep( .iratee , len = n) else
              stop("Need to input a value into argument 'iratee'")
      shape.init <- if (length( .ishape ))
                      rep( .ishape , len = n) else
                      (1/ratee.init - mean(y)) / ((y * 
                      exp(-ratee.init * y))/n)


      ratee.init <- rep(weighted.mean(ratee.init, w = w), len = n)
      
      etastart <-
        cbind(theta2eta(ratee.init, .lratee , earg = .eratee ),
              theta2eta(shape.init, .lshape , earg = .eshape ))
              

    }
  }), list( .lshape = lshape, .lratee = lratee, 
            .ishape = ishape, .iratee = iratee, 
            .eshape = eshape, .eratee = eratee))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )





    qexppois(p = 0.5, rate = ratee, shape = shape)
  }, list( .lshape = lshape, .lratee = lratee, 
           .eshape = eshape, .eratee = eratee))), 

  last = eval(substitute(expression({
    misc$link <-    c( rate = .lratee , shape = .lshape )
    misc$earg <- list( rate = .eratee , shape = .eshape )

    misc$expected <- TRUE
    misc$multipleResponses <- FALSE
  }), list( .lshape = lshape, .lratee = lratee,
            .eshape = eshape, .eratee = eratee))), 

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dexppois(x = y, shape = shape, rate = ratee,
                                 log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lratee = lratee , .lshape = lshape , 
           .eshape = eshape , .eratee = eratee ))), 

  vfamily = c("exppoisson"),

  deriv = eval(substitute(expression({
    ratee <- eta2theta(eta[, 1], .lratee , earg = .eratee )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    dl.dratee <- 1/ratee - y - y * shape * exp(-ratee * y)
    dl.dshape <- 1/shape - 1/expm1(shape) - 1 + exp(-ratee * y)
    dratee.deta <- dtheta.deta(ratee, .lratee , earg = .eratee )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    c(w) * cbind(dl.dratee * dratee.deta,
                 dl.dshape * dshape.deta)
  }), list( .lshape = lshape, .lratee = lratee,
            .eshape = eshape, .eratee = eratee ))), 

  weight = eval(substitute(expression({
    
    temp1 <- -expm1(-shape)
    
    ned2l.dshape2 <- (1 + exp(2 * shape) - shape^2 * exp(shape) - 2 *
                      exp(shape)) / (shape * temp1)^2


    ned2l.dratee2 <- 1 / ratee^2 - (shape^2 * exp(-shape) / (4 * 
                    ratee^2 * temp1)) * 
                    genhypergeo(c(2, 2, 2), c(3, 3, 3), shape) 

    ned2l.drateeshape <- (shape * exp(-shape) / (4 * ratee * temp1)) *
                           genhypergeo(c(2, 2), c(3, 3), shape)   

    wz <- matrix(0, n, dimm(M))
    wz[, iam(1, 1, M)] <- dratee.deta^2 * ned2l.dratee2
    wz[, iam(1, 2, M)] <- dratee.deta * dshape.deta * ned2l.drateeshape
    wz[, iam(2, 2, M)] <- dshape.deta^2 * ned2l.dshape2
    c(w) * wz
  }), list( .zero = zero ))))
}










dgenray <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)



  N <- max(length(x), length(shape), length(scale))
  x <- rep(x, len = N)
  shape <- rep(shape, len = N)
  scale <- rep(scale, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- x[xok] / scale[xok]
    logdensity[xok] <- log(2) + log(shape[xok]) + log(x[xok]) -
                       2 * log(scale[xok]) - temp1^2  +
                       (shape[xok] - 1) * log1p(-exp(-temp1^2))
  }
  logdensity[(shape <= 0) | (scale <= 0)] <- NaN
  logdensity[is.infinite(x)] <- log(0)  # 20141209 KaiH
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}






pgenray <- function(q, scale = 1, shape, 
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ans <- log((-expm1(-(q/scale)^2))^shape)
      ans[q <= 0 ] <- -Inf
    } else {
      ans <- (-expm1(-(q/scale)^2))^shape
      ans[q <= 0] <- 0
    }
  } else {
    if (log.p) {
      ans <- log(-expm1(shape*log(-expm1(-(q/scale)^2))))
      ans[q <= 0] <- 0
    } else {
      ans <- -expm1(shape*log(-expm1(-(q/scale)^2)))
      ans[q <= 0] <- 1
    }
  } 
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}









qgenray <- function(p, scale = 1, shape, 
                    lower.tail = TRUE, log.p = FALSE) {
  if (!is.logical(lower.tail) || length(lower.tail ) != 1)
    stop("bad input for argument 'lower.tail'")
  
  if (!is.logical(log.p) || length(log.p) != 1)
    stop("bad input for argument 'log.p'")
  
  if (lower.tail) {
    if (log.p) {
      ln.p <- p
      ans <- scale * sqrt(-log1p(-(exp(ln.p)^(1/shape))))
      ans[ln.p > 0] <- NaN
    } else {
      ans <- scale * sqrt(-log1p(-(p^(1/shape))))
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  } else {
    if (log.p) {
      ln.p <- p
      ans <- scale * sqrt(-log1p(-((-expm1(ln.p))^(1/shape))))
      ans[ln.p > 0] <- NaN
    } else { 
      ans <- scale * sqrt(-log1p(-exp((1/shape)*log1p(-p))))
      ans[p < 0] <- NaN
      ans[p > 1] <- NaN
    }
  }
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}






rgenray <- function(n, scale = 1, shape) {
  ans <- qgenray(runif(n), shape = shape, scale = scale)
  ans[(shape <= 0) | (scale <= 0)] <- NaN
  ans
}




genrayleigh.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}


 genrayleigh <-
  function(lscale = "loge", lshape = "loge",
           iscale = NULL,   ishape = NULL,
           tol12 = 1.0e-05, 
           nsimEIM = 300, zero = 2) {

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  if (length(ishape) &&
      !is.Numeric(ishape, positive = TRUE))
    stop("bad input for argument 'ishape'")
  if (length(iscale) &&
      !is.Numeric(iscale, positive = TRUE)) 
    stop("bad input for argument 'iscale'")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE) ||
      nsimEIM <= 50)
      stop("argument 'nsimEIM' should be an integer greater than 50")



  new("vglmff",
  blurb = c("Generalized Rayleigh distribution\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n"),
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
         nsimEIM = .nsimEIM ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lscale = lscale, .lshape = lshape,
           .nsimEIM = nsimEIM ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y





    predictors.names <- c(
      namesof("scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {
      genrayleigh.Loglikfun <- function(scale, y, x, w, extraargs) {
        temp1 <- y / scale
        shape <- -1 / weighted.mean(log1p(-exp(-temp1^2)), w = w)

        ans <- sum(c(w) * (log(2) + log(shape) + log(y) -
                           2 * log(scale) - temp1^2  +
                           (shape - 1) * log1p(-exp(-temp1^2))))
        ans
      }
      scale.grid <- seq(0.2 * stats::sd(c(y)),
                        5.0 * stats::sd(c(y)), len = 29)
      scale.init <- if (length( .iscale )) .iscale else
                    grid.search(scale.grid, objfun = genrayleigh.Loglikfun,
                                y = y, x = x, w = w)
      scale.init <- rep(scale.init, length = length(y))
 
      shape.init <- if (length( .ishape )) .ishape else
                    -1 / weighted.mean(log1p(-exp(-(y/scale.init)^2)),
                     w = w)
      shape.init <- rep(shape.init, length = length(y))

      etastart <- cbind(theta2eta(scale.init, .lscale , earg = .escale ),
                        theta2eta(shape.init, .lshape , earg = .eshape ))
                        
        }
    }), list( .lscale = lscale, .lshape = lshape,
              .iscale = iscale, .ishape = ishape,
              .escale = escale, .eshape = eshape))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    qgenray(p = 0.5, shape = shape, scale = Scale)
  }, list( .lshape = lshape, .lscale = lscale, 
           .eshape = eshape, .escale = escale ))),

  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale, shape = .lshape )

    misc$earg <- list(scale = .escale, shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lshape = lshape, .lscale = lscale,
            .eshape = eshape, .escale = escale,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dgenray(x = y, shape = shape,
                                scale = Scale, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lshape = lshape , .lscale = lscale , 
           .eshape = eshape , .escale = escale ))), 
      
  vfamily = c("genrayleigh"),

  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas <- cbind(dscale.deta, dshape.deta)

    temp1 <- y / Scale
    temp2 <- exp(-temp1^2)
    temp3 <- temp1^2 / Scale
    AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
    BBB   <- -expm1(-temp1^2)     # denominator
    dl.dshape <- 1/shape + log1p(-temp2)
    dl.dscale <- -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

    dl.dshape[!is.finite(dl.dshape)] =
      max(dl.dshape[is.finite(dl.dshape)])

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale ))),

  weight = eval(substitute(expression({


    run.varcov <- 0
    ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)
    for (ii in 1:( .nsimEIM )) {
      ysim <- rgenray(n = n, shape = shape, scale = Scale)

      temp1 <- ysim / Scale
      temp2 <- exp(-temp1^2)  # May be 1 if ysim is very close to 0.
      temp3 <- temp1^2 / Scale
      AAA   <- 2 * temp1^2 / Scale  # 2 * y^2 / Scale^3
      BBB   <- -expm1(-temp1^2)     # denominator
      dl.dshape <- 1/shape + log1p(-temp2)
      dl.dscale <- -2 / Scale + AAA * (1 - (shape - 1) * temp2 / BBB)

      dl.dshape[!is.finite(dl.dshape)] <- max(
      dl.dshape[is.finite(dl.dshape)])

      temp3 <- cbind(dl.dscale, dl.dshape)
      run.varcov <- run.varcov + temp3[, ind1$row.index] *
                                 temp3[, ind1$col.index]
    }
    run.varcov <- run.varcov / .nsimEIM

    wz <- if (intercept.only)
        matrix(colMeans(run.varcov, na.rm = FALSE),
               n, ncol(run.varcov), byrow = TRUE) else run.varcov
    wz <- wz * dthetas.detas[, ind1$row] * dthetas.detas[, ind1$col]
    c(w) * wz
  }), list( .lshape = lshape , .lscale = lscale,
            .eshape = eshape,  .escale = escale,
            .tol12 = tol12, .nsimEIM = nsimEIM ))))
}










dexpgeom <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(scale), length(shape))
  x <- rep(x, len = N)
  scale <- rep(scale, len = N)
  shape <- rep(shape, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- -x[xok] / scale[xok]
    logdensity[xok] <- -log(scale[xok]) + log1p(-shape[xok]) + 
                       temp1 - 2 * log1p(-shape[xok] * exp(temp1))
  }

  logdensity[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}


pexpgeom <- function(q, scale = 1, shape) {
  temp1 <- -q / scale
  ans <- -expm1(temp1) / (1 - shape * exp(temp1))
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}

 
qexpgeom <- function(p, scale = 1, shape) {
  ans <- (-scale) * log((p - 1) / (p * shape - 1))
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}


rexpgeom <- function(n, scale = 1, shape) {
  ans <- qexpgeom(runif(n), shape = shape, scale = scale)
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}






expgeometric.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 expgeometric <- function(lscale = "loge", lshape = "logit",
                          iscale = NULL,   ishape = NULL, 
                          tol12 = 1.0e-05, zero = 1,
                          nsimEIM = 400) {


  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE) || any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential geometric distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(shape - 1) * log(1 - ",
            "shape) / (shape / scale)"), 
                           
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
         nsimEIM = .nsimEIM ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lscale = lscale, .lshape = lshape,
           .nsimEIM = nsimEIM ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <- c(
      namesof("scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init <- if (is.Numeric( .iscale , positive = TRUE)) {
                      rep( .iscale , len = n)
                    } else {
                      stats::sd(c(y))  # The papers scale parameter beta
                    }

      shape.init <- if (is.Numeric( .ishape , positive = TRUE)) {
                      rep( .ishape , len = n)
                    } else {
                      rep(2 - exp(median(y)/scale.init), len = n)
                    }
      shape.init[shape.init >= 0.95] <- 0.95
      shape.init[shape.init <= 0.05] <- 0.05

      
      etastart <-
        cbind(theta2eta(scale.init, .lscale , earg = .escale ),
              theta2eta(shape.init, .lshape , earg = .eshape ))

    }
   }), list( .lscale = lscale, .lshape = lshape, 
             .iscale = iscale, .ishape = ishape, 
             .escale = escale, .eshape = eshape))), 

  linkinv = eval(substitute(function(eta, extra = NULL) {
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    
    (shape - 1) * log1p(-shape) / (shape / Scale)

  }, list( .lscale = lscale, .lshape = lshape, 
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale , shape = .lshape )

    misc$earg <- list(scale = .escale , shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dexpgeom(x = y, scale = Scale, shape = shape,
                                 log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale , .lshape = lshape , 
           .escale = escale , .eshape = eshape ))), 
      
  vfamily = c("expgeometric"),

  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

     temp2 <- exp(-y / Scale)
     temp3 <- shape * temp2
     temp4 <- y / Scale^2
     dl.dscale <-  -1 / Scale + temp4 + 2 * temp4 * temp3 / (1 - temp3)
     dl.dshape <- -1 / (1 - shape)    + 2 * temp2 / (1 - temp3)

    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )            
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas <- cbind(dscale.deta, dshape.deta)

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

  weight = eval(substitute(expression({
  








      run.varcov <- 0
      ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

      if (length( .nsimEIM )) {
        for (ii in 1:( .nsimEIM )) {
          ysim <- rexpgeom(n, scale=Scale, shape=shape)

          temp2 <- exp(-ysim / Scale)
          temp3 <- shape * temp2
          temp4 <- ysim / Scale^2
          dl.dscale <-  -1 / Scale + temp4 + 
                       2 * temp4 * temp3 / (1 - temp3)
          dl.dshape <- -1 / (1 - shape) + 
                       2 * temp2 / (1 - temp3)

          temp6 <- cbind(dl.dscale, dl.dshape)
          run.varcov <- run.varcov +
              temp6[,ind1$row.index] * temp6[,ind1$col.index]
      }

      run.varcov <- run.varcov / .nsimEIM

      wz <- if (intercept.only)
              matrix(colMeans(run.varcov),
                     n, ncol(run.varcov), byrow = TRUE) else run.varcov

      wz <- wz * dthetas.detas[, ind1$row] *
                 dthetas.detas[, ind1$col]
    }

    c(w) * wz      
  }), list( .nsimEIM = nsimEIM ))))
}











dexplog <- function(x, scale = 1, shape, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  N <- max(length(x), length(scale), length(shape))
  x     <- rep(x,     len = N)
  scale <- rep(scale, len = N)
  shape <- rep(shape, len = N)

  logdensity <- rep(log(0), len = N)
  if (any(xok <- (x > 0))) {
    temp1 <- -x[xok] / scale[xok]
    logdensity[xok] <- -log(-log(shape[xok])) - log(scale[xok]) + 
                       log1p(-shape[xok]) + temp1 - 
                       log1p(-(1-shape[xok]) * exp(temp1))
  }

  logdensity[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  if (log.arg) {
    logdensity
  } else {
     exp(logdensity)
  }
}


pexplog <- function(q, scale = 1, shape) {
  ans <- 1 - log1p(-(1-shape) * exp(-q / scale)) / log(shape)
  ans[q <= 0] <- 0
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}



qexplog <- function(p, scale = 1, shape) {


  ans <- -scale * (log1p(-shape^(1.0 - p)) - log1p(-shape))

  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans[p == 0] <- 0
  ans[p == 1] <- Inf
  ans
}



rexplog <- function(n, scale = 1, shape) {
  ans <- qexplog(runif(n), scale = scale, shape = shape)
  ans[(scale <= 0) | (shape <= 0) | (shape >= 1)] <- NaN
  ans
}









explogff.control <- function(save.weights = TRUE, ...) {
    list(save.weights = save.weights)
}


 explogff <- function(lscale = "loge", lshape = "logit",
                      iscale = NULL,   ishape = NULL,
                      tol12 = 1.0e-05, zero = 1,
                      nsimEIM = 400) {

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lshape <- as.list(substitute(lshape))
  eshape <- link2list(lshape)
  lshape <- attr(eshape, "function.name")


  if (length(ishape))
    if (!is.Numeric(ishape, positive = TRUE) ||
        any(ishape >= 1))
      stop("bad input for argument 'ishape'")

  if (length(iscale))
    if (!is.Numeric(iscale, positive = TRUE))
    stop("bad input for argument 'iscale'")



  if (!is.Numeric(nsimEIM, length.arg = 1,
                  integer.valued = TRUE))
      stop("bad input for argument 'nsimEIM'")
  if (nsimEIM <= 50)
      stop("argument 'nsimEIM' should be an integer greater than 50")


  new("vglmff",
  blurb = c("Exponential logarithmic distribution\n\n",
            "Links:    ",
            namesof("scale", lscale, earg = escale), ", ",
            namesof("shape", lshape, earg = eshape), "\n",
            "Mean:     ", "(-polylog(2, 1 - p) * scale) / log(shape)"),

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
         nsimEIM = .nsimEIM ,
         lscale = .lscale ,
         lshape = .lshape ,
         zero = .zero )
  }, list( .zero = zero, .lscale = lscale, .lshape = lshape,
           .nsimEIM = nsimEIM ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <- c(
      namesof("scale", .lscale , earg = .escale , short = TRUE),
      namesof("shape", .lshape , earg = .eshape , short = TRUE))

    if (!length(etastart)) {

      scale.init <- if (is.Numeric( .iscale , positive = TRUE)) {
                     rep( .iscale , len = n)
                   } else {
                     stats::sd(c(y))  
                   }

      shape.init <- if (is.Numeric( .ishape , positive = TRUE)) {
                     rep( .ishape , len = n)
                   } else {
                      rep((exp(median(y)/scale.init) - 1)^2, len = n)
                   }
      shape.init[shape.init >= 0.95] <- 0.95
      shape.init[shape.init <= 0.05] <- 0.05


      etastart <-
        cbind(theta2eta(scale.init, .lscale , earg = .escale ),
              theta2eta(shape.init, .lshape , earg = .eshape ))

    }
   }), list( .lscale = lscale, .lshape = lshape,
             .iscale = iscale, .ishape = ishape,
             .escale = escale, .eshape = eshape))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )



    qexplog(p = 0.5, shape = shape, scale = scale)  

  }, list( .lscale = lscale, .lshape = lshape,
           .escale = escale, .eshape = eshape ))),

  last = eval(substitute(expression({
    misc$link <-    c(scale = .lscale , shape = .lshape )

    misc$earg <- list(scale = .escale , shape = .eshape )

    misc$expected <- TRUE
    misc$nsimEIM <- .nsimEIM
    misc$multipleResponses <- FALSE
  }), list( .lscale = lscale, .lshape = lshape,
            .escale = escale, .eshape = eshape,
            .nsimEIM = nsimEIM ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )
    

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dexplog(x = y, scale = Scale,
                                shape = shape, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lscale = lscale , .lshape = lshape ,
           .escale = escale , .eshape = eshape ))),

  vfamily = c("explogff"),

  deriv = eval(substitute(expression({
    Scale <- eta2theta(eta[, 1], .lscale , earg = .escale )
    shape <- eta2theta(eta[, 2], .lshape , earg = .eshape )

     temp2 <- exp(-y / Scale)
     temp3 <- y / Scale^2
     temp4 <- 1 - shape
     dl.dscale <- (-1 / Scale) + temp3 + (temp4 * temp3 *
                  temp2) / (1 - temp4 * temp2)
     dl.dshape <- -1 / (shape * log(shape)) - 1 / temp4 -
                  temp2 / (1 - temp4 * temp2)

    dscale.deta <- dtheta.deta(Scale, .lscale , earg = .escale )
    dshape.deta <- dtheta.deta(shape, .lshape , earg = .eshape )
    dthetas.detas <- cbind(dscale.deta, dshape.deta)

    answer <- c(w) * cbind(dl.dscale, dl.dshape) * dthetas.detas
    answer
  }), list( .lscale = lscale , .lshape = lshape,
            .escale = escale,  .eshape = eshape ))),

  weight = eval(substitute(expression({



        run.varcov <- 0
        ind1 <- iam(NA, NA, M = M, both = TRUE, diag = TRUE)

        if (length( .nsimEIM )) {
          for (ii in 1:( .nsimEIM )) {
            ysim <- rexplog(n, scale = Scale, shape = shape)

            temp2 <- exp(-ysim / Scale)
            temp3 <- ysim / Scale^2
            temp4 <- 1 - shape
            dl.dscale <- (-1 / Scale) + temp3 + (temp4 * temp3 *
                         temp2) / (1 - temp4 * temp2)
            dl.dshape <- -1 / (shape * log(shape)) - 1 / temp4 -
                         temp2 / (1 - temp4 * temp2)

            temp6 <- cbind(dl.dscale, dl.dshape)
            run.varcov <- run.varcov +
                       temp6[,ind1$row.index] *
                       temp6[,ind1$col.index]
          }

          run.varcov <- run.varcov / .nsimEIM

          wz <- if (intercept.only)
                matrix(colMeans(run.varcov),
                       n, ncol(run.varcov), byrow = TRUE) else
                run.varcov

          wz <- wz * dthetas.detas[, ind1$row] *
                    dthetas.detas[, ind1$col]
        }

    c(w) * wz
  }), list( .nsimEIM = nsimEIM ))))
}











  
dweibull3 <- function(x, location = 0, scale = 1, shape,
                      log = FALSE) {

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  dweibull(x = x - location, shape = shape,
           scale = scale, log = log.arg)
}

pweibull3 <- function(q, location = 0, scale = 1, shape) {
  pweibull(q = q - location, scale = scale, shape = shape)
}


qweibull3 <- function(p, location = 0, scale = 1, shape) {
  location + qweibull(p = p, shape = shape, scale = scale)
}


rweibull3 <- function(n, location = 0, scale = 1, shape) {
  location + rweibull(n = n, shape = shape, scale = scale)
}









   ### Two-piece normal (TPN) family 


dtpn <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {


  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
           na.rm = TRUE))
    stop("some parameters out of bound")

  LLL <- max(length(x), length(location), length(scale),
            length(skewpar))
  if (length(x) != LLL) x <- rep(x, length = LLL)
  if (length(location) != LLL) location <- rep(location, length = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar <- rep(skewpar, length = LLL)
    
  zedd <- (x - location) / scale

  log.s1 <-  -zedd^2 / (8 * skewpar^2)
  log.s2 <-  -zedd^2 / (8 * (1 - skewpar)^2)
            
  logdensity <- log.s1
  logdensity[zedd > 0] <- log.s2[zedd > 0]
  
  logdensity <- logdensity -log(scale) - log(sqrt(2 * pi))

  if (log.arg) logdensity else exp(logdensity)
}

ptpn <- function(q, location = 0, scale = 1, skewpar = 0.5) {

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")


 zedd <- (q - location) / scale

  s1 <- 2 * skewpar * pnorm(zedd, sd = 2 * skewpar)  #/ scale
  s2 <- skewpar + (1 - skewpar) *
        pgamma(zedd^2 / (8 * (1-skewpar)^2), 0.5)
 
ans <- rep(0.0, length(zedd))
ans[zedd <= 0] <- s1[zedd <= 0]
ans[zedd > 0] <- s2[zedd > 0]

ans
}



pos <- function(x) ifelse(x > 0, x, 0.0)
 

qtpn <- function(p, location = 0, scale = 1, skewpar = 0.5) {

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
    stop("some parameters out of bound")
    # Recycle the vectors to equal lengths
  LLL <- max(length(pp), length(location), length(scale),
            length(skewpar))
  if (length(pp) != LLL) pp <- rep(pp, length = LLL)
  if (length(location) != LLL) location <- rep(location, length = LLL)
  if (length(scale) != LLL) scale <- rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar <- rep(skewpar, length = LLL)
       
  qtpn <- rep(NA_real_, length(LLL))
  qtpn <- qnorm(pp / (2 * skewpar), sd = 2 * skewpar)
  qtpn[pp > skewpar] <- sqrt(8 * ( 1 - skewpar)^2 * 
                        qgamma(pos( pp - skewpar) / ( 
                        1 - skewpar),.5))[pp > skewpar]
        
   qtpn * scale + location
  
}





rtpn <- function(n, location = 0, scale = 1, skewpar = 0.5) {


  qtpn(p = runif(n), location = location,
       scale = scale, skewpar = skewpar)
}





tpnff <- function(llocation = "identitylink", lscale = "loge",
                  pp = 0.5, method.init = 1,  zero = 2) {
  if (!is.Numeric(method.init, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
     stop("argument 'imethod' must be 1 or 2 or 3 or 4")

  if (!is.Numeric(pp, length.arg = 1, positive = TRUE))
    stop("bad input for argument 'pp'")


  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")



  new("vglmff",
  blurb = c("Two-piece normal distribution \n\n",
            "Links: ",
            namesof("location",  llocat,  earg = elocat), ", ",
            namesof("scale",     lscale,  earg = escale), "\n\n",
            "Mean: "),
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
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocat,
           .lscale = lscale ))),

  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    predictors.names <-
       c(namesof("location", .llocat , earg = .elocat , tag = FALSE),
         namesof("scale",    .lscale , earg = .escale , tag = FALSE))




    if (!length(etastart)) {
        junk <- lm.wfit(x = x, y = c(y), w = c(w))
        scale.y.est <-
          sqrt( sum(c(w) * junk$resid^2) / junk$df.residual )
        location.init <- if ( .llocat == "loge")
          pmax(1/1024, y) else {

        if ( .method.init == 3) {
          rep(weighted.mean(y, w), len = n)
        } else if ( .method.init == 2) {
          rep(median(rep(y, w)), len = n)
        } else if ( .method.init == 1) {
          junk$fitted
        } else {
          y
        }
      }
      etastart <- cbind(
           theta2eta(location.init,  .llocat , earg = .elocat ),
           theta2eta(scale.y.est,    .lscale , earg = .escale ))
    }
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .method.init = method.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat , earg = .elocat )
  }, list( .llocat = llocat,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link     <-    c("location" = .llocat , "scale" = .lscale )

    misc$earg     <- list("location" = .elocat , "scale" = .escale )

    misc$expected <- TRUE
    misc$pp       <- .pp
    misc$method.init <- .method.init
    misc$multipleResponses <- FALSE
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .pp     = pp,        .method.init = method.init ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    location <- eta2theta(eta[, 1], .llocat , earg = .elocat )
    myscale  <- eta2theta(eta[, 2], .lscale , earg = .escale )
    ppay     <- .pp
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtpn(y, skewpar = ppay, location = location,
                             scale = myscale, log.arg = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .llocat = llocat, .lscale = lscale,
           .elocat = elocat, .escale = escale,
           .pp      = pp ))),
  vfamily = c("tpnff"),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta[, 1], .llocat ,  earg = .elocat )
    myscale <- eta2theta(eta[, 2], .lscale ,  earg = .escale )
    mypp    <- .pp

    zedd <- (y - mylocat) / myscale
 #   cond1 <-    (zedd <= 0)
     cond2 <-    (zedd > 0)

    dl.dlocat        <-  zedd / (4 * mypp^2)  # cond1
    dl.dlocat[cond2] <- (zedd / (4 * (1 - mypp)^2))[cond2]
    dl.dlocat        <- dl.dlocat / myscale

    dl.dscale        <-  zedd^2 / (4 * mypp^2)
    dl.dscale[cond2] <- (zedd^2 / (4 * (1 - mypp)^2))[cond2]
    dl.dscale        <- (-1 + dl.dscale) / myscale

    #dl.dpp        <-  zedd^2 /  (4 * mypp^3)
    #dl.dpp[cond2] <- -zedd^2 /  (4 * (1 - mypp)^3)[cond2]
    
    dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
    dscale.deta <- dtheta.deta(myscale, .lscale, earg = .escale)

    ans <- c(w) * cbind(dl.dlocat * dlocat.deta,
                        dl.dscale * dscale.deta)
    ans
  }), list( .llocat = llocat, .lscale = lscale,
            .elocat = elocat, .escale = escale,
            .pp      = pp ))),
  weight = eval(substitute(expression({
    wz   <- matrix(0, n, M)  # diag matrix; y is one-col too
    temp10 <- mypp * (1 - mypp)
    ned2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
    ned2l.dscale2        <- 2 /  myscale^2
     

    wz[, iam(1, 1, M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2, M)] <- ned2l.dscale2 * dscale.deta^2
  # wz[, iam(3, 3, M)] <- ned2l.dskewpar2 * dskewpa.deta^2
  # wz[, iam(1, 3, M)] <- ned2l.dlocatdskewpar * dskewpar.deta * dlocat.deta
    c(w) * wz
  }))))
}



  ########################################################################


tpnff3 <- function(llocation = "identitylink",
                    lscale   = "loge",
                    lskewpar = "identitylink",
                    method.init = 1,  zero = 2)
{
  if (!is.Numeric(method.init, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      method.init > 4)
    stop("argument 'imethod' must be 1 or 2 or 3 or 4")



  llocat <- as.list(substitute(llocation))
  elocat <- link2list(llocat)
  llocat <- attr(elocat, "function.name")

  lscale <- as.list(substitute(lscale))
  escale <- link2list(lscale)
  lscale <- attr(escale, "function.name")

  lskewp <- as.list(substitute(lskewpar))
  eskewp <- link2list(lskewp)
  lskewp <- attr(eskewp, "function.name")






  new("vglmff",
  blurb = c("Two-piece normal distribution \n\n",
            "Links: ",
            namesof("location", llocat, earg = elocat), ", ",
            namesof("scale",    lscale, earg = escale),  ", ",
            namesof("skewpar",  lskewp, earg = eskewp),  "\n\n",
            "Mean: "),
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
         parameters.names = c("location", "scale", "skewpar"),
         llocation = .llocat ,
         lscale    = .lscale ,
         lskewpar  = .lskewp ,
         zero = .zero )
  }, list( .zero = zero,
           .llocat = llocat,
           .lscale = lscale,
           .lskewp = lskewp ))),


  initialize = eval(substitute(expression({

    temp5 <-
    w.y.check(w = w, y = y,
              out.wy = TRUE,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    predictors.names <-
       c(namesof("location", .llocat, earg = .elocat, tag = FALSE),
         namesof("scale",    .lscale, earg = .escale, tag = FALSE),
         namesof("skewpar",  .lskewp, earg = .eskewp, tag = FALSE))

    if (!length(etastart)) {
      junk = lm.wfit(x = x, y = c(y), w = c(w))
      scale.y.est <- sqrt(sum(c(w) * junk$resid^2) / junk$df.residual)
      location.init <- if ( .llocat == "loge") pmax(1/1024, y) else {
        if ( .method.init == 3) {
          rep(weighted.mean(y, w), len = n)
        } else if ( .method.init == 2) {
          rep(median(rep(y, w)), len = n)
        } else if ( .method.init == 1) {
          junk$fitted
        } else {
          y
        }
      }
      skew.l.in <- sum((y < location.init)) / length(y)
      etastart <- cbind(
           theta2eta(location.init, .llocat,   earg = .elocat),
           theta2eta(scale.y.est,   .lscale,   earg = .escale),
           theta2eta(skew.l.in,     .lskewp, earg = .escale))
    }
  }), list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
            .elocat = elocat, .escale = escale, .eskewp = eskewp,
            
            .method.init=method.init ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta[, 1], .llocat, earg = .elocat)
  }, list( .llocat = llocat,
           .elocat = elocat, .escale = escale ))),
  last = eval(substitute(expression({
    misc$link     <-     c("location" = .llocat,
                           "scale"    = .lscale, 
                           "skewpar"  = .lskewp)

    misc$earg     <-  list("location" = .elocat,
                           "scale"    = .escale,
                           "skewpar"  = .eskewp)

    misc$expected <- TRUE
         misc$method.init <- .method.init
  }), list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
            .elocat = elocat, .escale = escale, .eskewp = eskewp,
                    .method.init = method.init ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

   locat    <- eta2theta(eta[, 1], .llocat , earg = .elocat )
   myscale  <- eta2theta(eta[, 2], .lscale , earg = .escale )
   myskew   <- eta2theta(eta[, 3], .lskewp , earg = .eskewp )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dtpn(y, location = locat,  scale = myscale,
                             skewpar = myskew, log.arg = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
   }
  }, list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
           .elocat = elocat, .escale = escale, .eskewp = eskewp
           ))),
  vfamily = c("tpnff3"),
  deriv = eval(substitute(expression({
    mylocat <- eta2theta(eta[, 1], .llocat,   earg = .elocat)
    myscale <- eta2theta(eta[, 2], .lscale,   earg = .escale)
    myskew  <- eta2theta(eta[, 3], .lskewp, earg = .eskewp)
  

    zedd <- (y - mylocat) / myscale
   cond2 <-    (zedd > 0)

    dl.dlocat        <-  zedd / (4 * myskew^2)  # cond1
    dl.dlocat[cond2] <- (zedd / (4 * (1 - myskew)^2))[cond2]
    dl.dlocat        <- dl.dlocat / myscale

    dl.dscale        <-  zedd^2 / (4 * myskew^2)
    dl.dscale[cond2] <- (zedd^2 / (4 * (1 - myskew)^2))[cond2]
    dl.dscale        <- (-1 + dl.dscale) / myscale

    dl.dskewpar      <-     zedd^2 /  (4 * myskew^3)
    dl.dskewpar[cond2] <- (-zedd^2 /  (4 * (1 - myskew)^3))[cond2]
    


    dlocat.deta <- dtheta.deta(mylocat, .llocat, earg = .elocat)
    dscale.deta <- dtheta.deta(myscale, .lscale, earg = .escale)
    dskewpar.deta <- dtheta.deta(myskew, .lskewp, earg = .eskewp)
    ans <-
    c(w) * cbind(dl.dlocat * dlocat.deta,
              dl.dscale * dscale.deta,
              dl.dskewpar * dskewpar.deta
              )
    ans
  }), list( .llocat = llocat, .lscale = lscale, .lskewp = lskewp,
            .elocat = elocat, .escale = escale, .eskewp = eskewp
            ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, n, dimm(M))  # diag matrix; y is one-col too
   
    temp10 <- myskew * (1 - myskew)

    ned2l.dlocat2        <- 1 / ((4 * temp10) * myscale^2)
    ned2l.dscale2        <- 2 /  myscale^2
    ned2l.dskewpar2      <- 3 / temp10
    ned2l.dlocatdskewpar <- (-2 * sqrt(2)) / (temp10 * sqrt(pi) *
                             myscale)
     
    wz[, iam(1, 1,M)] <- ned2l.dlocat2 * dlocat.deta^2
    wz[, iam(2, 2,M)] <- ned2l.dscale2 * dscale.deta^2
    wz[, iam(3, 3,M)] <- ned2l.dskewpar2 * dskewpar.deta^2
    wz[, iam(1, 3,M)] <- ned2l.dlocatdskewpar * dskewpar.deta *
                         dlocat.deta
  
    ans
    c(w) * wz
  }))))
}




dozibeta <- function(x, shape1, shape2, pobs0 = 0,
                     pobs1 = 0, log = FALSE, tol = .Machine$double.eps) {
  log.arg <- log
  rm(log)
  LLL <- max(length(x), length(shape1),
             length(shape2), length(pobs0), length(pobs1))
  if (LLL != length(x))
    x <- rep(x, length = LLL)
  if (LLL != length(shape1))
    shape1 <- rep(shape1, length = LLL)
  if (LLL != length(shape2))
    shape2 <- rep(shape2, length = LLL)
  if (LLL != length(pobs0))
    pobs0 <- rep(pobs0, length = LLL)
  if (LLL != length(pobs1))
    pobs1 <- rep(pobs1, length = LLL)
  ans <- rep(NA, length = LLL)
  k1 <- (pobs0 < -tol | pobs1 < -tol |
    (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans[!k4 & !k1] <- dbeta(x[!k4 & !k1], 
                          shape1[!k4 & !k1], 
                          shape2[!k4 & !k1], log = TRUE) + 
                    log1p(-(pobs0[!k4 & !k1] + pobs1[!k4 & !k1]))
  k2 <- x == 0 & pobs0 > 0 & !is.na(x)
  k3 <- x == 1 & pobs1 > 0 & !is.na(x)
  ans[k2 & !k4 & !k1] <- log(pobs0[k2 & !k4 & !k1])
  ans[k3 & !k4 & !k1] <- log(pobs1[k3 & !k4 & !k1])
  if (!log.arg) ans <- exp(ans)
  if (any(k1 & !k4)) {
    ans[k1 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}


rozibeta <- function(n, shape1, shape2, pobs0 = 0, pobs1 = 0,
                     tol = .Machine$double.eps) {
  use.n <- if ((length.n <- length(n)) > 1) {
    length.n
  } else {
    if (!is.Numeric(n, integer.valued = TRUE, length.arg = 1, 
                    positive = TRUE)) {
      stop("bad input for argument 'n'")
    } else {
      n
    }
  }
  shape1 <- rep(shape1, length.out = use.n)
  shape2 <- rep(shape2, length.out = use.n)
  pobs0 <- rep(pobs0, length.out = use.n)
  pobs1 <- rep(pobs1, length.out = use.n)
  random.number <- runif(use.n)
  ans <- rep(NA, length = use.n)
  k5 <- (pobs0 < -tol | pobs1 < -tol |
           (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans[!k4] <- qozibeta(random.number[!k4], shape1 = shape1,
                       shape2 = shape2, pobs0 = pobs0,
                       pobs1 = pobs1)
  if (any(k5 & !k4)) {
    ans[k5 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}


pozibeta <- function(q, shape1, shape2, pobs0 = 0, pobs1 = 0,
                     lower.tail = TRUE, log.p = FALSE,
                     tol = .Machine$double.eps) {
  LLL <- max(length(q), length(shape1),
             length(shape2), length(pobs0), length(pobs1))
  if (LLL != length(q))
    q <- rep(q, length = LLL)
  if (LLL != length(shape1))
    shape1 <- rep(shape1, length = LLL)
  if (LLL != length(shape2))
    shape2 <- rep(shape2, length = LLL)
  if (LLL != length(pobs0))
    pobs0 <- rep(pobs0, length = LLL)
  if (LLL != length(pobs1))
    pobs1 <- rep(pobs1, length = LLL)
  k3 <- (pobs0 < -tol | pobs1 < -tol |
           (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans <- rep(NA, length = LLL)
  ans[!k3 & !k4] <- pbeta(q[!k3 & !k4],
                          shape1[!k3 & !k4], 
                          shape2[!k3 & !k4], log.p = TRUE) +
    log1p(-(pobs0[!k3 & !k4] + pobs1[!k3 & !k4]))
  ans <- exp(ans)
  k1 <- q >= 0 & !is.na(q)
  k2 <- q >= 1 & !is.na(q)
  ans[k1 & !k3 & !k4] <- ans[k1 & !k3 & !k4] + 
    pobs0[k1 & !k3 & !k4]
  ans[k2 & !k3 & !k4] <- ans[k2 & !k3 & !k4] + 
    pobs1[k2 & !k3 & !k4]
  if (!lower.tail & log.p) {
    ans <- log1p(-ans)
  } else {
    if (!lower.tail)
      ans <- 1 - ans
    if (log.p)
      ans <- log(ans)
  }
  if (any(k3 & !k4)) {
    ans[k3 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}


qozibeta <- function(p, shape1, shape2, pobs0 = 0, pobs1 = 0,
                     lower.tail = TRUE, log.p = FALSE,
                     tol = .Machine$double.eps) {
  LLL <- max(length(p), length(shape1),
             length(shape2), length(pobs0), length(pobs1))
  if (LLL != length(p))
    p <- rep(p, length = LLL)
  if (LLL != length(shape1))
    shape1 <- rep(shape1, length = LLL)
  if (LLL != length(shape2))
    shape2 <- rep(shape2, length = LLL)
  if (LLL != length(pobs0))
    pobs0 <- rep(pobs0, length = LLL)
  if (LLL != length(pobs1))
    pobs1 <- rep(pobs1, length = LLL)
  k0 <- (pobs0 < -tol | pobs1 < -tol |
           (pobs0 + pobs1) > (1 + tol))
  k4 <- is.na(pobs0) | is.na(pobs1)
  ans <- rep(NA, length = LLL)
  if (!lower.tail & log.p) {
    p <- -expm1(p)
  } else{
    if (!lower.tail)
      p <- 1 - p
    if (log.p) {
      p <- exp(p)
    }
  }
  k1 <- p >= 0 & p <= pobs0 & !is.na(p)
  k2 <- p > pobs0 & p < (1 - pobs1) & !is.na(p)
  k3 <- p >= (1 - pobs1) & p <= 1 & !is.na(p)
  ans[k1 & !k0 & !k4] <- 0
  ans[k2 & !k0 & !k4] <-
    qbeta((p[k2 & !k0 & !k4] -
           pobs0[k2 & !k0 & !k4]) / (1 - pobs0[k2 & !k0 & !k4] -
           pobs1[k2 & !k0 & !k4]),
           shape1 = shape1[k2 & !k0 & !k4], 
           shape2 = shape2[k2 & !k0 & !k4])
  ans[k3 & !k0 & !k4] <- 1
  if (any(k0 & !k4)) {
    ans[k3 & !k4] <- NaN
    warning("NaNs produced")
  }
  ans
}





log1mexp <- function(x) {
  if (any(x < 0 & !is.na(x)))
    stop("Inputs need to be non-negative!")
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}


log1pexp <- function(x){
  
  ifelse(x <= -37, exp(x),
         ifelse(x <= 18, log1p(exp(x)),
                ifelse(x <= 33, x + exp(-x), x)))
}






dozibetabinom.ab <- function(x, size, shape1, shape2, pstr0 = 0,
                             pstrsize = 0, log = FALSE) {
  log.arg <- log
  rm(log)
  LLL <- max(length(x), length(size), length(shape1),
             length(shape2), length(pstr0), length(pstrsize))
  if (LLL != length(x))
    x <- rep(x, length = LLL)
  if (LLL != length(size))
    size <- rep(size, length = LLL)
  if (LLL != length(shape1))
    shape1 <- rep(shape1, length = LLL)
  if (LLL != length(shape2))
    shape2 <- rep(shape2, length = LLL)
  if (LLL != length(pstr0))
    pstr0 <- rep(pstr0, length = LLL)
  if (LLL != length(pstrsize))
    pstrsize <- rep(pstrsize, length = LLL)
  ans <- rep(NA, length = LLL)
  k1 <- pstr0 < 0 | pstrsize < 0 |
           (pstr0 + pstrsize) > 1
  k <- is.na(size) | is.na(shape1) | is.na(shape2) |
    is.na(pstr0) | is.na(pstrsize) | is.na(x)
  if (sum(!k & !k1) > 0) {
    ans[!k & !k1] <-
      dbetabinom.ab(x[!k & !k1], size[!k & !k1], shape1[!k & !k1],
                    shape2[!k & !k1], log = TRUE) + 
      log1p(-(pstr0[!k & !k1]+pstrsize[!k & !k1]))
    if (!log.arg) ans <- exp(ans)
  }
  k2 <- x == 0 & pstr0 > 0 
  k3 <- x == size & pstrsize > 0 
  if (sum(k2 & !k & !k1) > 0)
    ans[k2 & !k & !k1] <- pstr0[k2 & !k & !k1] +
      ans[k2 & !k & !k1]
  if (sum(k3 & !k & !k1) > 0)
    ans[k3 & !k & !k1] <- pstrsize[k3 & !k & !k1] +
      ans[k3 & !k & !k1]
  if (any(k1 & !k)) {
    ans[k1 & !k] <- NaN
    warning("NaNs produced")
  }
  ans
}



dozibetabinom <- function(x, size, prob, rho = 0, pstr0 = 0,
                          pstrsize = 0, log = FALSE) {
  dozibetabinom.ab(x, size, shape1 = prob * (1 - rho) / rho,
                   shape2 = (1 - prob) * (1 - rho) / rho, 
                   pstr0 = pstr0, pstrsize = pstrsize, log = log)
}



rozibetabinom.ab <- function(n, size, shape1, shape2, 
                             pstr0 = 0, pstrsize = 0) {
  use.n <- if ((length.n <- length(n)) > 1) {
    length.n
  } else {
    if (!is.Numeric(n, integer.valued = TRUE, length.arg = 1, 
                    positive = TRUE)) {
      stop("bad input for argument 'n'")
    } else {
      n
    }
  }
  size <- rep(size, length.out = use.n)
  shape1 <- rep(shape1, length.out = use.n)
  shape2 <- rep(shape2, length.out = use.n)
  pstr0 <- rep(pstr0, length.out = use.n)
  pstrsize <- rep(pstrsize, length.out = use.n)
  k <- is.na(size) | is.na(shape1) | is.na(shape2) |
    is.na(pstr0) | is.na(pstrsize)
  ans <- rep(NA, length = use.n)
  k1 <- pstr0 < 0 | pstrsize < 0 |
    (pstr0 + pstrsize) > 1
  random.number <- runif(use.n)
  k2 <- random.number[!k] < pstr0[!k]
  k3 <- pstr0[!k] <= random.number[!k] & 
    random.number[!k] <= (1 - pstrsize[!k])
  k4 <- (1 - pstrsize[!k]) < random.number[!k]
  if (sum(k2 & !k1 & !k) > 0)
    ans[k2 & !k1 & !k] <- 0
  if (sum(k3 & !k1 & !k) > 0)
    ans[k3 & !k1 & !k] <- rbetabinom.ab(sum(k3 & !k1 & !k), 
                                        size =  size[k3 & !k1 & !k],
                                        shape1 = shape1[k3 & !k1 & !k], 
                                        shape2 = shape2[k3 & !k1 & !k])
  if (sum(k4 & !k1 & !k) > 0)
    ans[k4 & !k1 & !k] <- size[k4 & !k1 & !k]
  ans
}



rozibetabinom <- function(n, size, prob, rho = 0, pstr0 = 0,
                          pstrsize = 0) {
  rozibetabinom.ab(n, size, shape1 = prob * (1 - rho) / rho,
                   shape2 = (1 - prob) * (1 - rho) / rho, 
                   pstr0 = pstr0,
                   pstrsize = pstrsize)
}



pozibetabinom.ab <- function(q, size, shape1, shape2, pstr0 = 0,
                             pstrsize = 0, lower.tail = TRUE,
                             log.p = FALSE) {
  LLL <- max(length(q), length(size), length(shape1),
             length(shape2), length(pstr0), length(pstrsize))
  if (LLL != length(q))
    q <- rep(q, length = LLL)
  if (LLL != length(size))
    size <- rep(size, length = LLL)
  if (LLL != length(shape1))
    shape1 <- rep(shape1, length = LLL)
  if (LLL != length(shape2))
    shape2 <- rep(shape2, length = LLL)
  if (LLL != length(pstr0))
    pstr0 <- rep(pstr0, length = LLL)
  if (LLL != length(pstrsize))
    pstrsize <- rep(pstrsize, length = LLL)
  k <- is.na(size) | is.na(shape1) | is.na(shape2) |
    is.na(pstr0) | is.na(pstrsize) | is.na(q)
  ans <- rep(NA, length = LLL)
  k1 <- pstr0 < 0 | pstrsize < 0 |
    (pstr0 + pstrsize) > 1
  if (sum(!k1 & !k) > 0)
    ans[!k & !k1] <-
      pbetabinom.ab(q[!k & !k1], size[!k & !k1], 
                    shape1[!k & !k1], shape2[!k & !k1], log.p = TRUE) +
      log1p(-(pstr0[!k & !k1] + pstrsize[!k & !k1]))
  ans <- exp(ans)
  k2 <- q >= 0 
  k3 <- q >= size 
  if (sum(k2 & !k1 & !k) > 0)
    ans[k2 & !k & !k1] <- ans[k2 & !k & !k1] + 
      pstr0[k2 & !k & !k1]
  if (sum(k3 & !k1 & !k) > 0)
    ans[k3 & !k & !k1] <- ans[k3 & !k & !k1] + 
      pstrsize[k3 & !k & !k1]
  if (!lower.tail & log.p) {
    ans <- log1p(-ans)
  } else {
    if (!lower.tail)
      ans <- 1 - ans
    if (log.p)
      ans <- log(ans)
  }
  if (any(!k & k1)) {
    ans[!k & k1] <- NaN
    warning("NaNs produced")
  }
  ans
}


pozibetabinom <- function(q, size, prob, rho, 
                          pstr0 = 0, pstrsize = 0,
                        lower.tail = TRUE, log.p = FALSE) {
  pozibetabinom.ab(q, size, shape1 = prob * (1 - rho) / rho,
                 shape2 = (1 - prob) * (1 - rho) / rho, 
                 pstr0 = pstr0, pstrsize = pstrsize,
                 lower.tail = lower.tail, log.p = log.p)
}













