# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.













dzanegbin <- function(x, size, prob = NULL, munb = NULL, pobs0 = 0,
                      log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- 1 / (1 + munb/size)
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(pobs0), length(prob), length(size))
  if (length(x)     != LLL) x     <- rep(x,     len = LLL)
  if (length(pobs0) != LLL) pobs0 <- rep(pobs0, len = LLL)
  if (length(prob)  != LLL) prob  <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size  <- rep(size,  len = LLL)

  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  if (!is.Numeric(prob, positive = TRUE) ||
      max(prob, na.rm = TRUE) >= 1)
    stop("argument 'prob' must be in (0,1)")
  if (!is.Numeric(size, positive = TRUE))
    stop("argument 'size' must be in (0,Inf)")
  index0 <- x == 0

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                    dposnegbin(x[!index0], prob = prob[!index0],
                               size = size[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1 - pobs0[!index0]) * dposnegbin(x[!index0],
                      prob = prob[!index0], size = size[!index0])
  }
  ans
}



pzanegbin <- function(q, size, prob = NULL, munb = NULL, pobs0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- 1 / (1 + munb/size)
  }

  LLL <- max(length(q), length(pobs0), length(prob), length(size))
  if (length(q)     != LLL) q     <- rep(q,     len = LLL)
  if (length(pobs0) != LLL) pobs0 <- rep(pobs0, len = LLL)
  if (length(prob)  != LLL) prob  <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size  <- rep(size,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  qindex <- (q >  0)
  ans[ qindex] <- pobs0[qindex] + (1 - pobs0[qindex]) *
                  pposnegbin(q[qindex], size = size[qindex],
                                        prob = prob[qindex])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)

  ans
}


qzanegbin <- function(p, size, prob = NULL, munb = NULL, pobs0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- 1 / (1 + munb/size)
  }

  LLL <- max(length(p), length(pobs0), length(prob), length(size))
  if (length(p)     != LLL) p      <- rep(p,     len = LLL)
  if (length(pobs0) != LLL) pobs0  <- rep(pobs0, len = LLL)
  if (length(prob)  != LLL) prob   <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size   <- rep(size,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ans <- p
  ans[p <= pobs0] <- 0
  pindex <- (p > pobs0)
  ans[pindex] <-
    qposnegbin((p[pindex] - pobs0[pindex]) / (1 - pobs0[pindex]),
               prob = prob[pindex],
               size = size[pindex])
  ans
}



rzanegbin <- function(n, size, prob = NULL, munb = NULL, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- 1 / (1 + munb/size)
  }

  ans <- rposnegbin(n = use.n, prob = prob, size = size)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")

  ifelse(runif(use.n) < pobs0, 0, ans)
}






dzapois <- function(x, lambda, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pobs0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                    dpospois(x[!index0], lambda[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1 - pobs0[!index0]) *
                    dpospois(x[!index0], lambda[!index0])
  }
  ans
}



pzapois <- function(q, lambda, pobs0 = 0) {
  LLL <- max(length(q), length(lambda), length(pobs0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL)
  ans <- rep(0.0, len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  ans[q >  0] <-    pobs0[q > 0] +
                 (1-pobs0[q > 0]) * ppospois(q[q > 0], lambda[q > 0])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)

  ans
}



qzapois <- function(p, lambda, pobs0 = 0) {
  LLL <- max(length(p), length(lambda), length(pobs0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL)

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ans <- p
  ind4 <- (p > pobs0)
  ans[!ind4] <- 0
  ans[ ind4] <- qpospois((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         lambda = lambda[ind4])
  ans
}


rzapois <- function(n, lambda, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  ans <- rpospois(use.n, lambda)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, length = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must in [0,1]")

  ifelse(runif(use.n) < pobs0, 0, ans)
}





dzipois <- function(x, lambda, pstr0 = 0, log = FALSE) {



  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(lambda), length(pstr0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL)

  ans <- x + lambda + pstr0



  index0 <- (x == 0)
  if (log.arg) {
    ans[ index0] <- log(pstr0[ index0] + (1 - pstr0[ index0]) *
                       dpois(x[ index0], lambda[ index0]))
    ans[!index0] <- log1p(-pstr0[!index0]) +
                   dpois(x[!index0], lambda[!index0], log = TRUE)
  } else {
    ans[ index0] <-      pstr0[ index0] + (1 - pstr0[ index0]) *
                       dpois(x[ index0], lambda[ index0])
    ans[!index0] <- (1 - pstr0[!index0]) *
                    dpois(x[!index0], lambda[!index0])
  }


  deflat.limit <- -1 / expm1(lambda)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



pzipois <- function(q, lambda, pstr0 = 0) {

  LLL <- max(length(pstr0), length(lambda), length(q))
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(q)      != LLL) q      <- rep(q,      len = LLL)

  ans <- ppois(q, lambda)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)


  deflat.limit <- -1 / expm1(lambda)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}



qzipois <- function(p, lambda, pstr0 = 0) {

  LLL <- max(length(p), length(lambda), length(pstr0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL)
  if (length(lambda) != LLL) lambda <- rep(lambda, len = LLL)
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL)
  ans    <- p

  ans[p <= pstr0] <- 0 
  pindex <- (p > pstr0)
  ans[pindex] <-
    qpois((p[pindex] - pstr0[pindex]) / (1 - pstr0[pindex]),
          lambda = lambda[pindex])


  deflat.limit <- -1 / expm1(lambda)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * exp(-lambda[ind0])
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * exp(-lambda[pindex])
    ans[pindex] <- qpospois((p[pindex] - Pobs0) / (1 - Pobs0),
                            lambda = lambda[pindex])
  }


  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans[p < 0] <- NaN
  ans[p > 1] <- NaN
  ans
}



rzipois <- function(n, lambda, pstr0 = 0) {

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  if (length(pstr0)  != use.n) pstr0  <- rep(pstr0,  len = use.n)
  if (length(lambda) != use.n) lambda <- rep(lambda, len = use.n)
 
  ans <- rpois(use.n, lambda)
  ans <- ifelse(runif(use.n) < pstr0, 0, ans)



  prob0 <- exp(-lambda)
  deflat.limit <- -1 / expm1(lambda)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- rpospois(sum(ind0), lambda[ind0]) 
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}









 yip88 <- function(link = "loge", n.arg = NULL, imethod = 1) {








  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("Zero-inflated Poisson (based on Yip (1988))\n\n",
            "Link:     ",
            namesof("lambda", link, earg), "\n",
            "Variance: (1 - pstr0) * lambda"),
  first = eval(substitute(expression({
    zero <- y == 0
    if (any(zero)) {
      if (length(extra)) extra$sumw <- sum(w) else
        extra <- list(sumw=sum(w))
      if (is.numeric(.n.arg) && extra$sumw != .n.arg) 
        stop("value of 'n.arg' conflicts with data ",
             "(it need not be specified anyway)")
      warning("trimming out the zero observations")


      axa.save <-  attr(x, "assign")
      x <- x[!zero,, drop = FALSE]
      attr(x, "assign") <- axa.save    # Don't lose these!!
      w <- w[!zero]
      y <- y[!zero]
    } else {
      if (!is.numeric(.n.arg)) 
        stop("n.arg must be supplied")
    }
        
  }), list( .n.arg = n.arg ))),

  initialize = eval(substitute(expression({
    narg <- if (is.numeric(.n.arg)) .n.arg else extra$sumw
    if (sum(w) > narg)
      stop("sum(w) > narg")

    w.y.check(w = w, y = y,
              ncol.w.max = 1,
              ncol.y.max = 1)


    predictors.names <-
      namesof("lambda", .link, list(theta = NULL), tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             pos.only = FALSE)
      etastart <- theta2eta(lambda.init, .link , earg = .earg )
    }
    if (length(extra)) {
      extra$sumw <- sum(w)
      extra$narg <- narg   # For @linkinv
    } else {
      extra <- list(sumw = sum(w), narg = narg)
    }
  }), list( .link = link, .earg = earg,
            .n.arg = n.arg, .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    lambda <- eta2theta(eta, .link, .earg)
    temp5 <- exp(-lambda)
    pstr0 <- (1 - temp5 - extra$sumw/extra$narg) / (1 - temp5)
    if (any(pstr0 <= 0))
      stop("non-positive value(s) of pstr0")
    (1 - pstr0) * lambda
  }, list( .link = link, .earg = earg ))),

  last = eval(substitute(expression({
    misc$link <-    c(lambda = .link )

    misc$earg <- list(lambda = .earg )

    if (intercept.only) {
      suma <- extra$sumw
      pstr0 <- (1 - temp5[1] - suma / narg) / (1 - temp5[1])
      pstr0 <- if (pstr0 < 0 || pstr0 > 1) NA else pstr0
      misc$pstr0 <- pstr0
    }
  }), list( .link = link, .earg = earg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .link)
    temp5 <- exp(-lambda)
    pstr0 <- (1 - temp5 - extra$sumw / extra$narg) / (1 - temp5)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
                 dzipois(x = y, pstr0 = pstr0, lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),

  vfamily = c("yip88"),
  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link , earg = .earg )
    temp5 <- exp(-lambda)
    dl.dlambda <- -1 + y/lambda - temp5/(1-temp5)
    dlambda.deta <- dtheta.deta(lambda, .link , earg = .earg )
    w * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    d2lambda.deta2 <- d2theta.deta2(lambda, .link , earg = .earg )
    d2l.dlambda2 <- -y / lambda^2 + temp5 / (1 - temp5)^2
    -w * (d2l.dlambda2*dlambda.deta^2 + dl.dlambda*d2lambda.deta2)
  }), list( .link = link, .earg = earg ))))
}





 zapoisson <-
  function(lpobs0 = "logit", llambda = "loge",
           type.fitted = c("mean", "lambda", "pobs0", "onempobs0"),
           imethod = 1,
           ipobs0 = NULL, ilambda = NULL, ishrinkage = 0.95,
           probs.y = 0.35,
           zero = NULL) {



  lpobs.0 <- as.list(substitute(lpobs0))
  epobs.0 <- link2list(lpobs.0)
  lpobs.0 <- attr(epobs.0, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs0", "onempobs0"))[1]



  new("vglmff",
  blurb = c("Zero-altered Poisson ",
            "(Bernoulli and positive-Poisson conditional model)\n\n",
            "Links:    ",
            namesof("pobs0",  lpobs.0, earg = epobs.0, tag = FALSE), ", ",
            namesof("lambda", llambda, earg = elambda, tag = FALSE), "\n",
            "Mean:     (1 - pobs0) * lambda / (1 - exp(-lambda))"),

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
         parameters.names = c("pobs0", "lambda"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)
    extra$dimnamesy <- dimnames(y)
    extra$type.fitted      <- .type.fitted

    mynames1 <- param.names("pobs0",  ncoly)
    mynames2 <- param.names("lambda", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpobs.0, earg = .epobs.0, tag = FALSE),
          namesof(mynames2, .llambda, earg = .elambda, tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             imu = .ilambda,
                             ishrinkage = .ishrinkage,
                             pos.only = TRUE,
                             probs.y = .probs.y )

      etastart <-
        cbind(theta2eta(if (length( .ipobs0 )) .ipobs0 else
                        (0.5 + w * y0) / (1 + w),
                        .lpobs.0 , earg = .epobs.0 ),
              theta2eta(lambda.init, .llambda , earg = .elambda ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lpobs.0 = lpobs.0, .llambda = llambda,
            .epobs.0 = epobs.0, .elambda = elambda,
            .ipobs0 = ipobs0, .ilambda = ilambda,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs0", "onempobs0"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1


    pobs.0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                              .lpobs.0, earg = .epobs.0 ))
    lambda <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                              .llambda, earg = .elambda ))


    ans <- switch(type.fitted,
                  "mean"      = (1 - pobs.0) * lambda / (-expm1(-lambda)),
                  "lambda"    = lambda,
                  "pobs0"     =      pobs.0,  # P(Y=0)
                  "onempobs0" =  1 - pobs.0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))       
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpobs.0 = lpobs.0, .llambda = llambda,
           .epobs.0 = epobs.0, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

    temp.names <- c(rep( .lpobs.0 , len = NOS),
                    rep( .llambda , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs.0
      misc$earg[[M1*ii  ]] <- .elambda
    }
  }), list( .lpobs.0 = lpobs.0, .llambda = llambda,
            .epobs.0 = epobs.0, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    pobs0  <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                              .lpobs.0, earg = .epobs.0))
    lambda <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                              .llambda, earg = .elambda ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dzapois(x = y, pobs0 = pobs0, lambda = lambda,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs.0 = lpobs.0, .llambda = llambda,
           .epobs.0 = epobs.0, .elambda = elambda ))),
  vfamily = c("zapoisson"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    pobs0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpobs.0 , earg = .epobs.0 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    rzapois(nsim * length(lambda), lambda = lambda, pobs0 = pobs0)
  }, list( .lpobs.0 = lpobs.0, .llambda = llambda,
           .epobs.0 = epobs.0, .elambda = elambda ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    phimat <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                              .lpobs.0, earg = .epobs.0 ))
    lambda <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                              .llambda, earg = .elambda ))

    dl.dlambda <- y / lambda + 1 / expm1(-lambda)
    dl.dphimat <- -1 / (1 - phimat)  # For y > 0 obsns

    for (spp. in 1:NOS) {
      dl.dphimat[skip[, spp.], spp.] <- 1 / phimat[skip[, spp.], spp.]
      dl.dlambda[skip[, spp.], spp.] <- 0
    }
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    mu.phi0 <- phimat

    temp3 <- if (.lpobs.0 == "logit") {
      c(w) * (y0 - mu.phi0)
    } else {
      c(w) * dtheta.deta(mu.phi0, link = .lpobs.0 , earg = .epobs.0 ) *
             dl.dphimat
    }

    ans <- cbind(temp3,
                 c(w) * dl.dlambda * dlambda.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs.0 = lpobs.0, .llambda = llambda,
            .epobs.0 = epobs.0, .elambda = elambda ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1 * NOS)



    temp5 <- expm1(lambda)
    ned2l.dlambda2 <- (1 - phimat) * (temp5 + 1) *
                      (1 / lambda - 1 / temp5) / temp5
    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dlambda2 * dlambda.deta^2


    tmp100 <- mu.phi0 * (1 - mu.phi0)
    tmp200 <- if ( .lpobs.0 == "logit" && is.empty.list( .epobs.0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi0, link = .lpobs.0 , earg = .epobs.0 )^2)
    }


  if (FALSE)
    for (ii in 1:NOS) {
      index200 <- abs(tmp200[, ii]) < .Machine$double.eps
      if (any(index200)) {
        tmp200[index200, ii] <- 10.0 * .Machine$double.eps^(3/4)
      }
    }


    wz[, 1:NOS] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]



    wz
  }), list( .lpobs.0 = lpobs.0,
            .epobs.0 = epobs.0 ))))
}  # End of zapoisson





 zapoissonff <-
  function(llambda = "loge", lonempobs0 = "logit",
           type.fitted = c("mean", "lambda", "pobs0", "onempobs0"),
           imethod = 1,
           ilambda = NULL, ionempobs0 = NULL, ishrinkage = 0.95,
           probs.y = 0.35,
           zero = "onempobs0") {



  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs0", "onempobs0"))[1]


  new("vglmff",
  blurb = c("Zero-altered Poisson ",
            "(Bernoulli and positive-Poisson conditional model)\n\n",
            "Links:    ",
            namesof("lambda",     llambda,    earg = elambda,
                    tag = FALSE), ", ",
            namesof("onempobs0",  lonempobs0, earg = eonempobs0,
                    tag = FALSE), "\n",
            "Mean:     onempobs0 * lambda / (1 - exp(-lambda))"),

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
         parameters.names = c("lambda", "onempobs0"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.nonnegative.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    mynames1 <- param.names("lambda",    ncoly)
    mynames2 <- param.names("onempobs0", ncoly)
    predictors.names <-
        c(namesof(mynames1, .llambda,     earg = .elambda    , tag = FALSE),
          namesof(mynames2, .lonempobs0 , earg = .eonempobs0 , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                             imu = .ilambda,
                             ishrinkage = .ishrinkage,
                             pos.only = TRUE,
                             probs.y = .probs.y )

      etastart <-
        cbind(theta2eta(lambda.init, .llambda , earg = .elambda ),
              theta2eta(1 - (0.5 + w * y0) / (1 + w),
                        .lonempobs0 , earg = .eonempobs0 ))

      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lonempobs0 = lonempobs0, .llambda = llambda,
            .eonempobs0 = eonempobs0, .elambda = elambda,
                                      .ilambda = ilambda,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .type.fitted = type.fitted,
            .imethod = imethod ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "lambda", "pobs0", "onempobs0"))[1]

    M1 <- 2
    NOS <- ncol(eta) / M1

    lambda    <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                 .llambda    , earg = .elambda    ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))


    ans <- switch(type.fitted,
                  "mean"      = onempobs0 * lambda / (-expm1(-lambda)),
                  "lambda"    =    lambda,
                  "pobs0"     = 1 - onempobs0,  # P(Y=0)
                  "onempobs0" =     onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempobs0 = lonempobs0, .llambda = llambda,
           .eonempobs0 = eonempobs0, .elambda = elambda ))),
  last = eval(substitute(expression({
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

    temp.names <- c(rep( .llambda    , len = NOS),
                    rep( .lonempobs0 , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names
    names(misc$link) <-
      c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    misc$earg <- vector("list", M1 * NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .elambda
      misc$earg[[M1*ii  ]] <- .eonempobs0
    }
  }), list( .lonempobs0 = lonempobs0, .llambda = llambda,
            .eonempobs0 = eonempobs0, .elambda = elambda ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    lambda     <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                  .llambda    , earg = .elambda    ))
    onempobs0  <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                  .lonempobs0 , earg = .eonempobs0 ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzapois(x = y, lambda = lambda, pobs0 = 1 - onempobs0,
                       log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempobs0 = lonempobs0, .llambda = llambda,
           .eonempobs0 = eonempobs0, .elambda = elambda ))),
  vfamily = c("zapoissonff"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda    <- eta2theta(eta[, c(TRUE, FALSE)], .llambda    ,
                           earg = .elambda    )
    onempobs0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempobs0 ,
                           earg = .eonempobs0 )
    rzapois(nsim * length(lambda), lambda = lambda, pobs0 = 1 - onempobs0)
  }, list( .lonempobs0 = lonempobs0, .llambda = llambda,
           .eonempobs0 = eonempobs0, .elambda = elambda ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    lambda   <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                .llambda, earg = .elambda ))
    omphimat <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                .lonempobs0, earg = .eonempobs0 ))
    phimat <- 1 - omphimat


    dl.dlambda <- y / lambda + 1 / expm1(-lambda)
    dl.dPHImat <- +1 / (omphimat)  # For y > 0 obsns

    for (spp. in 1:NOS) {
      dl.dPHImat[skip[, spp.], spp.] <- -1 / phimat[skip[, spp.], spp.]
      dl.dlambda[skip[, spp.], spp.] <-  0
    }
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )
    mu.phi0 <- omphimat

    temp3 <- if ( FALSE && .lonempobs0 == "logit") {
    } else {
      c(w) * dtheta.deta(mu.phi0, link = .lonempobs0 , earg = .eonempobs0 ) *
            dl.dPHImat
    }

    ans <- cbind(c(w) * dl.dlambda * dlambda.deta,
                 temp3)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lonempobs0 = lonempobs0, .llambda = llambda,
            .eonempobs0 = eonempobs0, .elambda = elambda ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1 * NOS)

    temp5 <- expm1(lambda)

    ned2l.dlambda2 <- (1 - phimat) * (temp5 + 1) *
                      (1 / lambda - 1 / temp5) / temp5


    wz[, 0 * NOS + (1:NOS)] <- c(w) * ned2l.dlambda2 * dlambda.deta^2


    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lonempobs0 == "logit" && is.empty.list( .eonempobs0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi0, link = .lonempobs0, earg = .eonempobs0)^2)
    }


    wz[, 1 * NOS + (1:NOS)] <-  tmp200

    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]



    wz
  }), list( .lonempobs0 = lonempobs0,
            .eonempobs0 = eonempobs0 ))))
}  # End of zapoissonff







zanegbinomial.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 zanegbinomial <-
  function(
           zero = "size",
           type.fitted = c("mean", "munb", "pobs0"),
           nsimEIM = 500,
           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.support = 4000,  # 20160127; I have changed this
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           lpobs0 = "logit", lmunb = "loge", lsize = "loge",
           imethod = 1,
           ipobs0 = NULL,
           imunb = NULL,
           probs.y = 0.35,
           ishrinkage = 0.95,
           isize = NULL,

           gsize.mux = exp((-12:6)/2)) {





  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  if (length(ipobs0) && (!is.Numeric(ipobs0, positive = TRUE) ||
     max(ipobs0) >= 1))
    stop("If given, argument 'ipobs0' must contain values in (0,1) only")

  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("If given, argument 'isize' must contain positive values only")

  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "munb", "pobs0"))[1]

  ipobs0.small <- 1/64  # A number easily represented exactly

  new("vglmff",
  blurb = c("Zero-altered negative binomial (Bernoulli and\n",
            "positive-negative binomial conditional model)\n\n",
            "Links:    ",
            namesof("pobs0", lpobs0, earg = epobs0, tag = FALSE), ", ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), "\n",
            "Mean:     (1 - pobs0) * munb / (1 - (size / (size + ",
                                                  "munb))^size)"),
  constraints = eval(substitute(expression({

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         imethod = .imethod ,
         multipleResponses = TRUE,
         parameters.names = c("pobs0", "munb", "size"),
         nsimEIM = .nsimEIM ,
         eps.trig = .eps.trig ,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero, .imethod = imethod,
           .nsimEIM = nsimEIM, .eps.trig = eps.trig,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.nonnegative.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    M <- M1 * ncoly

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    mynames1 <- param.names("pobs0", NOS)
    mynames2 <- param.names("munb",  NOS)
    mynames3 <- param.names("size",  NOS)
    predictors.names <-
        c(namesof(mynames1, .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof(mynames2, .lmunb  , earg = .emunb  , tag = FALSE),
          namesof(mynames3, .lsize  , earg = .esize  , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)


    if (!length(etastart)) {
      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                           imu = .imunb , ishrinkage = .ishrinkage ,
                           pos.only = TRUE,
                           probs.y = .probs.y )



      pobs0.init <- matrix(if (length( .ipobs0 )) .ipobs0 else -1,
                           nrow = n, ncol = NOS, byrow = TRUE)
      for (jay in 1:NOS) {
        if (any(pobs0.init[, jay] < 0)) {
          index.y0 <- (y[, jay] < 0.5)
          pobs0.init[, jay] <- max(min(mean(index.y0), 1 - .ipobs0.small ),
                                   .ipobs0.small )
        }
      }


      if ( is.Numeric( .isize )) {
        size.init <- matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
         munb <- extraargs
         sum(c(w) * dposnegbin(y, munb = munb, size = kmat, log = TRUE))
        }
        size.init <- matrix(0, nrow = n, ncol = NOS) 
        for (jay in 1:NOS) {
          size.grid <- .gsize.mux * mean(munb.init[, jay])
          TFvec <- (y[, jay] > 0)
          size.init[, jay] <-
            grid.search(size.grid, objfun = posnegbinomial.Loglikfun,
                        y = y[TFvec, jay],  # x = x[TFvec, ],
                        w = w[TFvec, jay],
                        extraargs = munb.init[TFvec, jay])
        }
      }

      etastart <- cbind(theta2eta(pobs0.init, .lpobs0 , earg = .epobs0 ),
                        theta2eta(munb.init,  .lmunb  , earg = .emunb  ),
                        theta2eta(size.init,  .lsize  , earg = .esize  ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }  # End of if (!length(etastart))


  }), list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
            .epobs0 = epobs0, .emunb = emunb, .esize = esize,
            .ipobs0 = ipobs0,                 .isize = isize,
            .ipobs0.small = ipobs0.small,
                              .imunb = imunb, .gsize.mux = gsize.mux,
            .imethod = imethod, .ishrinkage = ishrinkage,
            .type.fitted = type.fitted, .probs.y = probs.y ))),


  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "munb", "pobs0"))[1]

    M1 <- 3
    NOS <- ncol(eta) / M1
    phi0 <- eta2theta(eta[, M1*(1:NOS)-2], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, M1*(1:NOS)-1], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)  ], .lsize  , earg = .esize  )



    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
    prob0  <- tempk^kmat  # p(0) from negative binomial
    oneminusf0  <- 1 - prob0

    smallval <- 1e-3  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
      prob0[big.size]  <- exp(-munb[big.size])  # The limit as kmat --> Inf
      oneminusf0[big.size] <- -expm1(-munb[big.size])
    }

    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) * munb / oneminusf0,
                  "munb"      = munb,
                  "pobs0"     = phi0)  # P(Y=0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpobs0 = lpobs0, .lsize = lsize, .lmunb = lmunb,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),


  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lpobs0 , length = NOS),
        rep( .lmunb  , length = NOS),
        rep( .lsize  , length = NOS))[interleave.VGAM(M1*NOS, M1 = M1)]
    temp.names <- c(mynames1,
                    mynames2,
                    mynames3)[interleave.VGAM(M1*NOS, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii - 2]] <- .epobs0
      misc$earg[[M1*ii - 1]] <- .emunb
      misc$earg[[M1*ii    ]] <- .esize
    }

    misc$nsimEIM <- .nsimEIM
    misc$ipobs0  <- .ipobs0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE
  }), list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
            .epobs0 = epobs0, .emunb = emunb, .esize = esize,
            .ipobs0 = ipobs0,                 .isize = isize,
            .nsimEIM = nsimEIM ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 3
    NOS <- ncol(eta) / M1
    phi0 <- eta2theta(eta[, M1*(1:NOS)-2], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, M1*(1:NOS)-1], .lmunb  , earg = .emunb  )
    size <- eta2theta(eta[, M1*(1:NOS)  ], .lsize  , earg = .esize  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzanegbin(x = y, pobs0 = phi0, munb = munb, size = size,
                         log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zanegbinomial"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    phi0 <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lsize  , earg = .esize  )
    rzanegbin(nsim * length(munb),
              pobs0 = phi0, munb = munb, size = kmat)
  }, list( .lpobs0 = lpobs0, .lmunb = lmunb, .lsize = lsize,
           .epobs0 = epobs0, .emunb = emunb, .esize = esize ))),





  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- ncol(eta) / M1
    y0 <- extra$y0

    phi0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                      .lpobs0 , earg = .epobs0 )
    munb <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                      .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                      .lsize , earg = .esize )
    skip <- extra$skip.these


    dphi0.deta <- dtheta.deta(phi0, .lpobs0 , earg = .epobs0 )
    dmunb.deta <- dtheta.deta(munb, .lmunb  , earg = .emunb  )
    dsize.deta <- dtheta.deta(kmat, .lsize  , earg = .esize  )



    smallval <- 1e-3  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
        warning("parameter 'size' has very large values; ",
                "try fitting a zero-altered Poisson ",
                "model instead")
        kmat[big.size] <- munb[big.size] / smallval
    }



    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    AA16 <- tempm + log(tempk)
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- prob0 * AA16
    df02.dmunb2 <- prob0 * tempk * (1 + 1/kmat) / (1 + munb/kmat)
    df02.dkmat2 <- prob0 * ((tempm^2) / kmat + AA16^2)
    df02.dkmat.dmunb <- -prob0 * (tempm/kmat + AA16) / (1 + munb/kmat)




    if (any(big.size)) {
      prob0[big.size]  <- exp(-munb[big.size])  # The limit as kmat --> Inf
      oneminusf0[big.size] <- -expm1(-munb[big.size])
      df0.dmunb[big.size] <- -tempk[big.size] * prob0[big.size]
      df0.dkmat[big.size] <-  prob0[big.size] * AA16[big.size]
      df02.dmunb2[big.size] <- prob0[big.size] * tempk[big.size] *
        (1 + 1/kmat[big.size]) / (1 + smallval)
      df02.dkmat2[big.size] <- prob0[big.size] *
        ((tempm[big.size])^2 / kmat[big.size] + AA16[big.size]^2)
      df02.dkmat.dmunb[big.size] <- -prob0[big.size] *
        (tempm[big.size]/kmat[big.size] + AA16[big.size]) / (1 + smallval)
    }


    mymu <- munb / oneminusf0  # E(Y) of Pos-NBD



    dl.dphi0 <- -1 / (1 - phi0)
    dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat) +
                df0.dmunb / oneminusf0
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y - munb) / (munb + kmat) + log(tempk) +
                df0.dkmat / oneminusf0


    if (any(big.size)) {
      dl.dsize[big.size] <- 1e-8  # A small number
    }


    dl.dphi0[y == 0] <-  1 / phi0[y == 0]  # Do it in one line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dsize[skip[, spp.], spp.] <-
      dl.dmunb[skip[, spp.], spp.] <- 0
    }

    dl.deta23 <- c(w) * cbind(dl.dmunb * dmunb.deta,
                              dl.dsize * dsize.deta)


    dl.deta1 <- if ( .lpobs0 == "logit") {
      c(w) * (y0 - phi0)
    } else {
      c(w) * dl.dphi0 * dphi0.deta
    }


    ans <- cbind(dl.deta1, dl.deta23)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs0 = lpobs0 , .lmunb = lmunb , .lsize = lsize ,
            .epobs0 = epobs0 , .emunb = emunb , .esize = esize  ))),



  weight = eval(substitute(expression({
    wz <- matrix(0, n, M + M-1)  # tridiagonal


    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB


    mu.phi0 <- phi0  # pobs0  # phi0
    tmp100 <- mu.phi0 * (1 - mu.phi0)
    wz[, (1:NOS)*M1 - 2] <-
    if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(mu.phi0, link = .lpobs0 , earg = .epobs0 )^2)
    }


    ned2l.dmunb2 <- mymu / munb^2 -
        ((1 + mymu/kmat) / kmat) / (1 + munb/kmat)^2 -
        df02.dmunb2 / oneminusf0 -
        (df0.dmunb / oneminusf0)^2
    wz[,     M1*(1:NOS) - 1] <- c(w) * (1 - phi0) *
                                ned2l.dmunb2 * dmunb.deta^2


    ned2l.dmunbsize <- (munb - mymu) / (munb + kmat)^2 -
      df02.dkmat.dmunb / oneminusf0 -
      df0.dmunb * df0.dkmat / oneminusf0^2
    wz[, M + M1*(1:NOS) - 1] <- c(w) * (1 - phi0) *
                                ned2l.dmunbsize * dmunb.deta * dsize.deta






    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 1
      Q.maxs <-      qposnegbin(p     = eff.p[2] ,
                                munb = munb[, jay],
                                size  = kmat[, jay]) + 10


      eps.trig <- .eps.trig
      Q.MAXS <-      pmax(10, ceiling(1 / sqrt(eps.trig)))
      Q.maxs <- pmin(Q.maxs, Q.MAXS)



      ind1 <- if (max.chunk.MB > 0) (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]

          wz[sind2, M1*jay] <-
            EIM.posNB.specialp(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.max = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only)
  if (FALSE)
          wz2[sind2, M1*jay] <-
            EIM.posNB.speciald(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.min       = min(Q.mins2[sind2]),
                               y.max       = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only)  # *



          if (any(eim.kk.TF <-       wz[sind2, M1*jay] <= 0 |
                               is.na(wz[sind2, M1*jay]))) {
            ind2[sind2[eim.kk.TF], jay] <- FALSE
          }
          
          
          lwr.ptr <- upr.ptr + 1
        }  # while
      }  # if
    }  # end of for (jay in 1:NOS)

    






    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        kkvec <- kmat[ii.TF, jay]
        muvec <- munb[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rzanegbin(sum(ii.TF), munb = muvec, size = kkvec,
                            pobs0 = phi0[ii.TF, jay])
          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) -
                   (ysim - muvec) / (muvec + kkvec) +
                   log1p(-muvec / (kkvec + muvec)) +
                   df0.dkmat[ii.TF, jay] / oneminusf0[ii.TF, jay]

          dl.dk[ysim == 0] <- 0

          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay] <- ned2l.dk2  # * (dsize.deta[ii.TF, jay])^2
      }
    }  # jay



    wz[, M1*(1:NOS)    ] <- wz[, M1*(1:NOS)    ] * dsize.deta^2



    save.weights <- !all(ind2)




    wz[,     M1*(1:NOS)    ] <- c(w) * (1 - phi0) *
                                wz[,     M1*(1:NOS)    ]



    wz
  }), list( .lpobs0 = lpobs0,
            .epobs0 = epobs0,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB = max.chunk.MB,
            .nsimEIM = nsimEIM ))))
}  # End of zanegbinomial()





zanegbinomialff.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 zanegbinomialff <-
  function(
           lmunb = "loge", lsize = "loge", lonempobs0 = "logit",
           type.fitted = c("mean", "munb", "pobs0", "onempobs0"),
           isize = NULL, ionempobs0 = NULL,
           zero = c("size", "onempobs0"),

           probs.y = 0.35,
           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.support = 4000,  # 20160127; I have changed this
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           gsize.mux = exp((-12:6)/2),

           imethod = 1,
           imunb = NULL,
           nsimEIM = 500,
           ishrinkage = 0.95) {




  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  if (length(ionempobs0) && (!is.Numeric(ionempobs0, positive = TRUE) ||
     max(ionempobs0) >= 1))
    stop("If given, argument 'ionempobs0' must contain values in (0,1) only")

  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("If given, argument 'isize' must contain positive values only")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")


  ipobs0.small <- 1/64  # A number easily represented exactly

  type.fitted <- match.arg(type.fitted,
                           c("mean", "munb", "pobs0", "onempobs0"))[1]


  new("vglmff",
  blurb = c("Zero-altered negative binomial (Bernoulli and\n",
            "positive-negative binomial conditional model)\n\n",
            "Links:    ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), ", ",
            namesof("onempobs0", lonempobs0, earg = eonempobs0,
                    tag = FALSE), "\n",
            "Mean:     onempobs0 * munb / (1 - (size / (size + ",
                                                 "munb))^size)"),
  constraints = eval(substitute(expression({

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         nsimEIM = .nsimEIM ,
         parameters.names = c("munb", "size", "onempobs0"),
         eps.trig = .eps.trig ,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM, .eps.trig = eps.trig,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.nonnegative.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    M <- M1 * ncoly

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    mynames1 <- param.names("munb",      NOS)
    mynames2 <- param.names("size",      NOS)
    mynames3 <- param.names("onempobs0", NOS)
    predictors.names <-
        c(namesof(mynames1, .lmunb  , earg = .emunb  , tag = FALSE),
          namesof(mynames2, .lsize  , earg = .esize  , tag = FALSE),
          namesof(mynames3, .lonempobs0 , earg = .eonempobs0 ,
                  tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]


    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)


    if (!length(etastart)) {
      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                           imu = .imunb , ishrinkage = .ishrinkage ,
                           pos.only = TRUE,
                           probs.y = .probs.y )


      pobs0.init <- matrix(if (length( .ionempobs0 )) 1 - .ionempobs0 else -1,
                           nrow = n, ncol = NOS, byrow = TRUE)
      for (jay in 1:NOS) {
        if (any(pobs0.init[, jay] < 0)) {
          index.y0 <- y[, jay] < 0.5
          pobs0.init[, jay] <- max(min(mean(index.y0), 1 - .ipobs0.small ),
                                   .ipobs0.small )
        }
      }


      if ( is.Numeric( .isize )) {
        size.init <- matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
         munb <- extraargs
         sum(c(w) * dposnegbin(x = y, munb = munb, size = kmat,
                               log = TRUE))
        }
        size.init <- matrix(0, nrow = n, ncol = NOS) 
        for (jay in 1:NOS) {
          size.grid <- .gsize.mux * mean(munb.init[, jay])
          TFvec <- (y[, jay] > 0)
          size.init[, jay] <-
            grid.search(size.grid, objfun = posnegbinomial.Loglikfun,
                        y = y[TFvec, jay],  # x = x[index.posy, ],
                        w = w[TFvec, jay],
                        extraargs = munb.init[TFvec, jay])
        }
      }

      etastart <-
        cbind(theta2eta(munb.init ,     .lmunb      , earg = .emunb      ),
              theta2eta(size.init ,     .lsize      , earg = .esize      ),
              theta2eta(1 - pobs0.init, .lonempobs0 , earg = .eonempobs0 ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }  # End of if (!length(etastart))


  }), list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
            .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize,
            .ionempobs0 = ionempobs0, .imunb = imunb, .isize = isize,
                                                      .gsize.mux = gsize.mux,
            .ipobs0.small = ipobs0.small,
            .imethod = imethod, .ishrinkage = ishrinkage,
            .probs.y = probs.y, .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "munb", "pobs0", "onempobs0"))[1]

    M1 <- 3
    NOS <- ncol(eta) / M1
    munb <- eta2theta(eta[, M1*(1:NOS)-2], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)-1], .lsize  , earg = .esize  )
    onempobs0 <- eta2theta(eta[, M1*(1:NOS)  ], .lonempobs0 ,
                           earg = .eonempobs0 )


    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb); NBD p(0)
    prob0  <- tempk^kmat  # p(0) from negative binomial
    oneminusf0  <- 1 - prob0

    smallval <- 1e-3  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
      prob0[big.size]  <- exp(-munb[big.size])  # The limit as kmat --> Inf
      oneminusf0[big.size] <- -expm1(-munb[big.size])
    }


    ans <- switch(type.fitted,
                  "mean"      =    onempobs0 * munb / oneminusf0,
                  "munb"      =    munb,
                  "pobs0"     = 1 - onempobs0,  # P(Y=0)
                  "onempobs0" =     onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))        
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempobs0 = lonempobs0, .lsize = lsize, .lmunb = lmunb,
           .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize ))),
  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lmunb      , length = NOS),
        rep( .lsize      , length = NOS),
        rep( .lonempobs0 , length = NOS))[
        interleave.VGAM(M1*NOS, M1 = M1)]
    temp.names <- c(mynames1,
                    mynames2,
                    mynames3)[interleave.VGAM(M1*NOS, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .emunb
      misc$earg[[M1*ii-1]] <- .esize
      misc$earg[[M1*ii  ]] <- .eonempobs0
    }

    misc$nsimEIM <- .nsimEIM
    misc$imethod <- .imethod
    misc$ionempobs0  <- .ionempobs0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE
  }), list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
            .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize,
            .ionempobs0 = ionempobs0, .isize = isize,
            .nsimEIM = nsimEIM,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 3
    NOS <- ncol(eta) / M1
    munb <- eta2theta(eta[, M1*(1:NOS)-2], .lmunb  , earg = .emunb  )
    kmat <- eta2theta(eta[, M1*(1:NOS)-1], .lsize  , earg = .esize  )
    onempobs0 <- eta2theta(eta[, M1*(1:NOS)  ], .lonempobs0 ,
                           earg = .eonempobs0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzanegbin(x = y, pobs0 = 1 - onempobs0,
                         munb = munb, size = kmat, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
           .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zanegbinomialff"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    munb      <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lmunb , earg = .emunb )
    kmat      <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lsize , earg = .esize )
    onempobs0 <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lonempobs0 ,
                           earg = .eonempobs0 )

    rzanegbin(nsim * length(munb),
              pobs0 = 1 - onempobs0, munb = munb, size = kmat)
  }, list( .lonempobs0 = lonempobs0, .lmunb = lmunb, .lsize = lsize,
           .eonempobs0 = eonempobs0, .emunb = emunb, .esize = esize ))),




  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- ncol(eta) / M1
    y0 <- extra$y0

    munb      <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb      , earg = .emunb )
    kmat      <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize      , earg = .esize )
    onempobs0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempobs0 , earg = .eonempobs0 )
    skip <- extra$skip.these
    phi0 <- 1 - onempobs0

    dmunb.deta      <- dtheta.deta(munb, .lmunb , earg = .emunb )
    dsize.deta      <- dtheta.deta(kmat, .lsize , earg = .esize )
    donempobs0.deta <- dtheta.deta(onempobs0, .lonempobs0 ,
                                   earg = .eonempobs0 )






    smallval <- 1e-3  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
        warning("parameter 'size' has very large values; ",
                "try fitting a zero-altered Poisson ",
                "model instead")
        kmat[big.size] <- munb[big.size] / smallval
    }



    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    AA16 <- tempm + log(tempk)
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- prob0 * AA16
    df02.dmunb2 <- prob0 * tempk * (1 + 1/kmat) / (1 + munb/kmat)
    df02.dkmat2 <- prob0 * ((tempm^2) / kmat + AA16^2)
    df02.dkmat.dmunb <- -prob0 * (tempm/kmat + AA16) / (1 + munb/kmat)



    mymu <- munb / oneminusf0  # E(Y) of Pos-NBD




    dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat) +
                df0.dmunb / oneminusf0
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y - munb) / (munb + kmat) + log(tempk) +
                df0.dkmat / oneminusf0
    dl.donempobs0 <- +1 / (onempobs0)



    if (any(big.size)) {
      dl.dsize[big.size] <- 1e-8  # A small number
    }



    dl.donempobs0[y == 0] <-
      -1 / (1 - onempobs0[y == 0])  # Do it in 1 line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dsize[skip[, spp.], spp.] <-
      dl.dmunb[skip[, spp.], spp.] <- 0
    }

    dl.deta12 <- c(w) * cbind(dl.dmunb * dmunb.deta,
                              dl.dsize * dsize.deta)



    dl.deta3 <- if ( .lonempobs0 == "logit") {
      -c(w) * (y0 - phi0)
    } else {
      -c(w) * dl.donempobs0 * donempobs0.deta
    }



    ans <- cbind(dl.deta12, dl.deta3)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lonempobs0 = lonempobs0 , .lmunb = lmunb , .lsize = lsize ,
            .eonempobs0 = eonempobs0 , .emunb = emunb , .esize = esize  ))),







  weight = eval(substitute(expression({

    wz <- matrix(0, n, M + M-1)  # tridiagonal


    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB



    tmp100 <- onempobs0 * (1 - onempobs0)
    wz[, (1:NOS)*M1    ] <-
    if ( .lonempobs0 == "logit" && is.empty.list( .eonempobs0 )) {
        cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (1 / tmp100) *
            dtheta.deta(onempobs0, link = .lonempobs0 , earg = .eonempobs0 )^2)
    }




    ned2l.dmunb2 <- mymu / munb^2 -
        ((1 + mymu/kmat) / kmat) / (1 + munb/kmat)^2 -
        df02.dmunb2 / oneminusf0 -
        (df0.dmunb / oneminusf0)^2
    wz[,     M1*(1:NOS) - 2] <- c(w) * (1 - phi0) *
                                ned2l.dmunb2 * dmunb.deta^2


    ned2l.dmunbsize <- (munb - mymu) / (munb + kmat)^2 -
      df02.dkmat.dmunb / oneminusf0 -
      df0.dmunb * df0.dkmat / oneminusf0^2
    wz[, M + M1*(1:NOS) - 2] <- c(w) * (1 - phi0) *
                                ned2l.dmunbsize * dmunb.deta * dsize.deta







    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 1
      Q.maxs <-      qposnegbin(p     = eff.p[2] ,
                                munb = munb[, jay],
                                size  = kmat[, jay]) + 10



      eps.trig <- .eps.trig
      Q.MAXS <-      pmax(10, ceiling(1 / sqrt(eps.trig)))
      Q.maxs <- pmin(Q.maxs, Q.MAXS)




      
      ind1 <- if (max.chunk.MB > 0) (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]

          wz[sind2, M1*jay - 1] <-
            EIM.posNB.specialp(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.max = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only)
  if (FALSE)
          wz2[sind2, M1*jay - 1] <-
            EIM.posNB.speciald(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.min       = min(Q.mins2[sind2]),
                               y.max       = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only) # *



          if (any(eim.kk.TF <-       wz[sind2, M1*jay - 1] <= 0 |
                               is.na(wz[sind2, M1*jay - 1]))) {
            ind2[sind2[eim.kk.TF], jay] <- FALSE
          }
          
          
          lwr.ptr <- upr.ptr + 1
        }  # while
      }  # if
    }  # end of for (jay in 1:NOS)














    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        kkvec <- kmat[ii.TF, jay]
        muvec <- munb[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rzanegbin(sum(ii.TF), munb = muvec, size = kkvec,
                            pobs0 = phi0[ii.TF, jay])
          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) -
                   (ysim - muvec) / (muvec + kkvec) +
                   log1p(-muvec / (kkvec + muvec)) +
                   df0.dkmat[ii.TF, jay] / oneminusf0[ii.TF, jay]

          dl.dk[ysim == 0] <- 0

          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay - 1] <- ned2l.dk2  # * (dsize.deta[ii.TF, jay])^2
      }
    }  # jay




    wz[, M1*(1:NOS) - 1] <- wz[, M1*(1:NOS) - 1] * dsize.deta^2






    save.weights <- !all(ind2)




    wz[,     M1*(1:NOS) - 1] <- c(w) * (1 - phi0) *
                                wz[,     M1*(1:NOS) - 1]



    wz
  }), list( .lonempobs0 = lonempobs0,
            .eonempobs0 = eonempobs0,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB = max.chunk.MB,
            .nsimEIM = nsimEIM ))))
}  # End of zanegbinomialff()











 zipoisson <-
  function(lpstr0 = "logit", llambda = "loge",
           type.fitted = c("mean", "lambda", "pobs0", "pstr0", "onempstr0"),
           ipstr0 = NULL,    ilambda = NULL,
           gpstr0 = NULL,  # (1:9) / 10,
           imethod = 1,
           ishrinkage = 0.95, probs.y = 0.35,
           zero = NULL) {
  ipstr00 <- ipstr0
  gpstr00 <- gpstr0
  ipstr0.small <- 1/64  # A number easily represented exactly


  lpstr0 <- as.list(substitute(lpstr0))
  epstr00 <- link2list(lpstr0)
  lpstr00 <- attr(epstr00, "function.name")

  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")



  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "pobs0", "pstr0", "onempstr0"))[1]


  if (length(ipstr00))
    if (!is.Numeric(ipstr00, positive = TRUE) ||
        any(ipstr00 >= 1))
      stop("argument 'ipstr0' values must be inside the interval (0,1)")
  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("argument 'ilambda' values must be positive")


  new("vglmff",
  blurb = c("Zero-inflated Poisson\n\n",
            "Links:    ",
            namesof("pstr0",  lpstr00, earg = epstr00 ), ", ",
            namesof("lambda", llambda, earg = elambda ), "\n",
            "Mean:     (1 - pstr0) * lambda"),

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
         parameters.names = c("pstr0", "lambda"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$M1 <- M1
    extra$dimnamesy <- dimnames(y)
    M <- M1 * ncoly
    extra$type.fitted      <- .type.fitted

    mynames1 <- param.names("pstr0",  ncoly)
    mynames2 <- param.names("lambda", ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpstr00 , earg = .epstr00 , tag = FALSE),
          namesof(mynames2, .llambda , earg = .elambda , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]



    if (!length(etastart)) {


      matL <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                      imu = .ilambda , ishrinkage = .ishrinkage ,
                      pos.only = TRUE,
                      probs.y = .probs.y )


      matP <- matrix(if (length( .ipstr00 )) .ipstr00 else 0,
                     n, ncoly, byrow = TRUE)
      phi.grid <- .gpstr00  # seq(0.02, 0.98, len = 21)
      ipstr0.small <- .ipstr0.small  # A number easily represented exactly

      if (!length( .ipstr00 ))
      for (jay in 1:ncoly) {

        zipois.Loglikfun <- function(phival, y, x, w, extraargs) {
          sum(c(w) * dzipois(x = y, pstr0 = phival,
                             lambda = extraargs$lambda, log = TRUE))
        }
        Phi.init <- if (length(phi.grid)) {
          grid.search(phi.grid, objfun = zipois.Loglikfun,
                      y = y[, jay], x = x, w = w[, jay],
                      extraargs = list(lambda = matL[, jay]))
        } else {
          pmax(ipstr0.small,
               weighted.mean(y[, jay] == 0, w[, jay]) -
               dpois(0, matL[, jay]))
        }
        if (mean(Phi.init == ipstr0.small) > 0.95)
          warning("from the initial values only, the data appears to ",
                  "have little or no 0-inflation")
        matP[, jay] <- Phi.init
      }  # for (jay)

      etastart <- cbind(theta2eta(matP, .lpstr00 , earg = .epstr00 ),
                        theta2eta(matL, .llambda , earg = .elambda ))[,
                        interleave.VGAM(M, M1 = M1)]
      mustart <- NULL  # Since etastart has been computed.
    }  # End of !length(etastart)
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda,
            .ipstr00 = ipstr00, .ilambda = ilambda,
            .gpstr00 = gpstr00,
            .imethod = imethod, .probs.y = probs.y,
            .ipstr0.small = ipstr0.small,
            .type.fitted = type.fitted,
            .ishrinkage = ishrinkage ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs0", "pstr0", "onempstr0"))[1]

    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )

    
    ans <- switch(type.fitted,
                  "mean"      = (1 - phimat) * lambda,
                  "lambda"    = lambda,
                  "pobs0"     = phimat + (1-phimat)*exp(-lambda),  # P(Y=0)
                  "pstr0"     =     phimat,
                  "onempstr0" = 1 - phimat)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans)) 
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda,
           .type.fitted = type.fitted
         ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .lpstr00 , length = ncoly),
        rep( .llambda , length = ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .epstr00
      misc$earg[[M1*ii  ]] <- .elambda
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

      misc$pobs0 <- phimat + (1 - phimat) * exp(-lambda)  # P(Y=0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pobs0) <- dimnames(y)

      misc$pstr0 <- phimat
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pstr0) <- dimnames(y)
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dzipois(x = y, pstr0 = phimat, lambda = lambda,
                                log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda ))),
  vfamily = c("zipoisson"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    phimat <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr00 , earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE)], .llambda , earg = .elambda )
    rzipois(nsim * length(lambda), lambda = lambda, pstr0 = phimat)
  }, list( .lpstr00 = lpstr00, .llambda = llambda,
           .epstr00 = epstr00, .elambda = elambda ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    phimat <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE], .lpstr00 ,
                        earg = .epstr00 )
    lambda <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE], .llambda ,
                        earg = .elambda )

    prob0 <- exp(-lambda)
    pobs0 <- phimat + (1 - phimat) * prob0
    index0 <- as.matrix(y == 0)

    dl.dphimat <- -expm1(-lambda) / pobs0
    dl.dphimat[!index0] <- -1 / (1 - phimat[!index0])

    dl.dlambda <- -(1 - phimat) * exp(-lambda) / pobs0
    dl.dlambda[!index0] <- (y[!index0] - lambda[!index0]) / lambda[!index0]

    dphimat.deta <- dtheta.deta(phimat, .lpstr00 , earg = .epstr00 )
    dlambda.deta <- dtheta.deta(lambda, .llambda , earg = .elambda )

    ans <- c(w) * cbind(dl.dphimat * dphimat.deta,
                        dl.dlambda * dlambda.deta)
    ans <- ans[, interleave.VGAM(M, M1 = M1)]


    if ( .llambda == "loge" && is.empty.list( .elambda ) &&
       any(lambda[!index0] < .Machine$double.eps)) {
      for (spp. in 1:(M / M1)) {
        ans[!index0[, spp.], M1 * spp.] <-
          w[!index0[, spp.]] *
         (y[!index0[, spp.], spp.] - lambda[!index0[, spp.], spp.])
      }
    }

    ans
  }), list( .lpstr00 = lpstr00, .llambda = llambda,
            .epstr00 = epstr00, .elambda = elambda ))),
  weight = eval(substitute(expression({

    ned2l.dphimat2 <- -expm1(-lambda) / ((1 - phimat) * pobs0)
    ned2l.dphimatlambda <- -exp(-lambda) / pobs0
    ned2l.dlambda2 <- (1 - phimat) / lambda -
                      phimat * (1 - phimat) * exp(-lambda) / pobs0




    wz <- array(c(c(w) * ned2l.dphimat2 * dphimat.deta^2,
                  c(w) * ned2l.dlambda2 * dlambda.deta^2,
                  c(w) * ned2l.dphimatlambda * dphimat.deta * dlambda.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)





    wz
  }), list( .llambda = llambda, .elambda = elambda ))))
}  # zipoisson









 zibinomial <-
  function(lpstr0 = "logit", lprob = "logit",
           type.fitted = c("mean", "prob", "pobs0", "pstr0", "onempstr0"),
           ipstr0 = NULL,
           zero = NULL,  # 20130917; was originally zero = 1,
           multiple.responses = FALSE, imethod = 1) {
  if (as.logical(multiple.responses))
    stop("argument 'multiple.responses' must be FALSE")

  lpstr0 <- as.list(substitute(lpstr0))
  epstr0 <- link2list(lpstr0)
  lpstr0 <- attr(epstr0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]


  if (is.Numeric(ipstr0))
    if (!is.Numeric(ipstr0, positive = TRUE) || any(ipstr0 >= 1))
      stop("'ipstr0' values must be inside the interval (0,1)")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Zero-inflated binomial\n\n",
            "Links:    ",
            namesof("pstr0", lpstr0, earg = epstr0), ", ",
            namesof("prob" , lprob , earg = eprob ), "\n",
            "Mean:     (1 - pstr0) * prob"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         type.fitted  = .type.fitted ,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("pstr0", "prob"),
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


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
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }


    if ( .imethod == 1)
      mustart <- (mustart + y) / 2


    extra$type.fitted <- .type.fitted
    extra$dimnamesy   <- dimnames(y)




    predictors.names <-
        c(namesof("pstr0", .lpstr0 , earg = .epstr0 , tag = FALSE),
          namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE))


    extra$w <- w  # Needed for @linkinv
    phi.init <- if (length( .ipstr0 )) .ipstr0 else {
        prob0.est <- sum(w[y == 0]) / sum(w)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^w) / (1 - (1 - mustart)^w)
        } else {
          prob0.est
        }
    }

    phi.init[phi.init <= -0.10] <- 0.10  # Lots of sample variation
    phi.init[phi.init <=  0.05] <- 0.15  # Last resort
    phi.init[phi.init >=  0.80] <- 0.80  # Last resort

    if ( length(mustart) && !length(etastart))
      mustart <- cbind(rep(phi.init, len = n),
                       mustart)  # 1st coln not a real mu
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob,
            .ipstr0 = ipstr0,
            .type.fitted = type.fitted,          
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pstr0 <- eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )


    orig.w <- if (length(tmp3 <- extra$orig.w)) tmp3 else
              rep(1, len = nrow(eta))
    priorw <- extra$w
    nvec <- priorw / orig.w


    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = (1 - pstr0) * mubin,
                  "prob"      = mubin,
                  "pobs0"     = pstr0 + (1-pstr0)*(1-mubin)^nvec,  # P(Y=0)
                  "pstr0"     =     pstr0,
                  "onempstr0" = 1 - pstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    misc$link <-    c("pstr0" = .lpstr0 , "prob" = .lprob )

    misc$earg <- list("pstr0" = .epstr0 , "prob" = .eprob )

    misc$imethod <- .imethod


  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    cbind(theta2eta(mu[, 1], .lpstr0 , earg = .epstr0 ),
          theta2eta(mu[, 2], .lprob  , earg = .eprob  ))
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr0 <- eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        dzibinom(x = round(w * y), size = w, prob = mubin,
                 log = TRUE, pstr0 = pstr0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr0 = lpstr0, .lprob = lprob,
           .epstr0 = epstr0, .eprob = eprob ))),
  vfamily = c("zibinomial"),
  deriv = eval(substitute(expression({
    phi   <- eta2theta(eta[, 1], .lpstr0 , earg = .epstr0 )
    mubin <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )

    prob0 <- (1 - mubin)^w  # Actually q^w
    pobs0 <- phi + (1 - phi) * prob0
    index <- (y == 0)
    dl.dphi <- (1 - prob0) / pobs0
    dl.dphi[!index] <- -1 / (1 - phi[!index])

    dl.dmubin <- -w * (1 - phi) * (1 - mubin)^(w - 1) / pobs0
    dl.dmubin[!index] <- w[!index] *
        (    y[!index]  /      mubin[!index]   -
        (1 - y[!index]) / (1 - mubin[!index]))

    dphi.deta   <- dtheta.deta(phi,   .lpstr0 , earg = .epstr0 )
    dmubin.deta <- dtheta.deta(mubin, .lprob  , earg = .eprob  )

    ans <- cbind(dl.dphi   * dphi.deta,
                 dl.dmubin * dmubin.deta)

      if ( .lprob == "logit") {
        ans[!index, 2] <- w[!index] * (y[!index] - mubin[!index])
      }

      ans
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, nrow = n, ncol = dimm(M))



    ned2l.dphi2 <- (1 - prob0) / ((1 - phi) * pobs0)


    ned2l.dphimubin <- -w * ((1 - mubin)^(w - 1)) / pobs0







    ned2l.dmubin2 <- (w * (1 - phi) / (mubin * (1 - mubin)^2)) *
                     (1 - mubin - w * mubin *
                     (1 - mubin)^w * phi / pobs0)





    wz[,iam(1, 1, M)] <- ned2l.dphi2     * dphi.deta^2
    wz[,iam(2, 2, M)] <- ned2l.dmubin2   * dmubin.deta^2
    wz[,iam(1, 2, M)] <- ned2l.dphimubin * dphi.deta * dmubin.deta
    if (TRUE) {
      ind6 <- (wz[, iam(2, 2, M)] < .Machine$double.eps)
      if (any(ind6))
        wz[ind6, iam(2, 2, M)] <- .Machine$double.eps
    }
    wz
  }), list( .lpstr0 = lpstr0, .lprob = lprob,
            .epstr0 = epstr0, .eprob = eprob ))))
}






 zibinomialff <-
  function(lprob = "logit", lonempstr0 = "logit",
           type.fitted = c("mean", "prob", "pobs0", "pstr0", "onempstr0"),
           ionempstr0 = NULL,
           zero = "onempstr0",
           multiple.responses = FALSE, imethod = 1) {






  if (as.logical(multiple.responses))
    stop("argument 'multiple.responses' must be FALSE")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")

  type.fitted <- match.arg(type.fitted,
                   c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]


  if (is.Numeric(ionempstr0))
    if (!is.Numeric(ionempstr0, positive = TRUE) || any(ionempstr0 >= 1))
      stop("'ionempstr0' values must be inside the interval (0,1)")
  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 2)
    stop("argument 'imethod' must be 1 or 2")



  new("vglmff",
  blurb = c("Zero-inflated binomial\n\n",
            "Links:    ",
            namesof("prob" ,     lprob     , earg = eprob     ), ", ",
            namesof("onempstr0", lonempstr0, earg = eonempstr0), "\n",
            "Mean:     onempstr0 * prob"),
  constraints = eval(substitute(expression({
   constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("prob", "onempstr0"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),
      
  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


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
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }


    if ( .imethod == 1)
      mustart <- (mustart + y) / 2


    extra$type.fitted <- .type.fitted
    extra$dimnamesy   <- dimnames(y)




    predictors.names <-
        c(namesof("prob"     , .lprob      , earg = .eprob      , tag = FALSE),
          namesof("onempstr0", .lonempstr0 , earg = .eonempstr0 , tag = FALSE))


    extra$w <- w  # Needed for @linkinv
    onemphi.init <- if (length( .ionempstr0 )) .ionempstr0 else {
        prob0.est <- sum(w[y == 0]) / sum(w)
        if ( .imethod == 1) {
          1 - (prob0.est - (1 - mustart)^w) / (1 - (1 - mustart)^w)
        } else {
          1 - prob0.est
        }
    }

    onemphi.init[onemphi.init <= -0.10] <- 0.10  # Lots of sample variation
    onemphi.init[onemphi.init <=  0.05] <- 0.15  # Last resort
    onemphi.init[onemphi.init >=  0.80] <- 0.80  # Last resort

    if ( length(mustart) && !length(etastart))
      mustart <- cbind(mustart,
                       rep(onemphi.init, len = n))  # 1st coln not a real mu

  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob,
            .ionempstr0 = ionempstr0,
            .type.fitted = type.fitted,          
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    mubin     <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempstr0 <- eta2theta(eta[, 2], .lonempstr0 , earg = .eonempstr0 )


    orig.w <- if (length(tmp3 <- extra$orig.w)) tmp3 else
              rep(1, len = nrow(eta))
    priorw <- extra$w
    nvec <- priorw / orig.w


    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = (onempstr0) * mubin,
                  "prob"      = mubin,
                  "pobs0"     = 1 - onempstr0 + (onempstr0)*(1-mubin)^nvec,  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempstr0 = lonempstr0, .lprob = lprob,
           .eonempstr0 = eonempstr0, .eprob = eprob,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    misc$link <-    c("prob" = .lprob , "onempstr0" = .lonempstr0 )

    misc$earg <- list("prob" = .eprob , "onempstr0" = .eonempstr0 )

    misc$imethod <- .imethod


      misc$pobs0 <- phi + (1 - phi) * (1 - mubin)^w  # [1]  # P(Y=0)
      misc$pstr0 <- phi
  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob,
            .imethod = imethod ))),
  linkfun = eval(substitute(function(mu, extra = NULL) {
    cbind(theta2eta(mu[, 1], .lprob      , earg = .eprob      ),
          theta2eta(mu[, 2], .lonempstr0 , earg = .eonempstr0 ))
  }, list( .lonempstr0 = lonempstr0, .lprob = lprob,
           .eonempstr0 = eonempstr0, .eprob = eprob ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    mubin     <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempstr0 <- eta2theta(eta[, 2], .lonempstr0 , earg = .eonempstr0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        dzibinom(x = round(w * y), size = w, prob = mubin,
                 log = TRUE, pstr0 = 1 - onempstr0)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempstr0 = lonempstr0, .lprob = lprob,
           .eonempstr0 = eonempstr0, .eprob = eprob ))),
  vfamily = c("zibinomialff"),
  deriv = eval(substitute(expression({
    mubin     <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempstr0 <- eta2theta(eta[, 2], .lonempstr0 , earg = .eonempstr0 )
    omphi     <-     onempstr0
    phi       <- 1 - onempstr0


    prob0 <- (1 - mubin)^w  # Actually q^w
    pobs0 <- phi + (omphi) * prob0
    index <- (y == 0)
    dl.domphi <- -(1 - prob0) / pobs0  # Note "-"
    dl.domphi[!index] <- +1 / (omphi[!index])  # Note "+"

    dl.dmubin <- -w * (omphi) * (1 - mubin)^(w - 1) / pobs0
    dl.dmubin[!index] <- w[!index] *
        (    y[!index]  /      mubin[!index]   -
        (1 - y[!index]) / (1 - mubin[!index]))

    dmubin.deta <- dtheta.deta(mubin, .lprob      , earg = .eprob      )
    domphi.deta <- dtheta.deta(omphi, .lonempstr0 , earg = .eonempstr0 )

    ans <- cbind(dl.dmubin * dmubin.deta,
                 dl.domphi * domphi.deta)

      if ( .lprob == "logit") {
        ans[!index, 1] <- w[!index] * (y[!index] - mubin[!index])
      }

      ans
  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob ))),
  weight = eval(substitute(expression({
    wz <- matrix(NA_real_, nrow = n, ncol = dimm(M))



    ned2l.domphi2 <- (1 - prob0) / ((omphi) * pobs0)


    ned2l.domphimubin <- +w * ((1 - mubin)^(w - 1)) / pobs0  # Note "+"






    ned2l.dmubin2 <- (w * (omphi) / (mubin * (1 - mubin)^2)) *
                     (1 - mubin - w * mubin *
                     (1 - mubin)^w * phi / pobs0)





    wz[,iam(1, 1, M)] <- ned2l.dmubin2     * dmubin.deta^2
    wz[,iam(2, 2, M)] <- ned2l.domphi2     * domphi.deta^2
    wz[,iam(1, 2, M)] <- ned2l.domphimubin * domphi.deta * dmubin.deta
    if (TRUE) {
      ind6 <- (wz[, iam(1, 1, M)] < .Machine$double.eps)
      if (any(ind6))
        wz[ind6, iam(1, 1, M)] <- .Machine$double.eps
    }
    wz
  }), list( .lonempstr0 = lonempstr0, .lprob = lprob,
            .eonempstr0 = eonempstr0, .eprob = eprob ))))
}










dzibinom <- function(x, size, prob, pstr0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(size), length(prob), length(pstr0))
  if (length(x)     != LLL) x     <- rep(x,     len = LLL);
  if (length(size)  != LLL) size  <- rep(size,  len = LLL);
  if (length(prob)  != LLL) prob  <- rep(prob,  len = LLL);
  if (length(pstr0) != LLL) pstr0 <- rep(pstr0, len = LLL);

  ans <- dbinom(x = x, size = size, prob = prob, log = TRUE)


  ans <- if (log.arg) {
    ifelse(x == 0, log(pstr0 + (1-pstr0) * exp(ans)), log1p(-pstr0) + ans)
  } else {
    ifelse(x == 0,     pstr0 + (1-pstr0) * exp(ans) ,
                    (1-pstr0) * exp(ans))
  }


  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


pzibinom <- function(q, size, prob, pstr0 = 0,
                    lower.tail = TRUE, log.p = FALSE) {

  LLL <- max(length(pstr0), length(size), length(prob), length(q))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);

  ans <- pbinom(q, size, prob, lower.tail = lower.tail, log.p = log.p)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)


  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}


qzibinom <- function(p, size, prob, pstr0 = 0,
                    lower.tail = TRUE, log.p = FALSE) {
  LLL <- max(length(p), length(size), length(prob), length(pstr0))
  p     <- rep(p,     length = LLL)
  size  <- rep(size,  length = LLL)
  prob  <- rep(prob,  length = LLL)
  pstr0 <- rep(pstr0, length = LLL)


  ans <- p 
  ans[p <= pstr0] <- 0 
  ans[p >  pstr0] <-
    qbinom((p[p > pstr0] - pstr0[p > pstr0]) / (1 - pstr0[p > pstr0]),
           size[p > pstr0],
           prob[p > pstr0],
           lower.tail = lower.tail, log.p = log.p)



  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- qposbinom((p[pindex] - Pobs0) / (1 - Pobs0),
                             size = size[pindex],
                             prob = prob[pindex])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


rzibinom <- function(n, size, prob, pstr0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n

  pstr0 <- rep(pstr0, len = use.n)
  size  <- rep(size,  len = use.n)
  prob  <- rep(prob,  len = use.n)

  ans <- rbinom(use.n, size, prob)
  ans[runif(use.n) < pstr0] <- 0



  prob0 <- (1 - prob)^size
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- rposbinom(sum(ind0), size = size[ind0], prob = prob[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}












dzinegbin <- function(x, size, prob = NULL, munb = NULL, pstr0 = 0,
                     log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  LLL <- max(length(pstr0), length(size), length(prob), length(x))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);


  ans <- dnbinom(x = x, size = size, prob = prob, log = log.arg)

  ans <- if (log.arg)
    ifelse(x == 0, log(pstr0+(1-pstr0)*exp(ans)), log1p(-pstr0) + ans) else
    ifelse(x == 0,     pstr0+(1-pstr0)*    ans,       (1-pstr0) * ans)



  prob0 <- prob^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}


pzinegbin <- function(q, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  LLL <- max(length(pstr0), length(size), length(prob), length(q))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);



  ans <- pnbinom(q = q, size = size, prob = prob)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)



  prob0 <- prob^size
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}


qzinegbin <- function(p, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }
  LLL <- max(length(p), length(prob), length(pstr0), length(size))
  if (length(p)     != LLL) p      <- rep(p,     len = LLL)
  if (length(pstr0) != LLL) pstr0  <- rep(pstr0, len = LLL);
  if (length(prob)  != LLL) prob   <- rep(prob,  len = LLL)
  if (length(size)  != LLL) size   <- rep(size,  len = LLL);

  ans <- p 
  ind4 <- (p > pstr0)
  ans[!ind4] <- 0
  ans[ ind4] <- qnbinom(p = (p[ind4] - pstr0[ind4]) / (1 - pstr0[ind4]),
                        size = size[ind4], prob = prob[ind4])





  prob0 <- prob^size
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- qposnegbin((p[pindex] - Pobs0) / (1 - Pobs0),
                              size = size[pindex],
                              prob = prob[pindex])
  }


  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN



  ans
}



rzinegbin <- function(n, size, prob = NULL, munb = NULL, pstr0 = 0) {
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n


  pstr0 <- rep(pstr0, len = use.n)
  size  <- rep(size,  len = use.n)
  prob  <- rep(prob,  len = use.n)


  ans <- rnbinom(n = use.n, size = size, prob = prob)
  ans <- ifelse(runif(use.n) < pstr0, rep(0, use.n), ans)



  prob0 <- rep(prob^size, len = use.n)
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0, na.rm = TRUE)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- rposnegbin(sum(ind0, na.rm = TRUE), size = size[ind0],
                    prob = prob[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}








zinegbinomial.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 zinegbinomial <-
  function(
           zero = "size",
           type.fitted = c("mean", "munb", "pobs0", "pstr0", "onempstr0"),
           nsimEIM = 500,
           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.support = 4000,  # 20160127; I have changed this
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           lpstr0 = "logit", lmunb = "loge", lsize = "loge",
           imethod = 1,
           ipstr0 = NULL,
           imunb =  NULL,
           probs.y = 0.35,
           ishrinkage = 0.95,
           isize = NULL,
           gsize.mux = exp((-12:6)/2)) {




  lpstr0 <- as.list(substitute(lpstr0))
  epstr0 <- link2list(lpstr0)
  lpstr0 <- attr(epstr0, "function.name")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")


  type.fitted <- match.arg(type.fitted,
                   c("mean", "munb", "pobs0", "pstr0", "onempstr0"))[1]



  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")

  ipstr0.small <- 1/64  # A number easily represented exactly
  if (length(ipstr0) &&
     (!is.Numeric(ipstr0, positive = TRUE) ||
      any(ipstr0 >= 1)))
    stop("argument 'ipstr0' must contain values in (0,1)")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("argument 'isize' must contain positive values only")

  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be greater than 50, say")


  new("vglmff",
  blurb = c("Zero-inflated negative binomial\n\n",
            "Links:    ",
            namesof("pstr0", lpstr0, earg = epstr0, tag = FALSE), ", ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), "\n",
            "Mean:     (1 - pstr0) * munb"),
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
         parameters.names = c("pstr0", "munb", "size"),
         eps.trig = .eps.trig ,
         type.fitted  = .type.fitted ,
         nsimEIM = .nsimEIM ,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM, .eps.trig = eps.trig,
           .type.fitted = type.fitted
         ))),

      
  initialize = eval(substitute(expression({
    M1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    
    mynames1 <- param.names("pstr0", NOS)
    mynames2 <- param.names("munb",  NOS)
    mynames3 <- param.names("size",  NOS)
    predictors.names <-
      c(namesof(mynames1, .lpstr0 , earg = .epstr0 , tag = FALSE),
        namesof(mynames2, .lmunb  , earg = .emunb  , tag = FALSE),
        namesof(mynames3, .lsize  , earg = .esize  , tag = FALSE))[
        interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {


      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                           imu = .imunb , ishrinkage = .ishrinkage ,
                           pos.only = TRUE,
                           probs.y = .probs.y )



      if ( is.Numeric( .isize )) {
        size.init <- matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
          munb <- extraargs
          sum(c(w) * dposnegbin(y, munb = munb, size = kmat, log = TRUE))
        }
        
        size.init <- matrix(0, nrow = n, ncol = NOS) 
        for (jay in 1:NOS) {
          size.grid <- .gsize.mux * mean(munb.init[, jay])
          TFvec <- (y[, jay] > 0)
          size.init[, jay] <-
            grid.search(size.grid, objfun = posnegbinomial.Loglikfun,
                        y = y[TFvec, jay],  # x = x[TFvec, ],
                        w = w[TFvec, jay],
                        extraargs = munb.init[TFvec, jay])
        }
      }


      
        if (length( .ipstr0 )) {
          pstr0.init <- matrix( .ipstr0 , n, ncoly, byrow = TRUE)
        } else {
          pstr0.init <- matrix(0, n, ncoly)
          ipstr0.small <- .ipstr0.small  # A number easily represented exactly
          for (jay in 1:NOS) {
            Phi.init <- pmax(ipstr0.small,
                             weighted.mean(y[, jay] == 0, w[, jay]) -
                             dnbinom(0, mu = munb.init[, jay],
                                      size = size.init[, jay]))
            if (mean(Phi.init == ipstr0.small) > 0.95)
              warning("from the initial values only, the data appears to ",
                      "have little or no 0-inflation")
            pstr0.init[, jay] <- Phi.init
          }  # for (jay)
        }
          





        etastart <-
          cbind(theta2eta(pstr0.init, .lpstr0 , earg = .epstr0 ),
                theta2eta(munb.init,  .lmunb  , earg = .emunb  ),
                theta2eta(size.init,  .lsize  , earg = .esize  ))
        etastart <-
          etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize,
            .ipstr0 = ipstr0, .imunb = imunb, .isize = isize,
                                              .gsize.mux = gsize.mux,
            .type.fitted = type.fitted,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .ipstr0.small = ipstr0.small,
            .imethod = imethod ))),
      
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "munb", "pobs0", "pstr0", "onempstr0"))[1]

    pstr0 <- eta2theta(eta[, c(TRUE, FALSE, FALSE)],
                       .lpstr0 , earg = .epstr0 )
    if (type.fitted %in% c("mean", "munb", "pobs0"))
      munb  <- eta2theta(eta[, c(FALSE, TRUE, FALSE)],
                         .lmunb  , earg = .emunb  )

    if (type.fitted %in% c("pobs0")) {
      kmat  <- eta2theta(eta[, c(FALSE, FALSE, TRUE)],
                         .lsize  , earg = .esize  )

      tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
      prob0  <- tempk^kmat  # p(0) from negative binomial

      smallval <- 1e-3  # Something like this is needed
      if (any(big.size <- munb / kmat < smallval)) {
        prob0[big.size]  <- exp(-munb[big.size])  # The limit as kmat --> Inf
      }
    }

    ans <- switch(type.fitted,
                  "mean"      = (1 - pstr0) * munb,
                  "munb"      = munb,
                  "pobs0"     = pstr0 + (1 - pstr0) * prob0,  # P(Y=0)
                  "pstr0"     =     pstr0,
                  "onempstr0" = 1 - pstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))        
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpstr0 = lpstr0, .lsize = lsize, .lmunb = lmunb,
           .epstr0 = epstr0, .esize = esize, .emunb = emunb,
           .type.fitted = type.fitted ))),
      
  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lpstr0 , length = NOS),
        rep( .lmunb  , length = NOS),
        rep( .lsize  , length = NOS))[interleave.VGAM(M1*NOS, M1 = M1)]
    temp.names <-
      c(mynames1,
        mynames2,
        mynames3)[interleave.VGAM(M1*NOS, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .epstr0
      misc$earg[[M1*ii-1]] <- .emunb
      misc$earg[[M1*ii  ]] <- .esize
    }

    misc$ipstr0  <- .ipstr0
    misc$isize <- .isize

    misc$max.chunk.MB <- .max.chunk.MB
    misc$cutoff.prob <- .cutoff.prob
    misc$imethod <- .imethod 
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$ishrinkage <- .ishrinkage
    misc$multipleResponses <- TRUE
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize,
            .ipstr0 = ipstr0,                 .isize = isize,
            .nsimEIM = nsimEIM, .imethod = imethod,
            .cutoff.prob = cutoff.prob,
            .max.chunk.MB = max.chunk.MB,
            .ishrinkage = ishrinkage
           ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
   pstr0 <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lpstr0 , earg = .epstr0 )
   munb  <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lmunb  , earg = .emunb  )
   kmat  <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lsize  , earg = .esize  )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzinegbin(x = y, size = kmat, munb = munb,
                         pstr0 = pstr0, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
           .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zinegbinomial"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
   pstr0 <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lpstr0 , earg = .epstr0 )
   munb  <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lmunb  , earg = .emunb  )
   kmat  <- eta2theta(eta[, c(FALSE, FALSE, TRUE)], .lsize  , earg = .esize  )
    rzinegbin(nsim * length(munb),
              size = kmat, munb = munb, pstr0 = pstr0)
  }, list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
           .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),





  validparams = eval(substitute(function(eta, extra = NULL) {
    M1 <- 3
    NOS <- ncol(eta) / M1

    pstr0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                       .lpstr0 , earg = .epstr0 )
    munb  <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                       .lmunb  , earg = .emunb  )
    size  <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                       .lsize  , earg = .esize  )

    smallval <- 1e-3
    ans <- all(is.finite(munb))  && all(munb  > 0) &&
           all(is.finite(size))  && all(size  > 0) &&
           all(is.finite(pstr0)) && all(pstr0 > 0) &&
                                    all(pstr0 < 1) &&
           (overdispersion <- all(munb / size > smallval))
    if (!overdispersion)
        warning("parameter 'size' has very large values; ",
                "replacing them by an arbitrary large value within ",
                "the parameter space. Try fitting ",
                "a zero-inflated Poisson ",
                "model instead.")
    ans
  }, list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
           .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),



    deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- ncol(eta) / M1

    pstr0 <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                       .lpstr0 , earg = .epstr0 )
    munb  <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                       .lmunb  , earg = .emunb  )
    kmat  <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                       .lsize  , earg = .esize  )

    dpstr0.deta <- dtheta.deta(pstr0, .lpstr0 , earg = .epstr0 )
    dmunb.deta  <- dtheta.deta(munb , .lmunb  , earg = .emunb  )
    dsize.deta  <- dtheta.deta(kmat , .lsize  , earg = .esize  )
    dthetas.detas <-
        (cbind(dpstr0.deta,
               dmunb.deta,
               dsize.deta))[, interleave.VGAM(M1*NOS, M1 = M1)]



    smallval <- 1e-2  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
        warning("parameter 'size' has very large values; ",
                "try fitting a zero-inflated Poisson ",
                "model instead")
        kmat[big.size] <- munb[big.size] / smallval
    }



    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    AA16 <- tempm + log(tempk)
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- prob0 * AA16
    df02.dmunb2 <- prob0 * tempk * (1 + 1/kmat) / (1 + munb/kmat)
    df02.dkmat2 <- prob0 * ((tempm^2) / kmat + AA16^2)
    df02.dkmat.dmunb <- -prob0 * (tempm/kmat + AA16) / (1 + munb/kmat)



    AA <- pobs0 <- cbind(pstr0 + (1 - pstr0) * prob0)







    dl.dpstr0 <- -1 / (1 - pstr0)
    dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat)
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y - munb) / (munb + kmat) + log(tempk)



    if (any(big.size)) {
      dl.dsize[big.size] <- 1e-7  # A small number
    }


    for (spp. in 1:NOS) {
      index0 <- (y[, spp.] == 0)
      if (all(index0) || all(!index0))
        stop("must have some 0s AND some positive counts in the data")

      pstr0. <- pstr0[index0, spp.]


      tempk. <- tempk[index0, spp.]  # kmat. / (kmat. + munb.)
      tempm. <- tempm[index0, spp.]  # munb. / (kmat. + munb.)
      prob0. <- prob0[index0, spp.]  # tempk.^kmat.
      df0.dmunb.  <- df0.dmunb[index0, spp.]  # -tempk.* prob0.
      df0.dkmat.  <- df0.dkmat[index0, spp.]  # prob0. * (tempm. + log(tempk.))

      denom. <- AA[index0, spp.]  # pstr0. + (1 - pstr0.) * prob0.
     dl.dpstr0[index0, spp.]  <- (1 - prob0.) / denom.
      dl.dmunb[index0, spp.]  <- (1 - pstr0.) * df0.dmunb. / denom.
      dl.dsize[index0, spp.]  <- (1 - pstr0.) * df0.dkmat. / denom.
    }  # of spp.


    dl.dthetas <-
      cbind(dl.dpstr0,
            dl.dmunb,
            dl.dsize)[, interleave.VGAM(M1*NOS, M1 = M1)]


    ans <- c(w) * dl.dthetas * dthetas.detas
    ans
  }), list( .lpstr0 = lpstr0, .lmunb = lmunb, .lsize = lsize,
            .epstr0 = epstr0, .emunb = emunb, .esize = esize ))),



  weight = eval(substitute(expression({


    wz <- matrix(0, n, M + M-1 + M-2)
    mymu <- munb / oneminusf0  # Is the same as 'mu', == E(Y)

    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB




    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 1
      Q.maxs <-      qposnegbin(p    = eff.p[2] ,
                                munb = munb[, jay],
                                size = kmat[, jay]) + 10




      eps.trig <- .eps.trig
      Q.MAXS <-      pmax(10, ceiling(1 / sqrt(eps.trig)))
      Q.maxs <- pmin(Q.maxs, Q.MAXS)



      ind1 <- if (max.chunk.MB > 0) (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]


          wz[sind2, M1*jay] <-
            EIM.posNB.specialp(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.max = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only,
                               second.deriv = FALSE)
  if (FALSE)
          wz2[sind2, M1*jay] <-
            EIM.posNB.speciald(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.min = min(Q.mins2[sind2]),
                               y.max = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only,
                               second.deriv = FALSE)




          wz[sind2, M1*jay] <-
          wz[sind2, M1*jay] * (1 - AA[sind2, jay]) -
          (1-pstr0[sind2, jay]) * (df02.dkmat2[sind2, jay] -
          (1-pstr0[sind2, jay]) * (df0.dkmat[sind2, jay]^2) / AA[sind2, jay])



          if (any(eim.kk.TF <-       wz[sind2, M1*jay] <= 0 |
                               is.na(wz[sind2, M1*jay]))) {
            ind2[sind2[eim.kk.TF], jay] <- FALSE
          }



          lwr.ptr <- upr.ptr + 1
        }  # while

      }
    }  # end of for (jay in 1:NOS)







    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        kkvec <-  kmat[ii.TF, jay]
        muvec <-  munb[ii.TF, jay]
        PSTR0 <- pstr0[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rzinegbin(sum(ii.TF), pstr0 = PSTR0,
                            mu = muvec, size = kkvec)

          index0 <- (ysim == 0)


          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) -
                   (ysim - muvec) / (muvec + kkvec) +
                   log1p(-muvec / (kkvec + muvec))  # +

          ans0 <- (1 - PSTR0) *
            df0.dkmat[ii.TF , jay] / AA[ii.TF , jay]
          dl.dk[index0] <- ans0[index0]

          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay] <- ned2l.dk2  # * (dsize.deta[ii.TF, jay])^2
      }
    }



    wz[, M1*(1:NOS)    ] <- wz[, M1*(1:NOS)    ] * dsize.deta^2




    save.weights <- !all(ind2)


    ned2l.dpstr02 <- oneminusf0 / (AA * (1 - pstr0))
    wz[,     M1*(1:NOS) - 2] <- ned2l.dpstr02 * dpstr0.deta^2


    ned2l.dpstr0.dmunb <- df0.dmunb / AA
    wz[, M + M1*(1:NOS) - 2] <- ned2l.dpstr0.dmunb *
                                dpstr0.deta * dmunb.deta

    ned2l.dpstr0.dsize <- df0.dkmat / AA
    wz[, M + M-1 + M1*(1:NOS) - 2] <- ned2l.dpstr0.dsize *
                                      dpstr0.deta * dsize.deta



    ned2l.dmunb2 <-
      (1 - AA) * (mymu / munb^2 -
                   ((1 + mymu/kmat) / kmat) / (1 + munb/kmat)^2) -
      (1-pstr0) * (df02.dmunb2 -
                  (1 - pstr0) * (df0.dmunb^2) / AA)

        wz[,     M1*(1:NOS) - 1] <- ned2l.dmunb2 * dmunb.deta^2


    dAA.dmunb <- (1 - pstr0) * df0.dmunb



    ned2l.dmunbsize <-
      (1 - AA) * (munb - mymu) / (munb + kmat)^2 -
      (1-pstr0) * (df02.dkmat.dmunb -
                   df0.dkmat * dAA.dmunb / AA)

    wz[, M +       M1*(1:NOS) - 1] <- ned2l.dmunbsize * dmunb.deta *
                                                        dsize.deta

    



    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .lpstr0 = lpstr0,
            .epstr0 = epstr0, .nsimEIM = nsimEIM,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB  = max.chunk.MB ))))
}  # End of zinegbinomial








zinegbinomialff.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}


 zinegbinomialff <-
  function(lmunb = "loge", lsize = "loge", lonempstr0 = "logit", 
           type.fitted = c("mean", "munb", "pobs0", "pstr0", "onempstr0"),
           imunb = NULL, isize = NULL, ionempstr0 = NULL,  
           zero = c("size", "onempstr0"),
           imethod = 1, ishrinkage = 0.95,

           probs.y = 0.35,
           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.support = 4000,  # 20160127; I have changed this
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           gsize.mux = exp((-12:6)/2),

           nsimEIM = 500) {


  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")

  ipstr0.small <- 1/64  # A number easily represented exactly

  type.fitted <- match.arg(type.fitted,
                   c("mean", "munb", "pobs0", "pstr0", "onempstr0"))[1]



  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")

  if (length(ionempstr0) &&
     (!is.Numeric(ionempstr0, positive = TRUE) ||
      any(ionempstr0 >= 1)))
    stop("argument 'ionempstr0' must contain values in (0,1)")
  if (length(isize) && !is.Numeric(isize, positive = TRUE))
    stop("argument 'isize' must contain positive values only")

  if (!is.Numeric(nsimEIM, length.arg = 1, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 50)
    warning("argument 'nsimEIM' should be greater than 50, say")



  new("vglmff",
  blurb = c("Zero-inflated negative binomial\n\n",
            "Links:    ",
            namesof("munb",  lmunb,  earg = emunb,  tag = FALSE), ", ",
            namesof("size",  lsize,  earg = esize,  tag = FALSE), ", ",
            namesof("onempstr0", lonempstr0, earg = eonempstr0, tag = FALSE),
            "\n",
            "Mean:     (1 - pstr0) * munb"),
  constraints = eval(substitute(expression({

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 3)
  }), list( .zero = zero ))),


  infos = eval(substitute(function(...) {
    list(M1 = 3,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("munb", "size", "onempstr0"),
         eps.trig = .eps.trig ,
         nsimEIM = .nsimEIM ,
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .nsimEIM = nsimEIM, .eps.trig = eps.trig,
           .type.fitted = type.fitted
         ))),

      
  initialize = eval(substitute(expression({
    M1 <- 3

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y


    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted <- .type.fitted
    extra$dimnamesy   <- dimnames(y)


    
    mynames1 <- param.names("munb",       NOS)
    mynames2 <- param.names("size",       NOS)
    mynames3 <- param.names("onempstr0",  NOS) 
    predictors.names <-
      c(namesof(mynames1, .lmunb  , earg = .emunb  , tag = FALSE),
        namesof(mynames2, .lsize  , earg = .esize  , tag = FALSE),
        namesof(mynames3, .lonempstr0 , earg = .eonempstr0 , tag = FALSE))[
        interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {

      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                           imu = .imunb , ishrinkage = .ishrinkage ,
                           pos.only = TRUE,
                           probs.y = .probs.y )



      if ( is.Numeric( .isize )) {
        size.init <- matrix( .isize , nrow = n, ncol = ncoly, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kmat, y, x, w, extraargs) {
          munb <- extraargs
          sum(c(w) * dposnegbin(y, munb = munb, size = kmat, log = TRUE))
        }
        
        size.init <- matrix(0, nrow = n, ncol = NOS) 
        for (jay in 1:NOS) {
          size.grid <- .gsize.mux * mean(munb.init[, jay])
          TFvec <- (y[, jay] > 0)
          size.init[, jay] <-
            grid.search(size.grid, objfun = posnegbinomial.Loglikfun,
                        y = y[TFvec, jay],  # x = x[TFvec, ],
                        w = w[TFvec, jay],
                        extraargs = munb.init[TFvec, jay])
        }
      }


      
        if (length( .ionempstr0 )) {
          onempstr0.init <- matrix( .ionempstr0 , n, ncoly, byrow = TRUE)
        } else {
          onempstr0.init <- matrix(0, n, ncoly)
          ipstr0.small <- .ipstr0.small  # Easily represented exactly
          for (jay in 1:NOS) {
            Phi.init <- pmax(ipstr0.small,
                             weighted.mean(y[, jay] == 0, w[, jay]) -
                             dnbinom(0, mu = munb.init[, jay],
                                      size = size.init[, jay]))
            if (mean(Phi.init == ipstr0.small) > 0.95)
              warning("from the initial values only, the data appears to ",
                      "have little or no 0-inflation")
            onempstr0.init[, jay] <- 1 - Phi.init
          }  # for (jay)
        }
          




        etastart <-
          cbind(theta2eta(munb.init,   .lmunb  , earg = .emunb  ),
                theta2eta(size.init,   .lsize  , earg = .esize  ),
                theta2eta(onempstr0.init, .lonempstr0 ,
                          earg = .eonempstr0 ))
        etastart <-
          etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
            .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize,
            .ionempstr0 = ionempstr0, .imunb = imunb, .isize = isize,
                                                      .gsize.mux = gsize.mux,
            .type.fitted = type.fitted,
            .ipstr0.small = ipstr0.small,
            .ishrinkage = ishrinkage, .probs.y = probs.y,
            .imethod = imethod ))),
      
  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "munb", "pobs0", "pstr0", "onempstr0"))[1]

    M1 <- 3
    NOS <- ncol(eta) / M1
    if (type.fitted %in% c("mean", "munb", "pobs0"))
      munb    <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb  , earg = .emunb  )

    if (type.fitted %in% c("pobs0")) {
      kmat    <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize , earg = .esize )

      tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
      prob0  <- tempk^kmat  # p(0) from negative binomial

      smallval <- 1e-3  # Something like this is needed
      if (any(big.size <- munb / kmat < smallval)) {
        prob0[big.size]  <- exp(-munb[big.size])  # The limit as kmat --> Inf
      }
    }

    onempstr0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempstr0 , earg = .eonempstr0 )


    ans <- switch(type.fitted,
                  "mean"      = onempstr0 * munb,
                  "munb"      = munb,
                  "pobs0"     = 1 - onempstr0 + onempstr0 * prob0,  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))        
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempstr0 = lonempstr0, .lsize = lsize, .lmunb = lmunb,
           .eonempstr0 = eonempstr0, .esize = esize, .emunb = emunb,
           .type.fitted = type.fitted ))),
      
  last = eval(substitute(expression({
    misc$link <-
      c(rep( .lmunb      , length = NOS),
        rep( .lsize      , length = NOS),
        rep( .lonempstr0 , length = NOS))[interleave.VGAM(M1*NOS, M1 = M1)]
    temp.names <-
      c(mynames1,
        mynames2,
        mynames3)[interleave.VGAM(M1*NOS, M1 = M1)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- temp.names
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-2]] <- .emunb
      misc$earg[[M1*ii-1]] <- .esize
      misc$earg[[M1*ii  ]] <- .eonempstr0
    }

    misc$imethod <- .imethod
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$M1 <- M1
    misc$ionempstr0  <- .ionempstr0
    misc$isize <- .isize
    misc$multipleResponses <- TRUE

  }), list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
            .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize,
            .ionempstr0 = ionempstr0,                 .isize = isize,
            .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    M1 <- 3
    NOS <- extra$NOS
    munb      <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb , earg = .emunb )
    kmat      <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize , earg = .esize )
    onempstr0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempstr0 , earg = .eonempstr0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzinegbin(x = y, size = kmat, munb = munb,
                         pstr0 = 1 - onempstr0, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
           .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize ))),
  vfamily = c("zinegbinomialff"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    munb <- eta2theta(eta[, c(TRUE, FALSE, FALSE)], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, c(FALSE, TRUE, FALSE)], .lsize , earg = .esize )
    onempstr0 <- eta2theta(eta[, c(FALSE, FALSE, TRUE)],
                       .lpstr0 , earg = .epstr0 )
    rzinegbin(nsim * length(munb),
              size = kmat, munb = munb, pstr0 = 1 - onempstr0)
  }, list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
           .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize ))),




  deriv = eval(substitute(expression({
    M1 <- 3
    NOS <- ncol(eta) / M1

    munb      <- eta2theta(eta[, M1*(1:NOS)-2, drop = FALSE],
                           .lmunb  , earg = .emunb  )
    kmat      <- eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                           .lsize  , earg = .esize  )
    onempstr0 <- eta2theta(eta[, M1*(1:NOS)  , drop = FALSE],
                           .lonempstr0 , earg = .eonempstr0 )

    donempstr0.deta <- dtheta.deta(onempstr0, .lonempstr0 ,
                                   earg = .eonempstr0 )
    dmunb.deta  <- dtheta.deta(munb , .lmunb , earg = .emunb )
    dsize.deta  <- dtheta.deta(kmat , .lsize , earg = .esize )
    dthetas.detas <-
        (cbind(dmunb.deta,
               dsize.deta,
               donempstr0.deta))[, interleave.VGAM(M1*NOS, M1 = M1)]




    smallval <- 1e-2  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
        warning("parameter 'size' has very large values; ",
                "try fitting a zero-inflated Poisson ",
                "model instead")
        kmat[big.size] <- munb[big.size] / smallval
    }



    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
    tempm <- munb / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0
    AA16 <- tempm + log(tempk)
    df0.dmunb   <- -tempk * prob0
    df0.dkmat   <- cbind(prob0 * AA16)
    df02.dmunb2 <- prob0 * tempk * (1 + 1/kmat) / (1 + munb/kmat)
    df02.dkmat2 <- prob0 * ((tempm^2) / kmat + AA16^2)
    df02.dkmat.dmunb <- -prob0 * (tempm/kmat + AA16) / (1 + munb/kmat)



    pstr0 <- 1 - onempstr0
    AA <- pobs0 <- cbind(pstr0 + (onempstr0) * prob0)


    dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat)
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y - munb) / (munb + kmat) + log(tempk)
    dl.donempstr0 <- +1 / (onempstr0)



    for (spp. in 1:NOS) {
      index0 <- (y[, spp.] == 0)
      if (all(index0) || all(!index0))
        stop("must have some 0s AND some positive counts in the data")

      kmat.      <-      kmat[index0, spp.]
      munb.      <-      munb[index0, spp.]
      onempstr0. <- onempstr0[index0, spp.]


      tempk. <- kmat. / (kmat. + munb.)
      tempm. <- munb. / (kmat. + munb.)
      prob0. <- tempk.^kmat.
      df0.dmunb.  <- -tempk.* prob0.
      df0.dkmat.  <- prob0. * (tempm. + log(tempk.))

      denom. <- 1 - onempstr0. + (onempstr0.) * prob0.
     dl.donempstr0[index0, spp.]  <- -(1 - prob0.) / denom.  # note "-"
          dl.dmunb[index0, spp.]  <- (onempstr0.) * df0.dmunb. / denom.
          dl.dsize[index0, spp.]  <- (onempstr0.) * df0.dkmat. / denom.
    }  # of spp.


    dl.dthetas <-
      cbind(dl.dmunb,
            dl.dsize,
            dl.donempstr0)[, interleave.VGAM(M1*NOS, M1 = M1)]


      c(w) * dl.dthetas * dthetas.detas
  }), list( .lonempstr0 = lonempstr0, .lmunb = lmunb, .lsize = lsize,
            .eonempstr0 = eonempstr0, .emunb = emunb, .esize = esize ))),




  weight = eval(substitute(expression({



    wz <- matrix(0, n, M + M-1 + M-2)
    mymu <- munb / oneminusf0  # Is the same as 'mu', == E(Y)

    max.support <- .max.support
    max.chunk.MB <- .max.chunk.MB




    ind2 <- matrix(FALSE, n, NOS)  # Used for SFS
    for (jay in 1:NOS) {
      eff.p <- sort(c( .cutoff.prob , 1 - .cutoff.prob ))
      Q.mins <- 1
      Q.maxs <-      qposnegbin(p    = eff.p[2] ,
                                munb = munb[, jay],
                                size = kmat[, jay]) + 10




      eps.trig <- .eps.trig
      Q.MAXS <- pmax(10, ceiling(1 / sqrt(eps.trig)))
      Q.maxs <- pmin(Q.maxs, Q.MAXS)




      ind1 <- if (max.chunk.MB > 0) (Q.maxs - Q.mins < max.support) else FALSE
      if ((NN <- sum(ind1)) > 0) {
        Object.Size <- NN * 8 * max(Q.maxs - Q.mins) / (2^20)
        n.chunks <- if (intercept.only) 1 else
                    max(1, ceiling( Object.Size / max.chunk.MB))
        chunk.rows <- ceiling(NN / n.chunks)
        ind2[, jay] <- ind1  # Save this
        wind2 <- which(ind1)


        upr.ptr <- 0
        lwr.ptr <- upr.ptr + 1
        while (lwr.ptr <= NN) {
          upr.ptr <- min(upr.ptr + chunk.rows, NN)
          sind2 <- wind2[lwr.ptr:upr.ptr]
          wz[sind2, M1*jay - 1] <-
            EIM.posNB.specialp(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.max = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only,
                               second.deriv = FALSE)
  if (FALSE)
          wz2[sind2, M1*jay - 1] <-
            EIM.posNB.speciald(munb        = munb[sind2, jay],
                               size        = kmat[sind2, jay],
                               y.min = min(Q.mins2[sind2]),
                               y.max = max(Q.maxs[sind2]),
                               cutoff.prob = .cutoff.prob ,
                               prob0       =       prob0[sind2, jay],
                               df0.dkmat   =   df0.dkmat[sind2, jay],
                               df02.dkmat2 = df02.dkmat2[sind2, jay],
                               intercept.only = intercept.only,
                               second.deriv = FALSE)






          wz[sind2, M1*jay - 1] <-
          wz[sind2, M1*jay - 1] * (1 - AA[sind2, jay]) -
          (1-pstr0[sind2, jay]) * (df02.dkmat2[sind2, jay] -
          (1-pstr0[sind2, jay]) * (df0.dkmat[sind2, jay]^2) / AA[sind2, jay])



          if (any(eim.kk.TF <-       wz[sind2, M1*jay - 1] <= 0 |
                               is.na(wz[sind2, M1*jay - 1]))) {
            ind2[sind2[eim.kk.TF], jay] <- FALSE
          }



          lwr.ptr <- upr.ptr + 1
        }  # while

      }
    }  # end of for (jay in 1:NOS)









    for (jay in 1:NOS) {
      run.varcov <- 0
      ii.TF <- !ind2[, jay]  # Not assigned above
      if (any(ii.TF)) {
        kkvec <-  kmat[ii.TF, jay]
        muvec <-  munb[ii.TF, jay]
        PSTR0 <- pstr0[ii.TF, jay]
        for (ii in 1:( .nsimEIM )) {
          ysim <- rzinegbin(sum(ii.TF), pstr0 = PSTR0,
                            mu = muvec, size = kkvec)

          index0 <- (ysim == 0)


          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) -
                   (ysim - muvec) / (muvec + kkvec) +
                   log1p(-muvec / (kkvec + muvec))  # +

          ans0 <- (1 - PSTR0) *
            df0.dkmat[ii.TF , jay] / AA[ii.TF , jay]
          dl.dk[index0] <- ans0[index0]

          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay - 1] <- ned2l.dk2  # * (dsize.deta[ii.TF, jay])^2
      }
    }



    wz[, M1*(1:NOS) - 1] <- wz[, M1*(1:NOS) - 1] * dsize.deta^2




    save.weights <- !all(ind2)


    ned2l.donempstr02 <- oneminusf0 / (AA * (onempstr0))
    wz[,     M1*(1:NOS)    ] <- ned2l.donempstr02 * donempstr0.deta^2

    ned2l.donempstr0.dmunb <- -df0.dmunb / AA  # Negated (1/2)
    wz[, M + M-1 + M1*(1:NOS) - 2] <- ned2l.donempstr0.dmunb *
                                      donempstr0.deta * dmunb.deta

    ned2l.donempstr0.dsize <- -df0.dkmat / AA  # Negated (2/2)
    wz[, M       + M1*(1:NOS) - 1] <- ned2l.donempstr0.dsize *
                                      donempstr0.deta * dsize.deta

    ned2l.dmunb2 <-
      (1 - AA) * (mymu / munb^2 -
                   ((1 + mymu/kmat) / kmat) / (1 + munb/kmat)^2) -
      (1-pstr0) * (df02.dmunb2 -
                  (1 - pstr0) * (df0.dmunb^2) / AA)
    wz[,     M1*(1:NOS) - 2] <- ned2l.dmunb2 * dmunb.deta^2


    dAA.dmunb <- (onempstr0) * df0.dmunb
    ned2l.dmunbsize <-
      (1 - AA) * (munb - mymu) / (munb + kmat)^2 -
      (onempstr0) * (df02.dkmat.dmunb -
                     df0.dkmat * dAA.dmunb / AA)

    wz[, M +       M1*(1:NOS) - 2] <- ned2l.dmunbsize * dmunb.deta *
                                                        dsize.deta


    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .lonempstr0 = lonempstr0,
            .eonempstr0 = eonempstr0, .nsimEIM = nsimEIM,
            .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB  = max.chunk.MB ))))
}  # End of zinegbinomialff










 zipoissonff <-
  function(llambda = "loge", lonempstr0 = "logit",
           type.fitted = c("mean", "lambda", "pobs0", "pstr0", "onempstr0"),
           ilambda = NULL,   ionempstr0 = NULL,
           gonempstr0 = NULL,  # (1:9) / 10, 
           imethod = 1,
           ishrinkage = 0.95, probs.y = 0.35,
           zero = "onempstr0") {


  type.fitted <- match.arg(type.fitted,
                   c("mean", "lambda", "pobs0", "pstr0", "onempstr0"))[1]



  llambda <- as.list(substitute(llambda))
  elambda <- link2list(llambda)
  llambda <- attr(elambda, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")

  ipstr0.small <- 1/64  # A number easily represented exactly


  if (length(ilambda))
    if (!is.Numeric(ilambda, positive = TRUE))
      stop("'ilambda' values must be positive")
  if (length(ionempstr0))
    if (!is.Numeric(ionempstr0, positive = TRUE) ||
      any(ionempstr0 >= 1))
      stop("'ionempstr0' values must be inside the interval (0,1)")


  new("vglmff",
  blurb = c("Zero-inflated Poisson\n\n",
            "Links:    ",
            namesof("lambda",    llambda,    earg = elambda), ", ",
            namesof("onempstr0", lonempstr0, earg = eonempstr0), "\n",
            "Mean:     onempstr0 * lambda"),
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
         parameters.names = c("lambda", "onempstr0"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    ncoly <- ncol(y)
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)

    mynames1 <- param.names("lambda",    ncoly)
    mynames2 <- param.names("onempstr0", ncoly)
    predictors.names <-
      c(namesof(mynames1, .llambda    , earg = .elambda    , tag = FALSE),
        namesof(mynames2, .lonempstr0 , earg = .eonempstr0 , tag = FALSE))[
          interleave.VGAM(M, M1 = M1)]


      if (!length(etastart)) {
      matL <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                      imu = .ilambda , ishrinkage = .ishrinkage ,
                      pos.only = TRUE,
                      probs.y = .probs.y )

        matP <- matrix(if (length( .ionempstr0 )) .ionempstr0 else 0,
                       n, ncoly, byrow = TRUE)
        phi0.grid <- .gonempstr0
        ipstr0.small <- .ipstr0.small  # A number easily represented exactly

        if (!length( .ionempstr0 ))
        for (jay in 1:ncoly) {
          zipois.Loglikfun <- function(phival, y, x, w, extraargs) {
            sum(c(w) * dzipois(x = y, pstr0 = phival,
                            lambda = extraargs$lambda,
                            log = TRUE))
          }
        Phi0.init <- if (length(phi0.grid)) {
          grid.search(phi0.grid,
                      objfun = zipois.Loglikfun,
                      y = y[, jay], x = x, w = w[, jay],
                      extraargs = list(lambda = matL[, jay]))
        } else {
          pmax(ipstr0.small,
               weighted.mean(y[, jay] == 0, w[, jay]) -
               dpois(0, matL[, jay]))
        }
        if (mean(Phi0.init == ipstr0.small) > 0.95)
          warning("from the initial values only, the data appears to ",
                  "have little or no 0-inflation")

          matP[, jay] <- Phi0.init
      }  # for (jay)

      etastart <-
        cbind(theta2eta(    matL, .llambda    , earg = .elambda    ),
              theta2eta(1 - matP, .lonempstr0 , earg = .eonempstr0 ))[,
                        interleave.VGAM(M, M1 = M1)]

      mustart <- NULL  # Since etastart has been computed.
    }
  }), list( .lonempstr0 = lonempstr0, .llambda = llambda,
            .eonempstr0 = eonempstr0, .elambda = elambda,
            .ionempstr0 = ionempstr0, .ilambda = ilambda,
            .gonempstr0 = gonempstr0,
            .type.fitted = type.fitted, .probs.y = probs.y,
            .ipstr0.small = ipstr0.small,
            .imethod = imethod, .ishrinkage = ishrinkage ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "pobs0", "pstr0", "onempstr0"))[1]

    M1 <- 2
    ncoly <- ncol(eta) / M1
    lambda    <- eta2theta(eta[, M1*(1:ncoly) - 1], .llambda ,
                           earg = .elambda )
    onempstr0 <- eta2theta(eta[, M1*(1:ncoly)    ], .lonempstr0 ,
                           earg = .eonempstr0 )


    ans <- switch(type.fitted,
                  "mean"      = onempstr0 * lambda,
                  "lambda"    = lambda,
                  "pobs0"     = 1 + onempstr0 * expm1(-lambda),  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      if (length(extra$dimnamesy[[1]]) == nrow(ans))
        dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempstr0 = lonempstr0, .llambda = llambda,
           .eonempstr0 = eonempstr0, .elambda = elambda,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    M1 <- extra$M1
    misc$link <-
      c(rep( .llambda    , length = ncoly),
        rep( .lonempstr0 , length = ncoly))[interleave.VGAM(M, M1 = M1)]
    temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
    names(misc$link) <- temp.names


    misc$earg <- vector("list", M1 * ncoly)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii-1]] <- .elambda
      misc$earg[[M1*ii  ]] <- .eonempstr0
    }

    misc$M1 <- M1
    misc$imethod <- .imethod
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE

      misc$pobs0 <- (1 - onempstr0) + onempstr0 * exp(-lambda)  # P(Y=0)
      misc$pobs0 <- as.matrix(misc$pobs0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pobs0) <- dimnames(y)

      misc$pstr0 <- (1 - onempstr0)
      misc$pstr0 <- as.matrix(misc$pstr0)
      if (length(dimnames(y)[[2]]) > 0)
        dimnames(misc$pstr0) <- dimnames(y)
  }), list( .lonempstr0 = lonempstr0, .llambda = llambda,
            .eonempstr0 = eonempstr0, .elambda = elambda,
            .imethod = imethod ))),
  loglikelihood = eval(substitute( 
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda    <- eta2theta(eta[, c(TRUE, FALSE)], .llambda    ,
                           earg = .elambda )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) *
                 dzipois(x = y, pstr0 = 1 - onempstr0, lambda = lambda,
                         log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempstr0 = lonempstr0, .llambda = llambda,
           .eonempstr0 = eonempstr0, .elambda = elambda ))),
  vfamily = c("zipoissonff"),



  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda    <- eta2theta(eta[, c(TRUE, FALSE)], .llambda ,
                           earg = .elambda    )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )
    rzipois(nsim * length(lambda), lambda = lambda, pstr0 = 1 - onempstr0)
  }, list( .lonempstr0 = lonempstr0, .llambda = llambda,
           .eonempstr0 = eonempstr0, .elambda = elambda ))),



  deriv = eval(substitute(expression({
    M1 <- 2
    ncoly <- ncol(eta) / M1  # extra$ncoly
    lambda    <- eta2theta(eta[, c(TRUE, FALSE)], .llambda ,
                           earg = .elambda    )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )


    dlambda.deta    <- dtheta.deta(lambda   , .llambda    ,
                                   earg = .elambda )
    donempstr0.deta <- dtheta.deta(onempstr0, .lonempstr0 ,
                                   earg = .eonempstr0 )

    denom <- 1 + onempstr0 * expm1(-lambda)
    ind0 <- (y == 0)
    dl.dlambda <- -onempstr0 * exp(-lambda) / denom
    dl.dlambda[!ind0] <- (y[!ind0] - lambda[!ind0]) / lambda[!ind0]
    dl.donempstr0 <- expm1(-lambda) / denom
    dl.donempstr0[!ind0] <- 1 / onempstr0[!ind0]

    ans <- c(w) * cbind(dl.dlambda    * dlambda.deta,
                        dl.donempstr0 * donempstr0.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]


    if ( .llambda == "loge" && is.empty.list( .elambda ) &&
       any(lambda[!ind0] < .Machine$double.eps)) {
      for (spp. in 1:ncoly) {
        ans[!ind0[, spp.], M1 * spp.] <-
          w[!ind0[, spp.]] *
         (y[!ind0[, spp.], spp.] - lambda[!ind0[, spp.], spp.])
      }
    }



    ans
  }), list( .lonempstr0 = lonempstr0, .llambda = llambda,
            .eonempstr0 = eonempstr0, .elambda = elambda ))),
  weight = eval(substitute(expression({


    ned2l.dlambda2 <-  (    onempstr0) / lambda -
                    onempstr0 * (1 - onempstr0) * exp(-lambda) / denom
    ned2l.donempstr0.2 <- -expm1(-lambda) / ((onempstr0) * denom)
    ned2l.dphilambda <- +exp(-lambda) / denom


    wz <- array(c(c(w) * ned2l.dlambda2 * dlambda.deta^2,
                  c(w) * ned2l.donempstr0.2 * donempstr0.deta^2,
                  c(w) * ned2l.dphilambda * donempstr0.deta * dlambda.deta),
                dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)

    wz
  }), list( .llambda = llambda ))))
}







dzigeom <- function(x, prob, pstr0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(prob), length(pstr0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);


  ans <- dgeom(x = x, prob = prob, log = TRUE)


  ans <- if (log.arg) {
    ifelse(x == 0, log(pstr0 + (1 - pstr0) * exp(ans)),
                   log1p(-pstr0) + ans)
  } else {
    ifelse(x == 0,     pstr0 + (1 - pstr0) * exp(ans) ,
                               (1 - pstr0) * exp(ans))
  }



  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



pzigeom <- function(q, prob, pstr0 = 0) {


  LLL <- max(length(q), length(prob), length(pstr0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pstr0)  != LLL) pstr0  <- rep(pstr0,  len = LLL);

  ans <- pgeom(q, prob)
  ans <- ifelse(q < 0, 0, pstr0 + (1-pstr0) * ans)


  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



qzigeom <- function(p, prob, pstr0 = 0) {
  LLL <- max(length(p), length(prob), length(pstr0))
  ans <- p <- rep(p,     len = LLL)
  prob     <- rep(prob,  len = LLL)
  pstr0    <- rep(pstr0, len = LLL)
  ans[p <= pstr0] <- 0 
  ind1 <- (p > pstr0)
  ans[ind1] <-
    qgeom((p[ind1] - pstr0[ind1]) / (1 - pstr0[ind1]),
          prob = prob[ind1])


  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[p[ind0] <= pobs0] <- 0 
    pindex <- (1:LLL)[ind0 & (p > pobs0)]
    Pobs0 <- pstr0[pindex] + (1 - pstr0[pindex]) * prob0[pindex]
    ans[pindex] <- 1 + qgeom((p[pindex] - Pobs0) / (1 - Pobs0),
                            prob = prob[pindex])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN

  ans
}



rzigeom <- function(n, prob, pstr0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
              stop("bad input for argument 'n'") else n


  pstr0 <- rep(pstr0, len = use.n)
  prob  <- rep(prob,  len = use.n)


  ans <- rgeom(use.n, prob)
  ans[runif(use.n) < pstr0] <- 0


  prob0 <- prob
  deflat.limit <- -prob0 / (1 - prob0)
  ind0 <- (deflat.limit <= pstr0) & (pstr0 <  0)
  if (any(ind0)) {
    pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * prob0[ind0]
    ans[ind0] <- 1 + rgeom(sum(ind0), prob = prob[ind0])
    ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
  }

  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN


  ans
}




 zigeometric <-
  function(
           lpstr0 = "logit",
           lprob  = "logit",
           type.fitted = c("mean", "prob", "pobs0", "pstr0", "onempstr0"),
           ipstr0  = NULL, iprob = NULL,
           imethod = 1,
           bias.red = 0.5,
           zero = NULL) {



  expected <- TRUE



  lpstr0 <- as.list(substitute(lpstr0))
  epstr0 <- link2list(lpstr0)
  lpstr0 <- attr(epstr0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  type.fitted <- match.arg(type.fitted,
                   c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]


  if (length(ipstr0))
    if (!is.Numeric(ipstr0, positive = TRUE) ||
        ipstr0 >= 1)
      stop("argument 'ipstr0' is out of range")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")

  if (!is.Numeric(bias.red, length.arg = 1, positive = TRUE) ||
     bias.red > 1)
    stop("argument 'bias.red' must be between 0 and 1")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-inflated geometric distribution,\n",
            "P[Y = 0] = pstr0 + (1 - pstr0) * prob,\n",
            "P[Y = y] = (1 - pstr0) * prob * (1 - prob)^y, ",
            "y = 1, 2, ...\n\n",
            "Link:     ",
            namesof("pstr0",  lpstr0,  earg = epstr0), ", ",
            namesof("prob",   lprob,   earg = eprob ), "\n",
            "Mean:     (1 - pstr0) * (1 - prob) / prob"),
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
         parameters.names = c("pstr0", "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted  = type.fitted ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    mynames1 <- param.names("pstr0", ncoly)
    mynames2 <- param.names("prob",  ncoly)
    predictors.names <-
            c(namesof(mynames1, .lpstr0,  earg = .epstr0, tag = FALSE),
              namesof(mynames2, .lprob,   earg = .eprob,  tag = FALSE))[
          interleave.VGAM(M1 * NOS, M1 = M1)]


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 3)
                       .bias.red / (1 + y + 1/8) else
                   if ( .imethod == 2)
                       .bias.red / (1 +
                   matrix(colMeans(y) + 1/8,
                          n, ncoly, byrow = TRUE)) else
                       .bias.red / (1 +
                   matrix(colSums(y * w) / colSums(w) + 1/8,
                          n, ncoly, byrow = TRUE))

      prob.init <- if (length( .iprob )) {
        matrix( .iprob , n, ncoly, byrow = TRUE)
      } else {
        prob.init # Already a matrix
      }


      prob0.est <- psze.init <- matrix(0, n, NOS)
      for (jlocal in 1:NOS) {
        prob0.est[, jlocal] <-
          sum(w[y[, jlocal] == 0, jlocal]) / sum(w[, jlocal])
        psze.init[, jlocal] <- if ( .imethod == 3)
                         prob0.est[, jlocal] / 2 else
                     if ( .imethod == 1)
                         pmax(0.05, (prob0.est[, jlocal] -
                                     median(prob.init[, jlocal]))) else
                         prob0.est[, jlocal] / 5
      }
      psze.init <- if (length( .ipstr0 )) {
        matrix( .ipstr0 , n, ncoly, byrow = TRUE)
      } else {
        psze.init # Already a matrix
      }



      etastart <-
        cbind(theta2eta(psze.init, .lpstr0, earg = .epstr0),
              theta2eta(prob.init, .lprob , earg = .eprob ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
            .iprob = iprob, .ipstr0 = ipstr0,
            .type.fitted = type.fitted,
            .bias.red = bias.red,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )

    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = (1 - pstr0) * (1 - prob) / prob,
                  "prob"      = prob,
                  "pobs0"     = pstr0 + (1 - pstr0) * prob,  # P(Y=0)
                  "pstr0"     =     pstr0,
                  "onempstr0" = 1 - pstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lpstr0 , len = NOS),
                    rep( .lprob  , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names


    misc$earg <- vector("list", M1 * NOS)
    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epstr0
      misc$earg[[M1*ii  ]] <- .eprob
    }


    misc$imethod <- .imethod
    misc$zero <- .zero
    misc$bias.red <- .bias.red
    misc$expected <- .expected
    misc$ipstr0 <- .ipstr0
    misc$type.fitted <- .type.fitted


    misc$pobs0 <- pobs0 
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$pobs0) <- dimnames(y)
    misc$pstr0 <- pstr0
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$pstr0) <- dimnames(y)
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
                            .ipstr0 = ipstr0,
            .zero = zero,
            .expected = expected,
            .type.fitted = type.fitted,
            .bias.red = bias.red,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzigeom(x = y, prob = prob, pstr0 = pstr0, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0 ))),
  vfamily = c("zigeometric"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )
    rzigeom(nsim * length(pstr0), prob = prob, pstr0 = pstr0)
  }, list( .lprob = lprob, .lpstr0 = lpstr0,
           .eprob = eprob, .epstr0 = epstr0 ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    pstr0  <- eta2theta(eta[, c(TRUE, FALSE)], .lpstr0 , earg = .epstr0 )
    prob   <- eta2theta(eta[, c(FALSE, TRUE)], .lprob  , earg = .eprob  )


    prob0 <- prob  # P(Y == 0) from parent distribution, aka f(0)
    pobs0 <- pstr0 + (1 - pstr0) * prob0  # P(Y == 0)
    index0 <- (y == 0)

    dl.dpstr0 <- (1 - prob0) / pobs0
    dl.dpstr0[!index0] <- -1 / (1 - pstr0[!index0])

    dl.dprob <- (1 - pstr0) / pobs0
    dl.dprob[!index0]   <- 1 / prob[!index0] -
                           y[!index0] / (1 - prob[!index0])

    dpstr0.deta  <- dtheta.deta(pstr0 , .lpstr0 , earg = .epstr0 )
    dprob.deta   <- dtheta.deta(prob,   .lprob  , earg = .eprob  )

    dl.deta12 <- c(w) * cbind(dl.dpstr0 * dpstr0.deta,
                              dl.dprob  * dprob.deta)

    dl.deta12 <- dl.deta12[, interleave.VGAM(ncol(dl.deta12), M1 = M1)]
    dl.deta12
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0 ))),
  weight = eval(substitute(expression({
    if ( .expected ) {


      ned2l.dprob2 <- (1 - pstr0)^2 / pobs0 +
                      (1 - pstr0) * ((1 - prob) / prob) *
                                    (1 / prob + 1 / (1 - prob)^2)


      ned2l.dpstr0.prob <- 1 / pobs0
      ned2l.dpstr02 <- (1 - prob0) / ((1 - pstr0) * pobs0)
    } else {
      od2l.dprob2 <- ((1 - pstr0) / pobs0)^2
      od2l.dprob2[!index0] <- 1 / (prob[!index0])^2 +
                              y[!index0] / (1 - prob[!index0])^2
      od2l.dpstr0.prob <- (pobs0 + (1 - prob0) * (1 - pstr0)) / pobs0^2
      od2l.dpstr0.prob[!index0] <- 0

      od2l.dpstr02 <- ((1 - prob0) / pobs0)^2
      od2l.dpstr02[!index0] <- 1 / (1 - pstr0[!index0])^2
    }


    allvals <- if ( .expected )
                 c(c(w) * ned2l.dpstr02 * dpstr0.deta^2,
                   c(w) * ned2l.dprob2  *  dprob.deta^2,
                   c(w) * ned2l.dpstr0.prob * dprob.deta * dpstr0.deta) else
                 c(c(w) *  od2l.dpstr02 * dpstr0.deta^2,
                   c(w) *  od2l.dprob2  *  dprob.deta^2,
                   c(w) *  od2l.dpstr0.prob * dprob.deta * dpstr0.deta)
    wz <- array(allvals, dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .lprob = lprob, .lpstr0 = lpstr0,
            .eprob = eprob, .epstr0 = epstr0,
            .expected = expected ))))
}




 zigeometricff <-
  function(lprob       = "logit",
           lonempstr0  = "logit",
           type.fitted = c("mean", "prob", "pobs0", "pstr0", "onempstr0"),
           iprob = NULL,   ionempstr0  = NULL,
           imethod = 1,
           bias.red = 0.5,
           zero = "onempstr0") {


  expected <- TRUE



  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempstr0 <- as.list(substitute(lonempstr0))
  eonempstr0 <- link2list(lonempstr0)
  lonempstr0 <- attr(eonempstr0, "function.name")


  type.fitted <- match.arg(type.fitted,
                   c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]


  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")

  if (length(ionempstr0))
    if (!is.Numeric(ionempstr0, positive = TRUE) ||
        ionempstr0 >= 1)
      stop("argument 'ionempstr0' is out of range")

  if (!is.Numeric(bias.red, length.arg = 1, positive = TRUE) ||
     bias.red > 1)
    stop("argument 'bias.red' must be between 0 and 1")


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-inflated geometric distribution,\n",
            "P[Y = 0] = 1 - onempstr0 + onempstr0 * prob,\n",
            "P[Y = y] = onempstr0 * prob * (1 - prob)^y, ",
            "y = 1, 2, ...\n\n",
            "Link:     ",
            namesof("prob",       lprob,       earg = eprob ), ", ",
            namesof("onempstr0",  lonempstr0,  earg = eonempstr0), "\n",
            "Mean:     onempstr0 * (1 - prob) / prob"),
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
         parameters.names = c("prob", "onempstr0"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted  = type.fitted ))),
  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    mynames1 <- param.names("prob",      ncoly)
    mynames2 <- param.names("onempstr0", ncoly)
    predictors.names <-
      c(namesof(mynames1, .lprob      , earg = .eprob      , tag = FALSE),
        namesof(mynames2, .lonempstr0 , earg = .eonempstr0 , tag = FALSE))[
        interleave.VGAM(M1*NOS, M1 = M1)]


    if (!length(etastart)) {
      prob.init <- if ( .imethod == 3)
                       .bias.red / (1 + y + 1/8) else
                   if ( .imethod == 2)
                       .bias.red / (1 +
                   matrix(colMeans(y) + 1/8,
                          n, ncoly, byrow = TRUE)) else
                       .bias.red / (1 +
                   matrix(colSums(y * w) / colSums(w) + 1/8,
                          n, ncoly, byrow = TRUE))

      prob.init <- if (length( .iprob )) {
        matrix( .iprob , n, ncoly, byrow = TRUE)
      } else {
        prob.init  # Already a matrix
      }


      prob0.est <- psze.init <- matrix(0, n, NOS)
      for (jlocal in 1:NOS) {
        prob0.est[, jlocal] <-
          sum(w[y[, jlocal] == 0, jlocal]) / sum(w[, jlocal])
        psze.init[, jlocal] <- if ( .imethod == 3)
                         prob0.est[, jlocal] / 2 else
                     if ( .imethod == 1)
                         pmax(0.05, (prob0.est[, jlocal] -
                                     median(prob.init[, jlocal]))) else
                         prob0.est[, jlocal] / 5
      }
      psze.init <- if (length( .ionempstr0 )) {
        matrix( 1 - .ionempstr0 , n, ncoly, byrow = TRUE)
      } else {
        psze.init # Already a matrix
      }



      etastart <-
        cbind(theta2eta(    prob.init, .lprob      , earg = .eprob      ),
              theta2eta(1 - psze.init, .lonempstr0 , earg = .eonempstr0 ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0,
            .iprob = iprob, .ionempstr0 = ionempstr0,
            .type.fitted = type.fitted,
            .bias.red = bias.red,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    prob      <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                           earg = .eprob  )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )

    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs0", "pstr0", "onempstr0"))[1]

    ans <- switch(type.fitted,
                  "mean"      = onempstr0 * (1 - prob) / prob,
                  "prob"      = prob,
                  "pobs0"     = 1 - onempstr0 + onempstr0 * prob,  # P(Y=0)
                  "pstr0"     = 1 - onempstr0,
                  "onempstr0" =     onempstr0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lonempstr0 = lonempstr0,
           .eprob = eprob, .eonempstr0 = eonempstr0,
           .type.fitted = type.fitted ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lprob  , len = NOS),
                    rep( .lonempstr0 , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names


    misc$earg <- vector("list", M1 * NOS)
    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .eprob
      misc$earg[[M1*ii  ]] <- .eonempstr0
    }


    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$bias.red <- .bias.red
    misc$expected <- .expected
    misc$ionempstr0   <- .ionempstr0


    misc$pobs0 <- pobs0 
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$pobs0) <- dimnames(y)
    misc$onempstr0 <- onempstr0
    if (length(dimnames(y)[[2]]) > 0)
      dimnames(misc$onempstr0) <- dimnames(y)
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0,
                            .ionempstr0 = ionempstr0,
            .zero = zero,
            .expected = expected,
            .bias.red = bias.red,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    prob       <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                            earg = .eprob )
    onempstr0  <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                            earg = .eonempstr0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzigeom(x = y, prob = prob, pstr0 = 1 - onempstr0,
                       log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lonempstr0 = lonempstr0,
           .eprob = eprob, .eonempstr0 = eonempstr0 ))),
  vfamily = c("zigeometricff"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    prob       <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                            earg = .eprob )
    onempstr0  <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                            earg = .eonempstr0 )
    rzigeom(nsim * length(onempstr0), prob = prob, pstr0 = 1 - onempstr0)
  }, list( .lprob = lprob, .lonempstr0 = lonempstr0,
           .eprob = eprob, .eonempstr0 = eonempstr0 ))),





  deriv = eval(substitute(expression({
    M1 <- 2
    prob      <- eta2theta(eta[, c(TRUE, FALSE)], .lprob      ,
                           earg = .eprob  )
    onempstr0 <- eta2theta(eta[, c(FALSE, TRUE)], .lonempstr0 ,
                           earg = .eonempstr0 )


    prob0 <- prob  # P(Y == 0) from the parent distribution
    pobs0 <- 1 - onempstr0 + (onempstr0) * prob0  # P(Y == 0)
    index0 <- (y == 0)


    dl.donempstr0 <- -(1 - prob0) / pobs0  # zz
    dl.donempstr0[!index0] <-  1 / (onempstr0[!index0])  # zz

    dl.dprob <- (onempstr0) / pobs0
    dl.dprob[!index0]   <- 1 / prob[!index0] -
                           y[!index0] / (1 - prob[!index0])

    dprob.deta       <- dtheta.deta(prob      , .lprob      ,
                                    earg = .eprob )
    donempstr0.deta  <- dtheta.deta(onempstr0 , .lonempstr0 ,
                                    earg = .eonempstr0 )

    dl.deta12 <- c(w) * cbind(dl.dprob      * dprob.deta,
                              dl.donempstr0 *  donempstr0.deta)

    dl.deta12 <- dl.deta12[, interleave.VGAM(ncol(dl.deta12), M1 = M1)]
    dl.deta12
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0 ))),
  weight = eval(substitute(expression({
    if ( .expected ) {

      ned2l.dprob2 <- (onempstr0)^2 / pobs0 +
                      (onempstr0) * ((1 - prob) / prob) *
                                    (1 / prob + 1 / (1 - prob)^2)


      ned2l.donempstr0.prob <- -1 / pobs0
      ned2l.donempstr02 <- (1 - prob0) / ((    onempstr0) * pobs0)
    } else {
      od2l.dprob2 <- ((    onempstr0) / pobs0)^2
      od2l.dprob2[!index0] <- 1 / (prob[!index0])^2 +
                              y[!index0] / (1 - prob[!index0])^2
      od2l.donempstr0.prob <- -(pobs0 + (1 - prob0) * (onempstr0)) / pobs0^2
      od2l.donempstr0.prob[!index0] <- 0

      od2l.donempstr02 <- ((1 - prob0) / pobs0)^2
      od2l.donempstr02[!index0] <- 1 / (    onempstr0[!index0])^2
    }


    allvals <- if ( .expected )
                 c(c(w) * ned2l.dprob2  *  dprob.deta^2,
                   c(w) * ned2l.donempstr02 * donempstr0.deta^2,
                   c(w) * ned2l.donempstr0.prob * dprob.deta *
                                                  donempstr0.deta) else
                 c(c(w) *  od2l.dprob2  *  dprob.deta^2,
                   c(w) *  od2l.donempstr02 * donempstr0.deta^2,
                   c(w) *  od2l.donempstr0.prob * dprob.deta *
                                                  donempstr0.deta)
    wz <- array(allvals, dim = c(n, M / M1, 3))
    wz <- arwz2wz(wz, M = M, M1 = M1)


    wz
  }), list( .lprob = lprob, .lonempstr0 = lonempstr0,
            .eprob = eprob, .eonempstr0 = eonempstr0,
            .expected = expected ))))
}




dzageom <- function(x, prob, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(prob), length(pobs0))
  if (length(x)      != LLL) x      <- rep(x,     len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,  len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0, len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                   dposgeom(x[!index0],
                            prob = prob[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1-pobs0[!index0]) *
                   dposgeom(x[!index0],
                            prob = prob[!index0])
  }
  ans
}



pzageom <- function(q, prob, pobs0 = 0) {

  LLL <- max(length(q), length(prob), length(pobs0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans[q >  0] <- pobs0[q > 0] +
                (1 - pobs0[q > 0]) *
                pposgeom(q[q > 0], prob = prob[q > 0])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)

  ans
}


qzageom <- function(p, prob, pobs0 = 0) {

  LLL <- max(length(p), length(prob), length(pobs0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans <- p
  ind4 <- (p > pobs0)
  ans[!ind4] <- 0.0
  ans[ ind4] <- qposgeom((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         prob = prob[ind4])
  ans
}


rzageom <- function(n, prob, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n

  ans <- rposgeom(use.n, prob)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}










dzabinom <- function(x, size, prob, pobs0 = 0, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)

  LLL <- max(length(x), length(size), length(prob), length(pobs0))
  if (length(x)      != LLL) x      <- rep(x,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")
  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(pobs0[index0])
    ans[!index0] <- log1p(-pobs0[!index0]) +
                   dposbinom(x[!index0], size = size[!index0],
                             prob = prob[!index0], log = TRUE)
  } else {
    ans[ index0] <- pobs0[index0]
    ans[!index0] <- (1-pobs0[!index0]) *
                   dposbinom(x[!index0], size = size[!index0],
                             prob = prob[!index0])
  }
  ans
}



pzabinom <- function(q, size, prob, pobs0 = 0) {

  LLL <- max(length(q), length(size), length(prob), length(pobs0))
  if (length(q)      != LLL) q      <- rep(q,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);
  ans <- rep(0.0, len = LLL)
  if (!is.Numeric(pobs0) ||
      any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans[q >  0] <- pobs0[q > 0] +
                (1 - pobs0[q > 0]) *
                pposbinom(q[q > 0], size = size[q > 0], prob = prob[q > 0])
  ans[q <  0] <- 0
  ans[q == 0] <- pobs0[q == 0]

  ans <- pmax(0, ans)
  ans <- pmin(1, ans)

  ans
}


qzabinom <- function(p, size, prob, pobs0 = 0) {

  LLL <- max(length(p), length(size), length(prob), length(pobs0))
  if (length(p)      != LLL) p      <- rep(p,      len = LLL);
  if (length(size)   != LLL) size   <- rep(size,   len = LLL);
  if (length(prob)   != LLL) prob   <- rep(prob,   len = LLL);
  if (length(pobs0)  != LLL) pobs0  <- rep(pobs0,  len = LLL);

  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be in [0,1]")

  ans <- p
  ind4 <- (p > pobs0)
  ans[!ind4] <- 0.0
  ans[ ind4] <- qposbinom((p[ind4] - pobs0[ind4]) / (1 - pobs0[ind4]),
                         size = size[ind4],
                         prob = prob[ind4])
  ans
}


rzabinom <- function(n, size, prob, pobs0 = 0) {
  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n

  ans <- rposbinom(use.n, size, prob)
  if (length(pobs0) != use.n)
    pobs0 <- rep(pobs0, len = use.n)
  if (!is.Numeric(pobs0) || any(pobs0 < 0) || any(pobs0 > 1))
    stop("argument 'pobs0' must be between 0 and 1 inclusive")
  ifelse(runif(use.n) < pobs0, 0, ans)
}






 zabinomial <-
  function(lpobs0 = "logit",
           lprob  = "logit",
           type.fitted = c("mean", "prob", "pobs0"),
           ipobs0 = NULL, iprob = NULL,
           imethod = 1,
           zero = NULL  # Was zero = 2 prior to 20130917
          ) {



  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")


  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "pobs0"))[1]

  if (length(ipobs0))
    if (!is.Numeric(ipobs0, positive = TRUE) ||
        ipobs0 >= 1)
      stop("argument 'ipobs0' is out of range")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-altered binomial distribution ",
            "(Bernoulli and positive-binomial conditional model)\n\n",
            "P[Y = 0] = pobs0,\n",
            "P[Y = y] = (1 - pobs0) * dposbinom(x = y, size, prob), ",
            "y = 1, 2, ..., size,\n\n",
            "Link:     ",
            namesof("pobs0",   lpobs0, earg = epobs0), ", ",
            namesof("prob" ,   lprob,  earg = eprob),  "\n",
            "Mean:     (1 - pobs0) * prob / (1 - (1 - prob)^size)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("pobs0", "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


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
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }
    if (!all(w == 1))
      extra$new.w <- w


    y <- as.matrix(y)
    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy <- dimnames(y)
    extra$type.fitted      <- .type.fitted


    predictors.names <-
        c(namesof("pobs0", .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof("prob" , .lprob  , earg = .eprob  , tag = FALSE))
          


    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi.init <- if (length( .ipobs0 )) .ipobs0 else {
        prob0.est <- sum(Size[y == 0]) / sum(Size)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^Size) / (1 - (1 - mustart)^Size)
        } else
        if ( .imethod == 2) {
          prob0.est
        } else {
          prob0.est * 0.5
        }
    }

    phi.init[phi.init <= -0.10] <- 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] <- 0.05  # Last resort
    phi.init[phi.init >=  0.99] <- 0.95  # Last resort




    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(phi.init, .lpobs0, earg = .epobs0 ),
              theta2eta( mustart, .lprob,  earg = .eprob  ))
              

      mustart <- NULL
    }
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .iprob = iprob, .ipobs0 = ipobs0,
            .imethod = imethod,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "prob", "pobs0"))[1]
    
    phi0  <- eta2theta(eta[, 1], .lpobs0, earg = .epobs0 )
    prob  <- eta2theta(eta[, 2], .lprob,  earg = .eprob  )
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) * prob / (1 - (1 - prob)^Size),
                  "prob"      = prob,
                  "pobs0"     = phi0)  # P(Y=0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),

  last = eval(substitute(expression({
    misc$link <-    c(prob = .lprob, pobs0 = .lpobs0 )
    misc$earg <- list(prob = .eprob, pobs0 = .epobs0 )

    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$expected <- TRUE
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0,
            .zero = zero,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    pobs0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob  <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        orig.w * dzabinom(x = round(y * Size), size = Size,
                          prob = prob, pobs0 = pobs0,
                          log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lpobs0 = lpobs0,
           .eprob = eprob, .epobs0 = epobs0 ))),
  vfamily = c("zabinomial"),

  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- if (length(extra$NOS)) extra$NOS else 1

    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi0 <- eta2theta(eta[, 1], .lpobs0 , earg = .epobs0 )
    prob <- eta2theta(eta[, 2], .lprob  , earg = .eprob  )

    dphi0.deta <- dtheta.deta(phi0, .lpobs0, earg = .epobs0 )
    dprob.deta <- dtheta.deta(prob, .lprob , earg = .eprob  )

    df0.dprob   <- -Size *              (1 -  prob)^(Size - 1)
    df02.dprob2 <-  Size * (Size - 1) * (1 -  prob)^(Size - 2)
    prob0  <- (1 -  prob)^(Size)
    oneminusf0  <- 1 - prob0


    dl.dphi0 <- -1 / (1 - phi0)
    dl.dprob <-  c(w)      * (y / prob - (1 - y) / (1 - prob)) +
                 c(orig.w) * df0.dprob / oneminusf0


    dl.dphi0[y == 0] <- 1 / phi0[y == 0]  # Do it in one line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dprob[skip[, spp.], spp.] <- 0
    }


    ans <- cbind(c(orig.w) * dl.dphi0 * dphi0.deta,
                             dl.dprob * dprob.deta)
                 
                 
    ans
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))),


  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M1)

    usualmeanY <-  prob
    meanY <- (1 - phi0) * usualmeanY / oneminusf0


    term1 <-  c(Size) * (meanY /      prob^2 -
                         meanY / (1 - prob)^2) +
             c(Size) * (1 - phi0) / (1 - prob)^2

    term2 <-  -(1 - phi0) * df02.dprob2 / oneminusf0
    term3 <-  -(1 - phi0) * (df0.dprob  / oneminusf0)^2
    ned2l.dprob2 <- term1 + term2 + term3
    wz[, iam(2, 2, M)] <- ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
      tmp100
    } else {
      (dphi0.deta^2) / tmp100
    }
    wz[, iam(1, 1, M)] <- tmp200


    c(orig.w) * wz
  }), list( .lprob = lprob, .lpobs0 = lpobs0,
            .eprob = eprob, .epobs0 = epobs0 ))))
}





 zabinomialff <-
  function(lprob  = "logit",
           lonempobs0 = "logit",
           type.fitted = c("mean", "prob", "pobs0", "onempobs0"),
           iprob = NULL, ionempobs0 = NULL,
           imethod = 1,
           zero = "onempobs0") {


  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")


  type.fitted <- match.arg(type.fitted,
                   c("mean", "prob", "pobs0", "onempobs0"))[1]

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
      iprob >= 1)
    stop("argument 'iprob' is out of range")
  if (length(ionempobs0))
    if (!is.Numeric(ionempobs0, positive = TRUE) ||
        ionempobs0 >= 1)
      stop("argument 'ionempobs0' is out of range")



  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")


  new("vglmff",
  blurb = c("Zero-altered binomial distribution ",
            "(Bernoulli and positive-binomial conditional model)\n\n",
            "P[Y = 0] = 1 - onempobs0,\n",
            "P[Y = y] = onempobs0 * dposbinom(x = y, size, prob), ",
            "y = 1, 2, ..., size,\n\n",
            "Link:     ",
            namesof("prob"     , lprob     , earg = eprob     ), ", ",
            namesof("onempobs0", lonempobs0, earg = eonempobs0), "\n",
            "Mean:     onempobs0 * prob / (1 - (1 - prob)^size)"),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 2)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         Q1 = NA,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("prob", "onempobs0"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    if (!all(w == 1))
      extra$orig.w <- w



    if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1]
      nn <- rep(1, n)
      if (!all(y >= 0 & y <= 1))
        stop("response values must be in [0, 1]")
      if (!length(mustart) && !length(etastart))
        mustart <- (0.5 + w * y) / (1.0 + w)


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
           "or a factor ",
           "(first level = fail, other levels = success),\n",
           "or a 2-column matrix where col 1 is the no. of ",
           "successes and col 2 is the no. of failures")
    }
    if (!all(w == 1))
      extra$new.w <- w


    y <- as.matrix(y)
    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted


    predictors.names <-
    c(namesof("prob"     , .lprob      , earg = .eprob      , tag = FALSE),
      namesof("onempobs0", .lonempobs0 , earg = .eonempobs0 , tag = FALSE))


    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    phi.init <- if (length( .ionempobs0 )) 1 - .ionempobs0 else {
        prob0.est <- sum(Size[y == 0]) / sum(Size)
        if ( .imethod == 1) {
          (prob0.est - (1 - mustart)^Size) / (1 - (1 - mustart)^Size)
        } else
        if ( .imethod == 2) {
          prob0.est
        } else {
          prob0.est * 0.5
        }
    }

    phi.init[phi.init <= -0.10] <- 0.50  # Lots of sample variation
    phi.init[phi.init <=  0.01] <- 0.05  # Last resort
    phi.init[phi.init >=  0.99] <- 0.95  # Last resort




    if (!length(etastart)) {
      etastart <-
        cbind(theta2eta(     mustart, .lprob      , earg = .eprob      ),
              theta2eta(1 - phi.init, .lonempobs0 , earg = .eonempobs0 ))

      mustart <- NULL
    }
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0,
            .iprob = iprob, .ionempobs0 = ionempobs0,
            .imethod = imethod,
            .type.fitted = type.fitted ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs0", "onempobs0"))[1]
    
    prob      <- eta2theta(eta[, 1], .lprob      , earg = .eprob  )
    onempobs0 <- eta2theta(eta[, 2], .lonempobs0 , earg = .eonempobs0 )
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    ans <- switch(type.fitted,
                  "mean"      = onempobs0 * prob / (1 - (1 - prob)^Size),
                  "prob"      = prob,
                  "pobs0"     = 1 - onempobs0,  # P(Y=0)
                  "onempobs0" =     onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lprob = lprob, .lonempobs0 = lonempobs0,
           .eprob = eprob, .eonempobs0 = eonempobs0 ))),

  last = eval(substitute(expression({
    misc$link <-    c(prob = .lprob, onempobs0 = .lonempobs0 )
    misc$earg <- list(prob = .eprob, onempobs0 = .eonempobs0 )

    misc$imethod  <- .imethod
    misc$zero     <- .zero
    misc$expected <- TRUE
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0,
            .zero = zero,
            .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w
    prob      <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempobs0 <- eta2theta(eta[, 2], .lonempobs0 , earg = .eonempobs0 )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        orig.w * dzabinom(x = round(y * Size), size = Size,
                          prob = prob, pobs0 = 1 - onempobs0,
                          log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lprob = lprob, .lonempobs0 = lonempobs0,
           .eprob = eprob, .eonempobs0 = eonempobs0 ))),
  vfamily = c("zabinomialff"),

  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- if (length(extra$NOS)) extra$NOS else 1

    orig.w <- if (length(extra$orig.w)) extra$orig.w else 1
    new.w  <- if (length(extra$new.w))  extra$new.w  else 1
    Size <- new.w / orig.w

    prob      <- eta2theta(eta[, 1], .lprob      , earg = .eprob      )
    onempobs0 <- eta2theta(eta[, 2], .lonempobs0 , earg = .eonempobs0 )
    phi0 <- 1 - onempobs0

    dprob.deta      <- dtheta.deta(prob     , .lprob      ,
                                   earg = .eprob      )
    donempobs0.deta <- dtheta.deta(onempobs0, .lonempobs0 ,
                                   earg = .eonempobs0 )

    df0.dprob   <- -Size *              (1 -  prob)^(Size - 1)
    df02.dprob2 <-  Size * (Size - 1) * (1 -  prob)^(Size - 2)
    prob0  <- (1 -  prob)^(Size)
    oneminusf0  <- 1 - prob0


    dl.dprob <-  c(w)      * (y / prob - (1 - y) / (1 - prob)) +
                 c(orig.w) * df0.dprob / oneminusf0
    dl.donempobs0 <- +1 / (onempobs0)


    dl.donempobs0[y == 0] <-
      -1 / (1 - onempobs0[y == 0])  # Do it in 1 line
    skip <- extra$skip.these
    for (spp. in 1:NOS) {
      dl.dprob[skip[, spp.], spp.] <- 0
    }


    ans <- cbind(            dl.dprob      * dprob.deta,
                 c(orig.w) * dl.donempobs0 * donempobs0.deta)
                 
    ans
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0 ))),


  weight = eval(substitute(expression({
    wz <- matrix(0.0, n, M1)

    usualmeanY <-  prob
    meanY <- (1 - phi0) * usualmeanY / oneminusf0


    term1 <-  c(Size) * (meanY /      prob^2 -
                         meanY / (1 - prob)^2) +
             c(Size) * (1 - phi0) / (1 - prob)^2

    term2 <-  -(1 - phi0) * df02.dprob2 / oneminusf0
    term3 <-  -(1 - phi0) * (df0.dprob  / oneminusf0)^2
    ned2l.dprob2 <- term1 + term2 + term3
    wz[, iam(1, 1, M)] <- ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if (FALSE &&
                  .lonempobs0 == "logit" &&
                  is.empty.list( .eonempobs0 )) {
      tmp100
    } else {
      (donempobs0.deta^2) / tmp100
    }
    wz[, iam(2, 2, M)] <- tmp200


    c(orig.w) * wz
  }), list( .lprob = lprob, .lonempobs0 = lonempobs0,
            .eprob = eprob, .eonempobs0 = eonempobs0 ))))
}






 zageometric <-
    function(lpobs0 = "logit", lprob = "logit",
             type.fitted = c("mean", "prob", "pobs0", "onempobs0"),
             imethod = 1,
             ipobs0 = NULL, iprob = NULL,
             zero = NULL) {



  lpobs0 <- as.list(substitute(lpobs0))
  epobs0 <- link2list(lpobs0)
  lpobs0 <- attr(epobs0, "function.name")

  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "prob", "pobs0", "onempobs0"))[1]


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")
  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
       max(iprob) >= 1)
    stop("argument 'iprob' out of range")
  if (length(ipobs0))
    if (!is.Numeric(ipobs0, positive = TRUE) ||
       max(ipobs0) >= 1)
      stop("argument 'ipobs0' out of range")


  new("vglmff",
  blurb = c("Zero-altered geometric ",
            "(Bernoulli and positive-geometric conditional model)\n\n",
            "Links:    ",
            namesof("pobs0", lpobs0, earg = epobs0, tag = FALSE), ", ",
            namesof("prob" , lprob , earg = eprob , tag = FALSE), "\n",
            "Mean:     (1 - pobs0) / prob"),

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
         parameters.names = c("pobs0", "prob"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy <- dimnames(y)
    extra$type.fitted      <- .type.fitted

    
    mynames1 <- param.names("pobs0", ncoly)
    mynames2 <- param.names("prob",  ncoly)
    predictors.names <-
        c(namesof(mynames1, .lpobs0 , earg = .epobs0 , tag = FALSE),
          namesof(mynames2, .lprob  , earg = .eprob  , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {

      foo <- function(x) mean(as.numeric(x == 0))
      phi0.init <- matrix(apply(y, 2, foo), n, ncoly, byrow = TRUE)
      if (length( .ipobs0 ))
        phi0.init <- matrix( .ipobs0 , n, ncoly, byrow = TRUE)


      prob.init <-
        if ( .imethod == 2)
          1 / (1 + y + 1/16) else
        if ( .imethod == 1)
          (1 - phi0.init) / (1 +
          matrix(colSums(y * w) / colSums(w) + 1/16,
                 n, ncoly, byrow = TRUE)) else
          (1 - phi0.init) / (1 +
          matrix(apply(y, 2, median), n, ncoly, byrow = TRUE) + 1/16)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, ncoly, byrow = TRUE)



      etastart <- cbind(theta2eta(phi0.init, .lpobs0 , earg = .epobs0 ),
                       theta2eta(prob.init, .lprob ,  earg = .eprob ))
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob,
            .ipobs0 = ipobs0, .iprob = iprob,
            .imethod = imethod,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                             c("mean", "prob", "pobs0", "onempobs0"))[1]
    M1 <- 2
    NOS <- ncol(eta) / M1

    phi0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                             .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                             .lprob  , earg = .eprob ))


    ans <- switch(type.fitted,
                  "mean"      = (1 - phi0) / prob,
                  "prob"      = prob,
                  "pobs0"     =      phi0,  # P(Y=0)
                  "onempobs0" =  1 - phi0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lpobs0 , len = NOS),
                    rep( .lprob  , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names

    misc$earg <- vector("list", M1 * NOS)

    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .epobs0
      misc$earg[[M1*ii  ]] <- .eprob
    }


    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$ipobs0  <- .ipobs0
    misc$iprob   <- .iprob
    misc$multipleResponses <- TRUE
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob,
            .ipobs0 = ipobs0, .iprob = iprob,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    phi0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                            .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                            .lprob  , earg = .eprob  ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzageom(x = y, pobs0 = phi0, prob = prob, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),
  vfamily = c("zageometric"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    phi0 <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                            .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                            .lprob  , earg = .eprob  ))
    rzageom(nsim * length(prob), prob = prob, pobs0 = phi0)
  }, list( .lpobs0 = lpobs0, .lprob = lprob,
           .epobs0 = epobs0, .eprob = eprob ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    phi0 <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                            .lpobs0 , earg = .epobs0 ))
    prob <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                            .lprob  , earg = .eprob  ))


    dl.dprob <-  1 / prob - (y - 1) / (1 - prob)
    dl.dphi0 <- -1 / (1 - phi0)


    for (spp. in 1:NOS) {
      dl.dphi0[skip[, spp.], spp.] <- 1 / phi0[skip[, spp.], spp.]
      dl.dprob[skip[, spp.], spp.] <- 0
    }
    dphi0.deta <- dtheta.deta(phi0, .lpobs0 , earg = .epobs0 )
    dprob.deta <- dtheta.deta(prob, .lprob  , earg = .eprob  )


    ans <- c(w) * cbind(dl.dphi0 * dphi0.deta,
                        dl.dprob * dprob.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lpobs0 = lpobs0, .lprob = lprob,
            .epobs0 = epobs0, .eprob = eprob ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1*NOS)


    ned2l.dprob2 <- (1 - phi0) / (prob^2 * (1 - prob))

    wz[, NOS+(1:NOS)] <- c(w) * ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( .lpobs0 == "logit" && is.empty.list( .epobs0 )) {
      cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (dphi0.deta^2) / tmp100)
    }
    wz[, 1:NOS] <-  tmp200


    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]


    wz
  }), list( .lpobs0 = lpobs0,
            .epobs0 = epobs0 ))))
}  # End of zageometric




 zageometricff <-
    function(lprob = "logit", lonempobs0 = "logit",
             type.fitted = c("mean", "prob", "pobs0", "onempobs0"),
             imethod = 1,
             iprob = NULL, ionempobs0 = NULL,
             zero = "onempobs0") {


  lprob <- as.list(substitute(lprob))
  eprob <- link2list(lprob)
  lprob <- attr(eprob, "function.name")

  lonempobs0 <- as.list(substitute(lonempobs0))
  eonempobs0 <- link2list(lonempobs0)
  lonempobs0 <- attr(eonempobs0, "function.name")

  type.fitted <- match.arg(type.fitted,
                   c("mean", "prob", "pobs0", "onempobs0"))[1]


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
     imethod > 3)
    stop("argument 'imethod' must be 1 or 2 or 3")

  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
       max(iprob) >= 1)
    stop("argument 'iprob' out of range")

  if (length(ionempobs0))
    if (!is.Numeric(ionempobs0, positive = TRUE) ||
       max(ionempobs0) >= 1)
      stop("argument 'ionempobs0' out of range")


  new("vglmff",
  blurb = c("Zero-altered geometric ",
            "(Bernoulli and positive-geometric conditional model)\n\n",
            "Links:    ",
            namesof("prob"     , lprob     , earg = eprob     , tag = FALSE), ", ",
            namesof("onempobs0", lonempobs0, earg = eonempobs0, tag = FALSE), "\n",
            "Mean:     onempobs0 / prob"),

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
         parameters.names = c("prob", "onempobs0"),
         type.fitted  = .type.fitted ,
         zero = .zero )
  }, list( .zero = zero,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.nonnegative.y = TRUE,
              Is.integer.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y




    extra$y0 <- y0 <- ifelse(y == 0, 1, 0)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$skip.these <- skip.these <- matrix(as.logical(y0), n, NOS)

    extra$dimnamesy   <- dimnames(y)
    extra$type.fitted <- .type.fitted

    
    mynames1 <- param.names("prob",       ncoly)
    mynames2 <- param.names("onempobs0",  ncoly)
    predictors.names <-
        c(namesof(mynames1, .lprob      , earg = .eprob      , tag = FALSE),
          namesof(mynames2, .lonempobs0 , earg = .eonempobs0 , tag = FALSE))[
          interleave.VGAM(M1*NOS, M1 = M1)]

    if (!length(etastart)) {

      foo <- function(x) mean(as.numeric(x == 0))
      phi0.init <- matrix(apply(y, 2, foo), n, ncoly, byrow = TRUE)
      if (length( .ionempobs0 ))
        phi0.init <- matrix( 1 - .ionempobs0 , n, ncoly, byrow = TRUE)


      prob.init <-
        if ( .imethod == 2)
          1 / (1 + y + 1/16) else
        if ( .imethod == 1)
          (1 - phi0.init) / (1 +
          matrix(colSums(y * w) / colSums(w) + 1/16,
                 n, ncoly, byrow = TRUE)) else
          (1 - phi0.init) / (1 +
          matrix(apply(y, 2, median), n, ncoly, byrow = TRUE) + 1/16)


      if (length( .iprob ))
        prob.init <- matrix( .iprob , n, ncoly, byrow = TRUE)



      etastart <-
        cbind(theta2eta(    prob.init, .lprob      , earg = .eprob      ),
              theta2eta(1 - phi0.init, .lonempobs0 , earg = .eonempobs0 ))
                        
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M1 = M1)]
    }
  }), list( .lonempobs0 = lonempobs0, .lprob = lprob,
            .eonempobs0 = eonempobs0, .eprob = eprob,
            .ionempobs0 = ionempobs0, .iprob = iprob,
            .imethod = imethod,
            .type.fitted = type.fitted ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "prob", "pobs0", "onempobs0"))[1]

    NOS <- extra$NOS
    M1 <- 2

    prob      <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                 .lprob  , earg = .eprob ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))


    ans <- switch(type.fitted,
                  "mean"      =  onempobs0 / prob,
                  "prob"      =  prob,
                  "pobs0"     =  1 - onempobs0,  # P(Y=0)
                  "onempobs0" =      onempobs0)  # P(Y>0)
    if (length(extra$dimnamesy) &&
        is.matrix(ans) &&
        length(extra$dimnamesy[[2]]) == ncol(ans) &&
        length(extra$dimnamesy[[2]]) > 0) {
      dimnames(ans) <- extra$dimnamesy
    } else
    if (NCOL(ans) == 1 &&
        is.matrix(ans)) {
      colnames(ans) <- NULL
    }
    ans
  }, list( .lonempobs0 = lonempobs0, .lprob = lprob,
           .eonempobs0 = eonempobs0, .eprob = eprob ))),
  last = eval(substitute(expression({
    temp.names <- c(rep( .lprob      , len = NOS),
                    rep( .lonempobs0 , len = NOS))
    temp.names <- temp.names[interleave.VGAM(M1*NOS, M1 = M1)]
    misc$link  <- temp.names

    misc$earg <- vector("list", M1 * NOS)

    names(misc$link) <-
    names(misc$earg) <-
        c(mynames1, mynames2)[interleave.VGAM(M1*NOS, M1 = M1)]

    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .eprob
      misc$earg[[M1*ii  ]] <- .eonempobs0
    }


    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$ionempobs0  <- .ionempobs0
    misc$iprob   <- .iprob
    misc$multipleResponses <- TRUE
  }), list( .lonempobs0 = lonempobs0, .lprob = lprob,
            .eonempobs0 = eonempobs0, .eprob = eprob,
            .ionempobs0 = ionempobs0, .iprob = iprob,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    NOS <- extra$NOS
    M1 <- 2

    prob      <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                                 .lprob      , earg = .eprob      ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dzageom(x = y, pobs0 = 1 - onempobs0, prob = prob,
                       log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lonempobs0 = lonempobs0, .lprob = lprob,
           .eonempobs0 = eonempobs0, .eprob = eprob ))),
  vfamily = c("zageometricff"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    onempobs0 <- cbind(eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                                 .lonempobs0 , earg = .eonempobs0 ))
    prob      <- cbind(eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                                 .lprob      , earg = .eprob      ))
    rzageom(nsim * length(prob), pobs0 = 1 - onempobs0, prob = prob)
  }, list( .lonempobs0 = lonempobs0, .lprob = lprob,
           .eonempobs0 = eonempobs0, .eprob = eprob ))),




  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- ncol(eta) / M1  # extra$NOS
    y0 <- extra$y0
    skip <- extra$skip.these

    prob      <- cbind(eta2theta(eta[, M1*(1:NOS)-1, drop = FALSE],
                       .lprob      , earg = .eprob      ))
    onempobs0 <- cbind(eta2theta(eta[, M1*(1:NOS)-0, drop = FALSE],
                       .lonempobs0 , earg = .eonempobs0 ))
    pobs0 <- 1 - onempobs0


    dl.dprob      <-  1 / prob - (y - 1) / (1 - prob)
    dl.donempobs0 <- +1 / (onempobs0)


    for (spp. in 1:NOS) {
      dl.donempobs0[skip[, spp.], spp.] <- -1 / pobs0[skip[, spp.], spp.]
      dl.dprob[skip[, spp.], spp.] <- 0
    }
    dprob.deta      <- dtheta.deta(prob,      .lprob  , earg = .eprob  )
    donempobs0.deta <- dtheta.deta(onempobs0, .lonempobs0 ,
                                   earg = .eonempobs0 )


    ans <- c(w) * cbind(dl.dprob      * dprob.deta,
                        dl.donempobs0 * donempobs0.deta)
    ans <- ans[, interleave.VGAM(ncol(ans), M1 = M1)]
    ans
  }), list( .lonempobs0 = lonempobs0, .lprob = lprob,
            .eonempobs0 = eonempobs0, .eprob = eprob ))),
  weight = eval(substitute(expression({

    wz <- matrix(0.0, n, M1*NOS)


    ned2l.dprob2 <- (1 - pobs0) / (prob^2 * (1 - prob))

    wz[, (1:NOS)] <- c(w) * ned2l.dprob2 * dprob.deta^2


    mu.phi0 <- pobs0  # phi0
    tmp100 <- mu.phi0 * (1.0 - mu.phi0)
    tmp200 <- if ( FALSE &&
                  .lonempobs0 == "logit" &&
                  is.empty.list( .eonempobs0 )) {

      cbind(c(w) * tmp100)
    } else {
      cbind(c(w) * (donempobs0.deta^2) / tmp100)
    }
    wz[, NOS+(1:NOS)] <- tmp200


    wz <- wz[, interleave.VGAM(ncol(wz), M1 = M1)]


    wz
  }), list( .lonempobs0 = lonempobs0,
            .eonempobs0 = eonempobs0 ))))
}  # End of zageometricff






