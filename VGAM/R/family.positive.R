# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.












N.hat.posbernoulli <-
  function(eta, link, earg = list(),
           R = NULL, w = NULL,
           X.vlm = NULL, Hlist = NULL,
           extra = list(),
           model.type = c("0", "b", "t", "tb")
          ) {





  if (!is.null(w) && !all(w[1] == w))
    warning("estimate of N may be wrong when prior weights ",
            "are not all the same")

  model.type <- match.arg(model.type, c("0", "b", "t", "tb"))[1]
  if (!is.matrix(eta))
    eta <- as.matrix(eta)  # May be needed for "0"
 
  tau <-
    switch(model.type,
           "0"  = extra$tau,
           "b"  = extra$tau,
           "t"  = ncol(eta),
           "tb" = (ncol(eta) + 1) / 2)
  if (length(extra$tau) && extra$tau != tau)
    warning("variable 'tau' is mistaken")  # Checking only


  jay.index <-
    switch(model.type,
           "0"  = rep(1, length = tau),
           "b"  = rep(1, length = tau),  # Subset: 2 out of 1:2
           "t"  = 1:tau,  # All of them
           "tb" = 1:tau)  # Subset: first tau of them out of M = 2*tau-2

  prc <- eta2theta(eta[, jay.index], link, earg = earg)  # cap.probs
  prc <- as.matrix(prc)  # Might be needed for Mtb(tau=2).


 
  if (FALSE && model.type == "tb") {
    if (tau == 2)
      prc <- cbind(prc, 1 - prc)
    if (tau >  3)
      stop("cannot handle tau > 3 yet")
    jay.index <- 1:tau  # 'Restore' it coz its used below. zz??
  }
  
  QQQ <- exp(rowSums(log1p(-prc)))
  pibbeta <- exp(log1p(-QQQ))  # One.minus.QQQ
  N.hat <- sum(1 / pibbeta)  # Point estimate
  ss2 <- sum(QQQ / pibbeta^2)  # Assumes bbeta is known


  if (length(extra$p.small) &&
      any(pibbeta < extra$p.small) &&
      !extra$no.warning)
    warning("The abundance estimation for this model can be unstable")


  if (length(R)) {

    dvect <- matrix(0, length(pibbeta), ncol = ncol(X.vlm))
    M <- nrow(Hlist[[1]])
    n.lm <- nrow(X.vlm) / M  # Number of rows of the LM matrix
    dprc.deta <- dtheta.deta(prc, link, earg = earg)
    Hmatrices <- matrix(c(unlist(Hlist)), nrow = M)
    for (jay in 1:tau) {
      linpred.index <- jay.index[jay]
      Index0 <- Hmatrices[linpred.index, ] != 0
      X.lm.jay <- X.vlm[(0:(n.lm - 1)) * M + linpred.index,
                        Index0,
                        drop = FALSE]

      dvect[, Index0] <-
      dvect[, Index0] +
        (QQQ / (1-prc[, jay])) * dprc.deta[, jay] * X.lm.jay
    }


   dvect <- dvect * (-1 / pibbeta^2)
   dvect <- colSums(dvect)  # Now a vector

    ncol.X.vlm <- nrow(R)
    rinv <- diag(ncol.X.vlm)
    rinv <- backsolve(R, rinv)
    rowlen <- drop(((rinv^2) %*% rep(1, ncol.X.vlm))^0.5)
    covun <- rinv %*% t(rinv)
    vecTF <- FALSE
    for (jay in 1:tau) {
      linpred.index <- jay.index[jay]
      vecTF <- vecTF | (Hmatrices[linpred.index, ] != 0)
    }
    vecTF.index <- (1:length(vecTF))[vecTF]
    covun <- covun[vecTF.index, vecTF.index, drop = FALSE]
    dvect <- dvect[vecTF.index, drop = FALSE]
  }
 
  list(N.hat    = N.hat,
       SE.N.hat = if (length(R))
                    c(sqrt(ss2 + t(dvect) %*% covun %*% dvect)) else
                    c(sqrt(ss2))
      )
}




 aux.posbernoulli.t <-
  function(y, check.y = FALSE,
           rename = TRUE,
           name = "bei") {






  y <- as.matrix(y)
  if ((tau <- ncol(y)) == 1)
    stop("argument 'y' needs to be a matrix with at least two columns")
  if (check.y) {
    if (!all(y == 0 | y == 1 | y == 1/tau | is.na(y)))
      stop("response 'y' must contain 0s and 1s only")
  }

  zeddij <- cbind(0, t(apply(y, 1, cumsum)))  # tau + 1 columns
  zij <- (0 + (zeddij > 0))[, 1:tau]  # 0 or 1.
  if (rename) {
    colnames(zij) <- paste(name, 1:ncol(y), sep = "")
  } else {
    if (length(colnames(y)))
      colnames(zij) <- colnames(y)
  }


  cp1 <- numeric(nrow(y))
  for (jay in tau:1)
    cp1[y[, jay] > 0] <- jay
  if (any(cp1 == 0))
    warning("some individuals were never captured!")

  yr1i <- zeddij[, tau + 1] - 1
  list(cap.hist1 = zij,  # A matrix of the same dimension as 'y'
       cap1      = cp1,  # Aka ti1
       y0i       = cp1 - 1,
       yr0i      = tau - cp1 - yr1i,
       yr1i      = yr1i)
}










rposbern <-
  function(n, nTimePts = 5, pvars = length(xcoeff),
           xcoeff = c(-2, 1, 2),
           Xmatrix = NULL,  # If is.null(Xmatrix) then it is created
           cap.effect =  1,
           is.popn = FALSE,
           link = "logit",
           earg.link = FALSE) {









  use.n <- if ((length.n <- length(n)) > 1) length.n else
           if (!is.Numeric(n, integer.valued = TRUE,
                           length.arg = 1, positive = TRUE))
               stop("bad input for argument 'n'") else n
  orig.n <- use.n
  if (!is.popn)
    use.n <- 1.50 * use.n + 100 # Bigger due to rejections

  if (pvars == 0)
    stop("argument 'pvars' must be at least one")
  if (pvars > length(xcoeff))
    stop("argument 'pvars' is too high")
  

  if (earg.link) {
    earg <- link
  } else {
    link <- as.list(substitute(link))
    earg <- link2list(link)
  }
  link <- attr(earg, "function.name")


  cap.effect.orig <- cap.effect


  Ymatrix <- matrix(0, use.n, nTimePts,
                    dimnames = list(as.character(1:use.n),
                                    paste("y", 1:nTimePts, sep = "")))

  CHmatrix <- matrix(0, use.n, nTimePts,
                     dimnames = list(as.character(1:use.n),
                                     paste("ch", 1:(nTimePts  ),
                                           sep = "")))


  if (is.null(Xmatrix)) {
    Xmatrix <- cbind(x1 = rep(1.0, len = use.n))
    if (pvars > 1)
      Xmatrix <- cbind(Xmatrix,
                       matrix(runif(n = use.n * (pvars-1)),
                              use.n, pvars - 1,
                              dimnames = list(as.character(1:use.n),
                                              paste("x", 2:pvars, sep = ""))))
  }


  lin.pred.baseline <- xcoeff[1]
  if (pvars > 1)
    lin.pred.baseline <- lin.pred.baseline +
                         Xmatrix[, 2:pvars, drop = FALSE] %*%
                         xcoeff[2:pvars]
  sumrowy <- rep(0, length = use.n)
  cap.effect <- rep(cap.effect.orig, length = use.n)

  for (jlocal in 1:nTimePts) {
    CHmatrix[, jlocal] <- as.numeric(sumrowy > 0)

    caught.before.TF <- (CHmatrix[, jlocal] >  0)
    lin.pred <- lin.pred.baseline + caught.before.TF * cap.effect

    Ymatrix[, jlocal] <-
      rbinom(use.n, size = 1,
             prob = eta2theta(lin.pred, link = link, earg = earg))

    sumrowy <- sumrowy + Ymatrix[, jlocal]
  }



  index0 <- (sumrowy == 0)
  if (all(!index0))
    stop("bug in this code: cannot handle no animals being caught")
   Ymatrix <-  Ymatrix[!index0, , drop = FALSE]
   Xmatrix <-  Xmatrix[!index0, , drop = FALSE]
  CHmatrix <- CHmatrix[!index0, , drop = FALSE]




  ans <- data.frame(Ymatrix, Xmatrix, CHmatrix  # zCHmatrix,
                   )


  if (!is.popn) {
    ans <- if (nrow(ans) >= orig.n) {
      ans[1:orig.n, ]
    } else {
      rbind(ans,
            Recall(n            = orig.n - nrow(ans),
                   nTimePts     = nTimePts,
                   pvars        = pvars,
                   xcoeff       = xcoeff,
                   cap.effect   = cap.effect.orig,
                   link         = earg,
                   earg.link    = TRUE))
    }
  }

  rownames(ans) <- as.character(1:nrow(ans))

  attr(ans, "pvars")      <- pvars
  attr(ans, "nTimePts")   <- nTimePts
  attr(ans, "cap.effect") <- cap.effect.orig
  attr(ans, "is.popn")    <- is.popn
  attr(ans, "n")          <- n

  ans
}





  

dposbern <- function(x, prob, prob0 = prob, log = FALSE) {


  x     <- as.matrix(x)
  prob  <- as.matrix(prob)
  prob0 <- as.matrix(prob0)

  if (!is.logical(log.arg <- log) ||
      length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  if (ncol(x) < 2)
    stop("columns of argument 'x' should be 2 or more")


  logAA0 <- rowSums(log1p(-prob0))
  AA0 <- exp(logAA0)
  ell1 <- x * log(prob) + (1 - x) * log1p(-prob) - log1p(-AA0) / ncol(x)
  if (log.arg) ell1 else exp(ell1)
}









dposnegbin <- function(x, size, prob = NULL, munb = NULL, log = FALSE) {
  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  LLL <- max(length(x), length(prob), length(size))
  if (length(x)    != LLL) x    <- rep(x,    len = LLL)
  if (length(prob) != LLL) prob <- rep(prob, len = LLL)
  if (length(size) != LLL) size <- rep(size, len = LLL)

  ans <- dnbinom(x = x, size = size, prob = prob, log = log.arg)
  index0 <- (x == 0)

  if (log.arg) {
    ans[ index0] <- log(0.0)
    ans[!index0] <- ans[!index0] - log1p(-dnbinom(x = 0 * x[!index0],
                    size = size[!index0], prob = prob[!index0]))
  } else {
    ans[ index0] <- 0.0
    ans[!index0] <- ans[!index0] / pnbinom(q = 0 * x[!index0],
                    size = size[!index0], prob = prob[!index0],
                    lower.tail = FALSE)
  }
  ans
}



pposnegbin <- function(q, size, prob = NULL, munb = NULL) {

  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }
  L <- max(length(q), length(prob), length(size))
  if (length(q)    != L)
    q    <- rep(q,    length.out = L)
  if (length(prob) != L)
    prob <- rep(prob, length.out = L)
  if (length(size) != L)
    size <- rep(size, length.out = L)

  ifelse(q < 1, 0,
        (pnbinom(q, size = size, prob = prob) -
         dnbinom(0, size = size, prob = prob))
       / pnbinom(0, size = size, prob = prob, lower.tail = FALSE))
}



qposnegbin <- function(p, size, prob = NULL, munb = NULL) {


  if (length(munb)) {
    if (length(prob))
      stop("'prob' and 'munb' both specified")
    prob <- size / (size + munb)
  }

  ans <- qnbinom(pnbinom(q = 0, size = size, prob = prob,
                         lower.tail = FALSE) * p +
                 dnbinom(x = 0, size = size, prob = prob),
                 size = size, prob = prob)
  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- Inf
  ans
}







    EIM.posNB.specialp <- function(munb, size,
                                   y.max = NULL,  # Must be an integer
                                   cutoff.prob = 0.995,
                                   prob0, df0.dkmat, df02.dkmat2,
                                   intercept.only = FALSE,
                                   second.deriv = TRUE) {


      if (intercept.only) {
        munb        <- munb[1]
        size        <- size[1]
        prob0       <- prob0[1]
        df0.dkmat   <- df0.dkmat[1]
        df02.dkmat2 <- df02.dkmat2[1]
      }

      y.min <- 0  # Same as negbinomial() actually. A fixed constant really

      if (!is.numeric(y.max)) {
        eff.p <- sort(c(cutoff.prob, 1 - cutoff.prob))
        y.max <- max(qposnegbin(p = eff.p[2], munb = munb, size = size)) + 10
      }

      Y.mat <- if (intercept.only) y.min:y.max else
               matrix(y.min:y.max, length(munb), y.max-y.min+1, byrow = TRUE)
  neff.row <- ifelse(intercept.only, 1, nrow(Y.mat))
  neff.col <- ifelse(intercept.only, length(Y.mat), ncol(Y.mat))

      if (FALSE) {
      Y.mat2 <- Y.mat + 1
      trigg.term0 <- if (intercept.only) {
         dposnegbin(Y.mat2, size=size, munb=munb) %*% trigamma(Y.mat2+size)
      } else {
         rowSums(dposnegbin(Y.mat2, size = size, munb = munb) *
                 trigamma(Y.mat2 + size))
      }
      }


  trigg.term <- 
  if (TRUE) {
    answerC <- .C("eimpnbinomspecialp",
      as.integer(intercept.only),
      as.double(neff.row), as.double(neff.col),
      as.double(size),
      as.double(1 - pposnegbin(Y.mat, size = size, munb = munb)),
      rowsums = double(neff.row))
      answerC$rowsums
  }



      mymu <- munb / (1 - prob0)  # E(Y)
      ned2l.dk2 <- trigg.term -
         munb / (size * (size + munb)) - (mymu - munb) / (munb + size)^2

      if (second.deriv)
        ned2l.dk2 <- ned2l.dk2 - df02.dkmat2 / (1 - prob0) -
         (df0.dkmat / (1 - prob0))^2
      ned2l.dk2
    }  # end of EIM.posNB.specialp()







    EIM.posNB.speciald <- function(munb, size,
                                   y.min = 1,  # 20160201; must be an integer
                                   y.max = NULL,  # Must be an integer
                                   cutoff.prob = 0.995,
                                   prob0, df0.dkmat, df02.dkmat2,
                                   intercept.only = FALSE,
                                   second.deriv = TRUE) {


      if (intercept.only) {
        munb        <- munb[1]
        size        <- size[1]
        prob0       <- prob0[1]
        df0.dkmat   <- df0.dkmat[1]
        df02.dkmat2 <- df02.dkmat2[1]
      }

      if (!is.numeric(y.max)) {
        eff.p <- sort(c(cutoff.prob, 1 - cutoff.prob))
        y.max <- max(qposnegbin(p = eff.p[2], munb = munb, size = size)) + 10
      }

      Y.mat <- if (intercept.only) y.min:y.max else
               matrix(y.min:y.max, length(munb), y.max-y.min+1, byrow = TRUE)
      trigg.term <- if (intercept.only) {
         dposnegbin(Y.mat, size = size, munb = munb) %*% trigamma(Y.mat + size)
      } else {
         rowSums(dposnegbin(Y.mat, size = size, munb = munb) *
                 trigamma(Y.mat + size))
      }

      mymu <- munb / (1 - prob0)  # E(Y)
      ned2l.dk2 <- trigamma(size) - munb / (size * (size + munb)) -
        (mymu - munb) / (munb + size)^2 - trigg.term
      if (second.deriv)
        ned2l.dk2 <- ned2l.dk2 - df02.dkmat2 / (1 - prob0) -
         (df0.dkmat / (1 - prob0))^2
      ned2l.dk2
    }  # end of EIM.posNB.speciald()





posnegbinomial.control <- function(save.weights = TRUE, ...) {
  list(save.weights = save.weights)
}



 posnegbinomial <-
  function(
           zero = "size",
           type.fitted = c("mean", "munb", "prob0"),
           nsimEIM = 500,
           cutoff.prob = 0.999,  # higher is better for large 'size'
           eps.trig = 1e-7,
           max.support = 4000,  # 20160201; I have changed this
           max.chunk.MB = 30,  # max.memory = Inf is allowed
           lmunb = "loge", lsize = "loge",
           imethod = 1,
           imunb = NULL,
           probs.y = 0.35,
           ishrinkage = 0.95,
           isize = NULL,
           gsize.mux = exp((-12:6)/2)) {



  if (length(isize) && !is.Numeric(isize, positive = TRUE))
      stop("bad input for argument 'isize'")

  lmunb <- as.list(substitute(lmunb))
  emunb <- link2list(lmunb)
  lmunb <- attr(emunb, "function.name")

  lsize <- as.list(substitute(lsize))
  esize <- link2list(lsize)
  lsize <- attr(esize, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "munb", "prob0"))[1]


  if (!is.Numeric(eps.trig, length.arg = 1,
                  positive = TRUE) || eps.trig > 0.001)
    stop("argument 'eps.trig' must be positive and smaller in value")

  if (!is.Numeric(nsimEIM, length.arg = 1,
                  positive = TRUE, integer.valued = TRUE))
    stop("argument 'nsimEIM' must be a positive integer")
  if (nsimEIM <= 30)
    warning("argument 'nsimEIM' should be greater than 30, say")


  new("vglmff",
  blurb = c("Positive-negative binomial distribution\n\n",
            "Links:    ",
            namesof("munb", lmunb, earg = emunb ), ", ",
            namesof("size", lsize, earg = esize ), "\n",
            "Mean:     munb / (1 - (size / (size + munb))^size)"),
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
         parameters.names = c("munb", "size"),
         nsimEIM = .nsimEIM ,
         eps.trig = .eps.trig ,
         lmunb = .lmunb ,
         emunb = .emunb ,
         type.fitted  = .type.fitted ,
         zero = .zero ,
         lsize = .lsize ,
         esize = .esize )
  }, list( .lmunb = lmunb, .lsize = lsize, .isize = isize,
           .emunb = emunb, .esize = esize,
           .zero = zero, .nsimEIM = nsimEIM,
           .ishrinkage = ishrinkage, .eps.trig = eps.trig,
           .imethod = imethod,
           .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    M1 <- 2

    temp5 <-
    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.positive.y = TRUE,
              ncol.w.max = Inf,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = 1,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y



    M <- M1 * ncol(y) 
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)

    predictors.names <- c(
      namesof(param.names("munb", NOS), .lmunb , earg = .emunb , tag = FALSE),
      namesof(param.names("size", NOS), .lsize , earg = .esize , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {
      munb.init <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                           imu = .imunb , ishrinkage = .ishrinkage ,
                           probs.y = .probs.y )


      if ( is.Numeric( .isize )) {
        size.init <- matrix( .isize , nrow = n, ncol = NOS, byrow = TRUE)
      } else {
        posnegbinomial.Loglikfun <- function(kval, y, x, w, extraargs) {
          munb <- extraargs
          sum(c(w) * dposnegbin(x = y, mu = munb, size = kval, log = TRUE))
        }
        size.init <- matrix(0, nrow = n, ncol = NOS)
        for (jay in 1:NOS) {
          size.grid <- .gsize.mux * mean(munb.init[, jay])
          size.init[, jay] <-
            grid.search(size.grid,
                        objfun = posnegbinomial.Loglikfun,
                        y = y[, jay],  # x = x,
                        w = w[, jay],
                        extraargs = munb.init[, jay])
        }
      }



      etastart <-
        cbind(
              theta2eta(munb.init            , .lmunb , earg = .emunb ),
              theta2eta(size.init,             .lsize , earg = .esize ))
      etastart <- etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
    }
  }), list( .lmunb = lmunb, .lsize  = lsize,
            .imunb = imunb, .isize = isize,
            .emunb = emunb, .esize  = esize, .gsize.mux = gsize.mux,
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
                     c("mean", "munb", "prob0"))[1]

    TF <- c(TRUE, FALSE)
    munb <- eta2theta(eta[,  TF, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, !TF, drop = FALSE], .lsize , earg = .esize )


    tempk <- 1 / (1 + munb / kmat)  # kmat / (kmat + munb)
    prob0  <- tempk^kmat
    oneminusf0  <- 1 - prob0

    smallval <- 1e-3  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
      prob0[big.size]  <- exp(-munb[big.size])  # The limit as kmat --> Inf
      oneminusf0[big.size] <- -expm1(-munb[big.size])
    }

    ans <- switch(type.fitted,
                  "mean"      = munb / oneminusf0,
                  "munb"      = munb,
                  "prob0"     = prob0)  # P(Y=0)
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
  }, list( .lsize = lsize, .lmunb = lmunb,
           .esize = esize, .emunb = emunb ))),
  last = eval(substitute(expression({
    temp0303 <- c(rep( .lmunb , length = NOS),
                  rep( .lsize , length = NOS))
    names(temp0303) <- c(param.names("munb", NOS),
                         param.names("size", NOS))
    temp0303  <- temp0303[interleave.VGAM(M, M1 = M1)]
    misc$link <- temp0303  # Already named

    misc$earg <- vector("list", M1*NOS)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:NOS) {
      misc$earg[[M1*ii-1]] <- .emunb
      misc$earg[[M1*ii  ]] <- .esize
    }

    misc$max.chunk.MB <- .max.chunk.MB
    misc$cutoff.prob <- .cutoff.prob
    misc$imethod <- .imethod 
    misc$nsimEIM <- .nsimEIM
    misc$expected <- TRUE
    misc$ishrinkage <- .ishrinkage
    misc$multipleResponses <- TRUE
   }), list( .lmunb = lmunb, .lsize = lsize,
             .emunb = emunb, .esize = esize,
             .cutoff.prob = cutoff.prob,
             .max.chunk.MB = max.chunk.MB,
             .ishrinkage = ishrinkage,
             .nsimEIM = nsimEIM, .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    TFvec <- c(TRUE, FALSE)
    munb <- eta2theta(eta[,  TFvec, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, !TFvec, drop = FALSE], .lsize , earg = .esize )
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <-
        c(w) * dposnegbin(x = y, size = kmat, munb = munb, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize ))),

  vfamily = c("posnegbinomial"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    munb <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                      .lmunb, earg = .emunb )
    kmat <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                      .lsize, earg = .esize )
    rposnegbin(nsim * length(munb), size = kmat, munb = munb)
  }, list( .lmunb = lmunb, .lsize = lsize,
           .emunb = emunb, .esize = esize ))),


  validparams = eval(substitute(function(eta, extra = NULL) {
    munb <- eta2theta(eta[, c(TRUE, FALSE), drop = FALSE],
                     .lmunb , earg = .emunb )
    size <- eta2theta(eta[, c(FALSE, TRUE), drop = FALSE],
                     .lsize , earg = .esize )

    smallval <- 1e-2
    ans <- all(is.finite(munb)) && all(munb > 0) &&
           all(is.finite(size)) && all(size > 0) &&
           (overdispersion <- all(munb / size > smallval))
    if (!overdispersion)
        warning("parameter 'size' has very large values; ",
                "replacing them by an arbitrary large value within ",
                "the parameter space. Try fitting a positive-Poisson ",
                "model instead.")
    ans
  }, list( .lmunb = lmunb, .emunb = emunb,
           .lsize = lsize, .esize = esize))),



  deriv = eval(substitute(expression({
    M1 <- 2
    NOS <- extra$NOS

    TFvec <- c(TRUE, FALSE)
    munb <- eta2theta(eta[,  TFvec, drop = FALSE], .lmunb , earg = .emunb )
    kmat <- eta2theta(eta[, !TFvec, drop = FALSE], .lsize , earg = .esize )


    smallval <- 1e-3  # Something like this is needed
    if (any(big.size <- munb / kmat < smallval)) {
        warning("parameter 'size' has very large values; ",
                "try fitting a positive-Poisson ",
                "model instead")
        kmat[big.size] <- munb[big.size] / smallval
    }


    dmunb.deta <- dtheta.deta(munb, .lmunb , earg = .emunb )
    dsize.deta <- dtheta.deta(kmat, .lsize , earg = .esize )


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




    smallno <- 1e-6
    if (FALSE && all(near.boundary <- oneminusf0 < smallno)) {
        warning("solution near the boundary; either there is no need ",
                "to fit a positive NBD or the distribution is centred ",
                "on the value 1")
        oneminusf0[near.boundary] <- smallno
        prob0[near.boundary] <- 1 - oneminusf0[near.boundary]
    }




    dl.dmunb <- y / munb - (1 + y/kmat) / (1 + munb/kmat) +
                df0.dmunb / oneminusf0
    dl.dsize <- digamma(y + kmat) - digamma(kmat) -
                (y - munb) / (munb + kmat) + log(tempk) +
                df0.dkmat / oneminusf0


    if (any(big.size)) {
      dl.dsize[big.size] <- 1e-8  # A small number
    }


    
    myderiv <- c(w) * cbind(dl.dmunb * dmunb.deta,
                            dl.dsize * dsize.deta)
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lmunb = lmunb, .lsize = lsize,
            .emunb = emunb, .esize = esize ))),


  weight = eval(substitute(expression({
    wz <- matrix(0, n, M+M-1)
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
      Q.MAXS <-      pmax(10, ceiling(1 / sqrt(eps.trig)))  #
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
          ysim <- rposnegbin(sum(ii.TF), munb = muvec, size = kkvec)
          dl.dk <- digamma(ysim + kkvec) - digamma(kkvec) -
                   (ysim - muvec) / (muvec + kkvec) +
                   log1p(-muvec / (kkvec + muvec)) +
                   df0.dkmat[ii.TF, jay] / oneminusf0[ii.TF, jay]
          run.varcov <- run.varcov + dl.dk^2
        }  # end of for loop

        run.varcov <- c(run.varcov / .nsimEIM )
        ned2l.dk2 <- if (intercept.only) mean(run.varcov) else run.varcov

        wz[ii.TF, M1*jay] <- ned2l.dk2  # * (dsize.deta[ii.TF, jay])^2
      }
    }  # jay



    wz[, M1*(1:NOS)    ] <- wz[, M1*(1:NOS)    ] * dsize.deta^2






    save.weights <- !all(ind2)


    ned2l.dmunb2 <- mymu / munb^2 -
        ((1 + mymu/kmat) / kmat) / (1 + munb/kmat)^2 -
        df02.dmunb2 / oneminusf0 -
        (df0.dmunb / oneminusf0)^2
    wz[,     M1*(1:NOS) - 1] <- ned2l.dmunb2 * dmunb.deta^2


    ned2l.dmunbsize <- (munb - mymu) / (munb + kmat)^2 -
      df02.dkmat.dmunb / oneminusf0 -
      df0.dmunb * df0.dkmat / oneminusf0^2
    wz[, M + M1*(1:NOS) - 1] <- ned2l.dmunbsize * dmunb.deta * dsize.deta




    w.wz.merge(w = w, wz = wz, n = n, M = M, ndepy = NOS)
  }), list( .cutoff.prob = cutoff.prob, .eps.trig = eps.trig,
            .max.support = max.support,
            .max.chunk.MB = max.chunk.MB,
            .nsimEIM = nsimEIM ))))

}






dposgeom <- function(x, prob, log = FALSE) {
  dgeom(x - 1, prob = prob, log = log)
}



pposgeom <- function(q, prob) {
  if (!is.Numeric(prob, positive = TRUE))
    stop("bad input for argument 'prob'")
  L <- max(length(q), length(prob))
  if (length(q)    != L) q    <- rep(q,    length.out = L)
  if (length(prob) != L) prob <- rep(prob, length.out = L)

  ifelse(q < 1, 0,
        (pgeom(q, prob) -
         dgeom(0, prob))
       / pgeom(0, prob, lower.tail = FALSE))
}



qposgeom <- function(p, prob) {




  ans <- qgeom(pgeom(0, prob, lower.tail = FALSE) * p +
               dgeom(0, prob),
               prob = prob)
  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- Inf
  ans
}



rposgeom <- function(n, prob) {
  qgeom(p = runif(n, min = dgeom(0, prob)), prob)
}









dpospois <- function(x, lambda, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  L <- max(length(x), length(lambda))
  if (length(x)      != L) x      <- rep(x,      len = L)
  if (length(lambda) != L) lambda <- rep(lambda, len = L)

  ans <- if (log.arg) {
    ifelse(x == 0, log(0.0), dpois(x, lambda, log = TRUE) -
           log1p(-exp(-lambda)))
  } else {
    ifelse(x == 0, 0, -dpois(x, lambda) / expm1(-lambda))
  }
  ans
}


ppospois <- function(q, lambda) {
  if (!is.Numeric(lambda, positive = TRUE))
    stop("bad input for argument 'lambda'")
  L <- max(length(q), length(lambda))
  if (length(q)      != L) q      <- rep(q,      length.out = L)
  if (length(lambda) != L) lambda <- rep(lambda, length.out = L)

  ifelse(q < 1, 0,
        (ppois(q, lambda) -
         dpois(0, lambda))
       / ppois(0, lambda, lower.tail = FALSE))
}


qpospois <- function(p, lambda) {


  ans <- qpois(ppois(0, lambda, lower.tail = FALSE) * p +
               dpois(0, lambda),
               lambda = lambda)

  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- Inf
  ans
}




rpospois <- function(n, lambda) {
  qpois(p = runif(n, min = dpois(0, lambda)), lambda)
}



rposnegbin <- function(n, size, prob = NULL, munb = NULL) {
  if (!is.null(munb)) {
    if (!is.null(prob))
        stop("'prob' and 'mu' both specified")
    qnbinom(p = runif(n,
                      min = dnbinom(0, size,              mu = munb)),
            size,              mu = munb)
  } else {
    qnbinom(p = runif(n,
                      min = dnbinom(0, size, prob = prob           )),
            size, prob = prob           )
  }
}




 pospoisson <- function(link = "loge",
                        type.fitted = c("mean", "lambda", "prob0"),
                        expected = TRUE,
                        ilambda = NULL, imethod = 1, zero = NULL) {

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (!is.logical(expected) || length(expected) != 1)
    stop("bad input for argument 'expected'")
  if (length( ilambda) && !is.Numeric(ilambda, positive = TRUE))
    stop("bad input for argument 'ilambda'")

  type.fitted <- match.arg(type.fitted,
                           c("mean", "lambda", "prob0"))[1]




  new("vglmff",
  blurb = c("Positive-Poisson distribution\n\n",
            "Links:    ",
            namesof("lambda", link, earg = earg, tag = FALSE)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .zero = zero ))),

  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("lambda"),
         link = .link ,
         type.fitted  = .type.fitted ,
         expected = .expected ,
         earg = .earg)
  }, list( .link = link, .earg = earg,
          .expected = expected,
          .type.fitted = type.fitted ))),

  initialize = eval(substitute(expression({
    temp5 <-
    w.y.check(w = w, y = y,
              Is.positive.y = TRUE,
              Is.integer.y = TRUE,
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
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    mynames1 <- param.names("lambda", ncoly)
    predictors.names <- namesof(mynames1, .link , earg = .earg, tag = FALSE)

    if (!length(etastart)) {
      lambda.init <- Init.mu(y = y, w = w, imethod = .imethod ,
                             imu = .ilambda )

      etastart <- theta2eta(lambda.init, .link , earg = .earg)
    }
  }), list( .link = link, .earg = earg,
            .ilambda = ilambda, .imethod = imethod,
            .type.fitted = type.fitted ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
   type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning the 'mean'.")
                     "mean"
                   }

    type.fitted <- match.arg(type.fitted,
                     c("mean", "lambda", "prob0"))[1]

    lambda <- eta2theta(eta, .link , earg = .earg )
    ans <- switch(type.fitted,
                  "mean"      = -lambda / expm1(-lambda),
                  "lambda"    = lambda,
                  "prob0"     = exp(-lambda))  # P(Y=0)
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
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <- rep( .link , len = M)
    names(misc$link) <- mynames1

    misc$earg <- vector("list", M)
    names(misc$earg) <- mynames1
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  }), list( .link = link, .earg = earg, .expected = expected ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {
    lambda <- eta2theta(eta, .link , earg = .earg ) 
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      ll.elts <- c(w) * dpospois(x = y, lambda = lambda, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("pospoisson"),


  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")
    if (any(pwts != 1)) 
      warning("ignoring prior weights")
    eta <- predict(object)
    lambda <- eta2theta(eta, .link , earg = .earg ) 
    rpospois(nsim * length(lambda), lambda)
  }, list( .link = link, .earg = earg ))),




  deriv = eval(substitute(expression({
    lambda <- eta2theta(eta, .link , earg = .earg ) 

    temp6 <- expm1(lambda)
    dl.dlambda <- y / lambda - 1 - 1 / temp6

    dlambda.deta <- dtheta.deta(lambda, .link , earg = .earg )

    c(w) * dl.dlambda * dlambda.deta
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({
    if ( .expected ) {
      ned2l.dlambda2 <- (1 + 1 / temp6) * (1/lambda - 1/temp6)
      wz <-  ned2l.dlambda2 * dlambda.deta^2
    } else {
      d2l.dlambda2 <- y / lambda^2 - (1 + 1 / temp6 + 1) / temp6
      d2lambda.deta2 <- d2theta.deta2(lambda, .link , earg = .earg)
      wz <- (dlambda.deta^2) * d2l.dlambda2 - dl.dlambda * d2lambda.deta2
    }
    c(w) * wz
  }), list( .link = link, .earg = earg, .expected = expected ))))
}








pposbinom <- function(q, size, prob 
                     ) {


  if (!is.Numeric(prob, positive = TRUE)) 
    stop("no zero or non-numeric values allowed for argument 'prob'")
  L <- max(length(q), length(size), length(prob))
  if (length(q)      != L) q      <- rep(q,      length.out = L)
  if (length(size)   != L) size   <- rep(size,   length.out = L)
  if (length(prob)   != L) prob   <- rep(prob,   length.out = L)

  ifelse(q < 1, 0,
        (pbinom(q = q, size = size, prob = prob) -
         dbinom(x = 0, size = size, prob = prob))
       / pbinom(q = 0, size = size, prob = prob, lower.tail = FALSE))
}


qposbinom <- function(p, size, prob
                     ) {




  ans <- qbinom(pbinom(0, size, prob, lower.tail = FALSE) * p +
                dbinom(0, size, prob),
                size = size, prob = prob)

  ans[p >  1] <- NaN
  ans[p <  0] <- NaN
  ans[p == 1] <- size[p == 1]
  ans
}



rposbinom <- function(n, size, prob) {
  qbinom(p = runif(n, min = dbinom(0, size, prob)), size, prob)
}



dposbinom <- function(x, size, prob, log = FALSE) {
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)


  L <- max(length(x), length(size), length(prob))
  if (length(x)      != L) x    <- rep(x,    len = L)
  if (length(size)   != L) size <- rep(size, len = L)
  if (length(prob)   != L) prob <- rep(prob, len = L)

  answer <- NaN * x
  is0 <- (x == 0)
  ok2 <- (prob > 0) & (prob <= 1) &
         (size == round(size)) & (size > 0)

  answer <-        dbinom(x = x, size = size, prob = prob, log = TRUE) -
            log1p(-dbinom(x = 0, size = size, prob = prob))
  answer[!ok2] <- NaN
  if (log.arg) {
    answer[is0 & ok2] <- log(0.0)
  } else {
    answer <- exp(answer)
    answer[is0 & ok2] <- 0.0
  }
  answer
}







 posbinomial <-
  function(link = "logit",
           multiple.responses = FALSE, parallel = FALSE,
           omit.constant = FALSE,

           p.small = 1e-4, no.warning = FALSE,

           zero = NULL) {




  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")



  if (!is.logical(multiple.responses) || length(multiple.responses) != 1)
    stop("bad input for argument 'multiple.responses'")

  if (!is.logical(omit.constant) || length(omit.constant) != 1)
    stop("bad input for argument 'omit.constant'")



  if (!is.Numeric(p.small, positive = TRUE, length.arg = 1))
    stop("bad input for argument 'p.small'")


  new("vglmff",
  blurb = c("Positive-binomial distribution\n\n",
            "Links:    ",
            if (multiple.responses)
            c(namesof("prob1", link, earg = earg, tag = FALSE),
              ",...,",
              namesof("probM", link, earg = earg, tag = FALSE)) else
            namesof("prob", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x, 
                           bool = .parallel , 
                           constraints = constraints)

    constraints <- cm.zero.VGAM(constraints, x = x, .zero , M = M,
                                predictors.names = predictors.names,
                                M1 = 1)
  }), list( .parallel = parallel, .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = 1,
         expected = TRUE,
         multipleResponses = .multiple.responses ,
         parameters.names = c("prob"),
         p.small    = .p.small ,
         no.warning = .no.warning ,
         zero = .zero )
  }, list( .zero = zero,
           .p.small    = p.small,
           .multiple.responses = multiple.responses,
           .no.warning = no.warning ))),

  initialize = eval(substitute(expression({

    mustart.orig <- mustart
    if ( .multiple.responses ) {
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
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    extra$p.small    <- .p.small
    extra$no.warning <- .no.warning

      extra$orig.w <- w
      mustart <- matrix(colSums(y) / colSums(w),  # Not colSums(y * w)...
                        n, ncoly, byrow = TRUE)

    } else {
      eval(binomialff(link = .earg ,  # earg = .earg ,
                      earg.link = TRUE)@initialize)
    }


    if ( .multiple.responses ) {

      dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
      dn2 <- if (length(dn2)) {
        paste("E[", dn2, "]", sep = "")
      } else {
        paste("prob", 1:M, sep = "")
      }
      predictors.names <-
        namesof(if (M > 1) dn2 else "prob",
                .link , earg = .earg, short = TRUE)

      w <- matrix(w, n, ncoly)
      y <- y / w  # Now sample proportion
    } else {
      predictors.names <-
        namesof("prob", .link , earg = .earg , tag = FALSE)
    }

    if (length(extra)) extra$w <- w else extra <- list(w = w)

    if (!length(etastart)) {
      mustart.use <- if (length(mustart.orig)) mustart.orig else mustart
      etastart <- cbind(theta2eta(mustart.use, .link , earg = .earg ))
    }
    mustart <- NULL



    nvec <- if (ncol(as.matrix(y)) > 1) {
              NULL
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
            }
    extra$tau <- if (length(nvec) && length(unique(nvec) == 1))
                   nvec[1] else NULL
  }), list( .link = link,
            .p.small    = p.small,
            .no.warning = no.warning,
            .earg = earg, .multiple.responses = multiple.responses ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    w <- extra$w
    binprob <- eta2theta(eta, .link , earg = .earg )
    nvec <- if ( .multiple.responses ) {
             w
           } else {
             if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
               round(w)
           }
    binprob / (1.0 - (1.0 - binprob)^nvec)
  },

  list( .link = link, .earg = earg,
        .multiple.responses = multiple.responses ))),
  last = eval(substitute(expression({
    extra$w <- NULL  # Kill it off 


    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else "prob"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg

    misc$expected <- TRUE
    misc$omit.constant <- .omit.constant
    misc$needto.omit.constant <- TRUE  # Safety mechanism
    
    
    misc$multiple.responses   <- .multiple.responses
    w <- as.numeric(w)



    if (length(extra$tau)) {
      R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
      R[lower.tri(R)] <- 0
      tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                                 R = R, w = w,
                                 X.vlm = X.vlm.save,
                                 Hlist = Hlist,  # 20150428; bug fixed here
                                 extra = extra, model.type = "0")
      extra$N.hat    <- tmp6$N.hat
      extra$SE.N.hat <- tmp6$SE.N.hat
    }

    
  }), list( .link = link, .earg = earg,
            .multiple.responses = multiple.responses,
            .omit.constant = omit.constant ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

      ycounts <- if ( .multiple.responses ) {
                  round(y * extra$orig.w)
                 } else {
                   if (is.numeric(extra$orig.w)) y * w / extra$orig.w else
                   y * w  # Convert proportions to counts
                 }
      nvec <- if ( .multiple.responses ) {
                w
              } else {
                if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                  round(w)
              }
      use.orig.w <- if (is.numeric(extra$orig.w)) extra$orig.w else 1
    binprob <- eta2theta(eta, .link , earg = .earg )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      answer <- c(use.orig.w) * dposbinom(x = ycounts, size = nvec,
                                          prob = binprob, log = TRUE)
      if ( .omit.constant ) {
        answer <- answer - c(use.orig.w) * lchoose(n = nvec, k = ycounts)
      }
      ll.elts <- answer
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg,
          .multiple.responses = multiple.responses,
          .omit.constant = omit.constant ))),

  vfamily = c("posbinomial"),




  simslot = eval(substitute(
  function(object, nsim) {

    pwts <- if (length(pwts <- object@prior.weights) > 0)
              pwts else weights(object, type = "prior")

    if ( .multiple.responses )
      stop("cannot run simulate() when 'multiple.responses = TRUE'")

    eta <- predict(object)
    binprob <- eta2theta(eta, .link , earg = .earg )

    extra <- object@extra
    w <- extra$w  # Usual code
    w <- pwts  # 20140101


    nvec <- if ( .multiple.responses ) {
              w
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
                round(w)
            }
    rposbinom(nsim * length(eta), size = nvec, prob = binprob)
  }, list( .link = link, .earg = earg,
          .multiple.responses = multiple.responses,
          .omit.constant = omit.constant ))),





  deriv = eval(substitute(expression({
    use.orig.w <- if (is.numeric(extra$orig.w)) extra$orig.w else
                  rep(1, n)

    nvec <- if ( .multiple.responses ) {
              w
            } else {
              if (is.numeric(extra$orig.w)) round(w / extra$orig.w) else
              round(w)
            }
    binprob <- eta2theta(eta, .link , earg = .earg )
    dmu.deta <- dtheta.deta(binprob, .link , earg = .earg )

    temp1 <- 1 - (1 - binprob)^nvec
    temp2 <-     (1 - binprob)^2
    temp3 <-     (1 - binprob)^(nvec-2)

    dl.dmu <- y / binprob - (1 - y) / (1 - binprob) -
             (1 - binprob) * temp3 / temp1

    c(w) * dl.dmu * dmu.deta
  }), list( .link = link, .earg = earg,
            .multiple.responses = multiple.responses ))),
  weight = eval(substitute(expression({

    ned2l.dmu2 <- 1 / (binprob * temp1) +
                  (1 - mu) / temp2 -
                  (nvec-1) * temp3 / temp1 -
                  nvec * (temp2^(nvec-1)) / temp1^2



    wz <- c(w) * ned2l.dmu2 * dmu.deta^2
    wz
  }), list( .link = link, .earg = earg,
            .multiple.responses = multiple.responses ))))
}







 posbernoulli.t <-
  function(link = "logit",

           parallel.t = FALSE ~ 1,



           iprob = NULL,

           p.small = 1e-4, no.warning = FALSE) {







  apply.parint <- FALSE



  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  if (length(iprob))
  if (!is.Numeric(iprob, positive = TRUE) ||
        max(iprob) >= 1)
    stop("argument 'iprob' must have values in (0, 1)")

  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")

  if (!is.Numeric(p.small, positive = TRUE, length.arg = 1))
    stop("bad input for argument 'p.small'")


  new("vglmff",
  blurb = c("Positive-Bernoulli (capture-recapture) model ",
            "with temporal effects (M_{t}/M_{th})\n\n",
            "Links:    ",
            namesof("prob1", link, earg = earg, tag = FALSE), ", ",
            namesof("prob2", link, earg = earg, tag = FALSE), ", ..., ",
            namesof("probM", link, earg = earg, tag = FALSE),
            "\n"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x, 
                           bool = .parallel.t , 
                           constraints = constraints,
                           apply.int = .apply.parint ,  #  TRUE,
                           cm.default = diag(M),
                           cm.intercept.default = diag(M))
  }), list( .parallel.t = parallel.t,
            .apply.parint = apply.parint ))),
  infos = eval(substitute(function(...) {
    list(M1 = 1,
         Q1 = NA,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = c("prob"),
         p.small    = .p.small ,
         no.warning = .no.warning ,
         apply.parint = .apply.parint ,
         parallel.t = .parallel.t )
  }, list( .parallel.t   = parallel.t,
           .p.small    = p.small,
           .no.warning = no.warning,          
           .apply.parint = apply.parint ))),

  initialize = eval(substitute(expression({
    M1 <- 1

    mustart.orig <- mustart
    y <- as.matrix(y)
    M <- ncoly <- ncol(y)
    extra$ncoly       <- ncoly <- ncol(y)
    extra$tau <- tau <- ncol(y)
    extra$orig.w <- w

    extra$p.small    <- .p.small
    extra$no.warning <- .no.warning
    

    w <- matrix(w, n, ncoly)
    mustart <- matrix(colSums(y) / colSums(w),
                    n, ncol(y), byrow = TRUE)
    mustart[mustart == 0] <- 0.05
    mustart[mustart == 1] <- 0.95

    if (ncoly == 1)
      stop("the response is univariate, therefore use posbinomial()")






    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")
    if (!all(w == 1))
      stop("argument 'weight' must contain 1s only")



    dn2 <- if (is.matrix(y)) dimnames(y)[[2]] else NULL
    dn2 <- if (length(dn2)) {
      paste("E[", dn2, "]", sep = "")
    } else {
      paste("prob", 1:M, sep = "")
    }


    predictors.names <- namesof(dn2, .link , earg = .earg, short = TRUE)


    if (length(extra)) extra$w <- w else extra <- list(w = w)

    if (!length(etastart)) {
      mustart.use <- if (length(mustart.orig)) {
        mustart.orig
      } else {
        mustart
      }
      etastart <- cbind(theta2eta(mustart.use, .link , earg = .earg ))
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg,
            .p.small    = p.small,
            .no.warning = no.warning
           ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    tau <- extra$ncoly
    probs <- eta2theta(eta, .link , earg = .earg )
    logAA0 <- rowSums(log1p(-probs))
    AA0 <- exp(logAA0)
    AAA <- exp(log1p(-AA0))  # 1 - AA0



    fv <- probs / AAA
    fv
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    extra$w   <- NULL   # Kill it off 


    misc$link <- rep( .link , length = M)
    names(misc$link) <- if (M > 1) dn2 else "prob"

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M) misc$earg[[ii]] <- .earg


    misc$multiple.responses  <- TRUE
    misc$iprob               <- .iprob


    R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                               R = R, w = w,
                               X.vlm = X.vlm.save,
                               Hlist = Hlist,  # 20150428; bug fixed here
                               extra = extra, model.type = "t")
    extra$N.hat    <- tmp6$N.hat
    extra$SE.N.hat <- tmp6$SE.N.hat




    misc$parallel.t   <- .parallel.t
    misc$apply.parint <- .apply.parint
  }), list( .link = link, .earg = earg,
            .parallel.t = parallel.t,
            .apply.parint = apply.parint,
            .iprob = iprob ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    ycounts <- y
    use.orig.w <- if (length(extra$orig.w)) extra$orig.w else 1

    probs <- eta2theta(eta, .link , earg = .earg )

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {


      ll.elts <-
        c(use.orig.w) *
          dposbern(x = ycounts,  # size = 1,  # Bernoulli trials
                   prob = probs, prob0 = probs, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("posbernoulli.t"),
  deriv = eval(substitute(expression({
    probs <- eta2theta(eta, .link , earg = .earg )
    dprobs.deta <- dtheta.deta(probs, .link , earg = .earg )

    logAA0 <- rowSums(log1p(-probs))
    AA0 <- exp(logAA0)
    AAA <- exp(log1p(-AA0))  # 1 - AA0

    B.s <- AA0 / (1 - probs)
    B.st <- array(AA0, c(n, M, M))
    for (slocal in 1:(M-1))
      for (tlocal in (slocal+1):M)
        B.st[, slocal, tlocal] <-
        B.st[, tlocal, slocal] <- B.s[, slocal] / (1 - probs[, tlocal])

    temp2 <-     (1 - probs)^2
    dl.dprobs <- y / probs - (1 - y) / (1 - probs) - B.s / AAA

    deriv.ans <- c(w) * dl.dprobs * dprobs.deta
    deriv.ans
  }), list( .link = link, .earg = earg ))),
  weight = eval(substitute(expression({

    ned2l.dprobs2 <- 1 / (probs * AAA) + 1 / temp2 -
                     probs / (AAA * temp2) - (B.s / AAA)^2

    wz <- matrix(NA_real_, n, dimm(M))
    wz[, 1:M] <- ned2l.dprobs2 * (dprobs.deta^2)

    for (slocal in 1:(M-1))
      for (tlocal in (slocal+1):M)
        wz[, iam(slocal, tlocal, M = M)] <- dprobs.deta[, slocal] *
                                            dprobs.deta[, tlocal] *
                                            (B.st[,slocal,tlocal] +
                                             B.s [,slocal] *
                                             B.s [,tlocal] / AAA) / (-AAA)



    wz
  }), list( .link = link, .earg = earg ))))
}






 posbernoulli.b <-
  function(link = "logit",


           drop.b = FALSE ~ 1,


           type.fitted = c("likelihood.cond", "mean.uncond"),

           I2 = FALSE,
           ipcapture = NULL,
           iprecapture = NULL,
           p.small = 1e-4, no.warning = FALSE
           ) {




  type.fitted <- match.arg(type.fitted,
                           c("likelihood.cond", "mean.uncond"))[1]

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  apply.parint.b <- FALSE


  if (length(ipcapture))
  if (!is.Numeric(ipcapture, positive = TRUE) ||
        max(ipcapture) >= 1)
    stop("argument 'ipcapture' must have values in (0, 1)")
  if (length(iprecapture))
  if (!is.Numeric(iprecapture, positive = TRUE) ||
        max(iprecapture) >= 1)
    stop("argument 'iprecapture' must have values in (0, 1)")

  if (!is.logical(I2) ||
      length(I2) != 1)
    stop("argument 'I2' must be a single logical")


  if (!is.Numeric(p.small, positive = TRUE, length.arg = 1))
    stop("bad input for argument 'p.small'")

 


  new("vglmff",
  blurb = c("Positive-Bernoulli (capture-recapture) model ",
            "with behavioural effects (M_{b}/M_{bh})\n\n",
            "Links:    ",
            namesof("pcapture",   link, earg = earg, tag = FALSE), ", ",
            namesof("precapture", link, earg = earg, tag = FALSE),
            "\n"),

  constraints = eval(substitute(expression({

    cm.intercept.default <- if ( .I2 ) diag(2) else cbind(0:1, 1)

    constraints <- cm.VGAM(matrix(1, 2, 1), x = x,
                           bool = .drop.b ,
                           constraints = constraints,
                           apply.int = .apply.parint.b ,  # TRUE, 
                           cm.default = cm.intercept.default,  # diag(2),
                           cm.intercept.default = cm.intercept.default)
  }), list( .drop.b = drop.b,
            .I2 = I2,
            .apply.parint.b = apply.parint.b ))),

  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses = FALSE,
         parameters.names = c("pcapture", "precapture"),
         p.small    = .p.small ,
         no.warning = .no.warning ,
         type.fitted = .type.fitted ,
         apply.parint.b = .apply.parint.b )
  }, list(
           .apply.parint.b = apply.parint.b,
           .p.small    = p.small,
           .no.warning = no.warning,
           .type.fitted = type.fitted
         ))),

  initialize = eval(substitute(expression({
    M1 <- 2
    if (!is.matrix(y) || ncol(y) == 1)
      stop("the response appears to be univariate")

    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")

    orig.y <- y
    extra$orig.w <- w
    extra$tau     <- tau   <- ncol(y)
    extra$ncoly   <- ncoly <- ncol(y)
    extra$type.fitted      <- .type.fitted
    extra$dimnamesy <- dimnames(y)


    extra$p.small    <- .p.small
    extra$no.warning <- .no.warning


    

    mustart.orig <- mustart
    M <- 2


    tmp3 <- aux.posbernoulli.t(y, rename = FALSE)
    y0i        <- extra$y0i  <-       tmp3$y0i
    yr0i       <- extra$yr0i <-       tmp3$yr0i
    yr1i       <- extra$yr1i <-       tmp3$yr1i
    cap1       <- extra$cap1 <-       tmp3$cap1
    cap.hist1  <- extra$cap.hist1  <- tmp3$cap.hist1


    temp5 <-
    w.y.check(w = w, y = y,
              Is.integer.y = TRUE,
              Is.nonnegative.y = TRUE,
              ncol.w.max = 1,
              ncol.y.min = 2,
              ncol.y.max = Inf,
              out.wy = TRUE,
              colsyperw = ncol(y),
              maximize = TRUE)
    w <- temp5$w  # Retain the 0-1 response
    y <- temp5$y  # Retain the 0-1 response

    mustart <- matrix(colMeans(y), n, tau, byrow = TRUE)
    mustart <- (mustart + orig.y) / 2




    predictors.names <-
      c(namesof(  "pcapture",  .link , earg = .earg, short = TRUE),
        namesof("precapture",  .link , earg = .earg, short = TRUE))

    if (!length(etastart)) {
      mustart.use <- if (length(mustart.orig)) {
        mustart.orig
      } else {
        mustart
      }

      etastart <-
        cbind(theta2eta(rowMeans(mustart.use), .link , earg = .earg ),
              theta2eta(rowMeans(mustart.use), .link , earg = .earg ))

      if (length(   .ipcapture ))
        etastart[, 1] <- theta2eta(   .ipcapture , .link , earg = .earg )
      if (length( .iprecapture ))
        etastart[, 2] <- theta2eta( .iprecapture , .link , earg = .earg )
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg,
            .type.fitted = type.fitted,
            .p.small    = p.small,
            .no.warning = no.warning,
            .ipcapture =   ipcapture,
            .iprecapture = iprecapture
          ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    cap.probs <- eta2theta(eta[, 1], .link , earg = .earg )
    rec.probs <- eta2theta(eta[, 2], .link , earg = .earg )
    tau <- extra$tau
    prc <- matrix(cap.probs, nrow(eta), tau)
    prr <- matrix(rec.probs, nrow(eta), tau)
    logQQQ <- rowSums(log1p(-prc))
    QQQ <- exp(logQQQ)
    AAA <- exp(log1p(-QQQ))  # 1 - QQQ


    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning 'likelihood.cond'.")
                     "likelihood.cond"
                   }


    type.fitted <- match.arg(type.fitted,
                             c("likelihood.cond", "mean.uncond"))[1]


    if ( type.fitted == "likelihood.cond") {
      probs.numer <- prr 
      mat.index <- cbind(1:nrow(prc), extra$cap1)
      probs.numer[mat.index] <- prc[mat.index]
      probs.numer[extra$cap.hist1 == 0] <- prc[extra$cap.hist1 == 0]
      fv <- probs.numer / AAA

    } else {


      fv <- prc - prr
      for (jay in 2:tau)
        fv[, jay] <- fv[, jay-1] * (1 - cap.probs)
      fv <- (fv + prr) / AAA
    }



    ans <- fv
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
  }, list( .link = link,
           .type.fitted = type.fitted,
           .earg = earg ))),
  last = eval(substitute(expression({

    misc$link <- c( .link , .link )
    names(misc$link) <- predictors.names

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    misc$earg[[1]] <- .earg
    misc$earg[[2]] <- .earg

    misc$expected           <- TRUE
    misc$multiple.responses <- TRUE
    misc$ipcapture   <- .ipcapture
    misc$iprecapture <- .iprecapture
    misc$drop.b      <- .drop.b
    misc$multipleResponses <- FALSE
    misc$apply.parint.b <- .apply.parint.b



    R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                               R = R, w = w,
                               X.vlm = X.vlm.save,
                               Hlist = Hlist,  # 20150428; bug fixed here
                               extra = extra, model.type = "b")
    extra$N.hat    <- tmp6$N.hat
    extra$SE.N.hat <- tmp6$SE.N.hat


  }), list( .link = link, .earg = earg,
            .drop.b = drop.b,
            .ipcapture =   ipcapture,
            .iprecapture = iprecapture,
            .apply.parint.b = apply.parint.b
          ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    tau <- extra$ncoly
    ycounts <- y
    use.orig.w <- if (length(extra$orig.w)) extra$orig.w else 1

    cap.probs <- eta2theta(eta[, 1], .link , earg = .earg )
    rec.probs <- eta2theta(eta[, 2], .link , earg = .earg )
    prc <- matrix(cap.probs, nrow(eta), tau)
    prr <- matrix(rec.probs, nrow(eta), tau)

    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      probs.numer <- prr 
      mat.index <- cbind(1:nrow(prc), extra$cap1)
      probs.numer[mat.index] <- prc[mat.index]
      probs.numer[extra$cap.hist1 == 0] <- prc[extra$cap.hist1 == 0]

      ll.elts <-
        c(use.orig.w) *
          dposbern(x = ycounts,  # Bernoulli trials
                   prob = probs.numer, prob0 = prc, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("posbernoulli.b"),












  deriv = eval(substitute(expression({
    cap.probs <- eta2theta(eta[, 1], .link , earg = .earg )
    rec.probs <- eta2theta(eta[, 2], .link , earg = .earg )
    y0i  <- extra$y0i
    yr0i <- extra$yr0i
    yr1i <- extra$yr1i
    cap1 <- extra$cap1
    tau  <- extra$tau

    dcapprobs.deta <- dtheta.deta(cap.probs, .link , earg = .earg )
    drecprobs.deta <- dtheta.deta(rec.probs, .link , earg = .earg )

    QQQ <- (1 - cap.probs)^tau
    dl.dcap <-   1  /      cap.probs -
               y0i  / (1 - cap.probs) -
               tau * ((1 - cap.probs)^(tau - 1)) / (1 - QQQ)

    dl.drec <- yr1i /      rec.probs -
               yr0i / (1 - rec.probs)


    deriv.ans <- c(w) * cbind(dl.dcap * dcapprobs.deta,
                              dl.drec * drecprobs.deta)
    deriv.ans
  }), list( .link = link, .earg = earg ))),

  weight = eval(substitute(expression({

    wz <- matrix(0, n, M)  # Diagonal EIM



    dA.dcapprobs <- -tau * ((1 - QQQ) * (tau-1) * (1 - cap.probs)^(tau-2) +
                     tau * (1 - cap.probs)^(2*tau -2)) / (1 - QQQ)^2





    prc <- matrix(cap.probs, n, tau)
    prr <- matrix(rec.probs, n, tau)

    dQ.dprc   <- -QQQ / (1 - prc)
    QQQcummat <- exp(t( apply(log1p(-prc), 1, cumsum)))



    GGG <- (1 - QQQ - cap.probs * (1 + (tau-1) * QQQ)) / (
            cap.probs * (1-cap.probs)^2)
    wz.pc <- GGG / (1 - QQQ) + 1 / cap.probs^2 + dA.dcapprobs
    wz[, iam(1, 1, M = M)] <- wz.pc * dcapprobs.deta^2  # Efficient





    wz.pr <- (tau - (1 - QQQ) / cap.probs) / (
              rec.probs * (1 - rec.probs) * (1 - QQQ))
    wz[, iam(2, 2, M = M)] <- wz.pr * drecprobs.deta^2

  


    wz <- c(w) * wz
    wz
  }), list( .link = link, .earg = earg ))))
}





 posbernoulli.tb <-
  function(link = "logit",
           parallel.t = FALSE ~  1,
           parallel.b = FALSE ~  0,
           drop.b     = FALSE ~  1,
           type.fitted = c("likelihood.cond", "mean.uncond"),
           imethod = 1,
           iprob = NULL,
           p.small = 1e-4, no.warning = FALSE,  
           ridge.constant = 0.01,
           ridge.power = -4) {



  apply.parint.t <- FALSE
  apply.parint.b <- TRUE
  apply.parint.d <- FALSE  # For 'drop.b' actually.

  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")

  type.fitted <- match.arg(type.fitted,
                           c("likelihood.cond", "mean.uncond"))[1]


  if (!is.Numeric(imethod, length.arg = 1,
                  integer.valued = TRUE, positive = TRUE) ||
      imethod > 2)
    stop("argument 'imethod' must be 1 or 2")


  if (!is.Numeric(ridge.constant) ||
      ridge.constant < 0)
    warning("argument 'ridge.constant' should be non-negative")
  if (!is.Numeric(ridge.power) ||
      ridge.power > 0)
    warning("argument 'ridge.power' should be non-positive")


  if (length(iprob))
    if (!is.Numeric(iprob, positive = TRUE) ||
          max(iprob) >= 1)
      stop("argument 'iprob' must have values in (0, 1)")


  if (!is.Numeric(p.small, positive = TRUE, length.arg = 1))
    stop("bad input for argument 'p.small'")


  
  new("vglmff",
  blurb = c("Positive-Bernoulli (capture-recapture) model\n",
            "with temporal and behavioural effects (M_{tb}/M_{tbh})\n\n",
            "Links:    ",
            namesof("pcapture.1",     link, earg = earg, tag = FALSE),
            ", ..., ",
            namesof("pcapture.tau",   link, earg = earg, tag = FALSE), ", ",
            namesof("precapture.2",   link, earg = earg, tag = FALSE),
            ", ..., ",
            namesof("precapture.tau", link, earg = earg, tag = FALSE)),
  constraints = eval(substitute(expression({
 

    constraints.orig <- constraints
    cm1.d <-
    cmk.d <- matrix(0, M, 1)  # All 0s inside
    con.d <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .drop.b ,
                           constraints = constraints.orig,
                           apply.int = .apply.parint.d ,  # FALSE,  
                           cm.default           = cmk.d,
                           cm.intercept.default = cm1.d)
   


    cm1.t <-
    cmk.t <- rbind(diag(tau), diag(tau)[-1, ])  # More readable
    con.t <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel.t ,  # Same as .parallel.b
                           constraints = constraints.orig,
                           apply.int = .apply.parint.t ,  # FALSE,  
                           cm.default           = cmk.t,
                           cm.intercept.default = cm1.t)
   
    

    cm1.b <-
    cmk.b <- rbind(matrix(0, tau, tau-1), diag(tau-1))
    con.b <- cm.VGAM(matrix(c(rep(0, len = tau  ),
                              rep(1, len = tau-1)), M, 1), x = x,
                           bool = .parallel.b ,  # Same as .parallel.b
                           constraints = constraints.orig,
                           apply.int = .apply.parint.b ,  # FALSE,  
                           cm.default           = cmk.b,
                           cm.intercept.default = cm1.b)
   
    con.use <- con.b
    for (klocal in 1:length(con.b)) {
      con.use[[klocal]] <-
        cbind(if (any(con.d[[klocal]] == 1)) NULL else con.b[[klocal]],
              con.t[[klocal]])

    }

    
    constraints <- con.use
    
  }), list( .parallel.t = parallel.t,
            .parallel.b = parallel.b,
            .drop.b     = drop.b,
            .apply.parint.b = apply.parint.b,
            .type.fitted    = type.fitted,
            .apply.parint.d = apply.parint.d,
            .apply.parint.t = apply.parint.t ))),
  infos = eval(substitute(function(...) {
    list(M1 = 2,
         expected = TRUE,
         multipleResponses  = TRUE,
         parameters.names = as.character(NA),
         ridge.constant     = .ridge.constant ,
         ridge.power        = .ridge.power ,
         drop.b             = .drop.b,
         imethod            = .imethod ,
         type.fitted        = .type.fitted ,
         p.small    = .p.small ,
         no.warning = .no.warning ,
         apply.parint.b     = .apply.parint.b ,
         apply.parint.t     = .apply.parint.t ,
         parallel.t         = .parallel.t ,
         parallel.b         = .parallel.b )
  }, list( .parallel.t         = parallel.t,
           .parallel.b         = parallel.b,
           .drop.b             = drop.b,
           .type.fitted        = type.fitted,
           .p.small    = p.small,
           .no.warning = no.warning,
           .imethod            = imethod,
           .ridge.constant     = ridge.constant,
           .ridge.power        = ridge.power,
           .apply.parint.b     = apply.parint.b,
           .apply.parint.t     = apply.parint.t ))),

  initialize = eval(substitute(expression({
    M1 <- 2  # Not quite true


    if (ncol(cbind(w)) > 1)
      stop("variable 'w' should be a vector or one-column matrix")
    w <- c(w)  # Make it a vector

    mustart.orig <- mustart
    y <- as.matrix(y)
    extra$tau     <- tau   <- ncol(y)
    extra$ncoly   <- ncoly <- ncol(y)
    extra$orig.w  <- w
    extra$ycounts <- y
    extra$type.fitted <- .type.fitted
    extra$dimnamesy <- dimnames(y)
    M <- M1 * tau - 1  # recap.prob.1 is unused


    mustart <- (y + matrix(apply(y, 2, weighted.mean, w = w),
                           n, tau, byrow = TRUE)) / 2
    mustart[mustart < 0.01] <- 0.01
    mustart[mustart > 0.99] <- 0.99

    mustart <- cbind(mustart, mustart[, -1])



 
    extra$p.small    <- .p.small
    extra$no.warning <- .no.warning

   



    if (!all(y == 0 | y == 1))
      stop("response must contain 0s and 1s only")


    tmp3 <- aux.posbernoulli.t(y)
    cap.hist1  <- extra$cap.hist1  <- tmp3$cap.hist1
    

    dn2.cap   <- paste("pcapture.",   1:ncoly, sep = "")
    dn2.recap <- paste("precapture.", 2:ncoly, sep = "")

    predictors.names <- c(
      namesof(dn2.cap,   .link , earg = .earg, short = TRUE),
      namesof(dn2.recap, .link , earg = .earg, short = TRUE))


    if (length(extra)) extra$w <- w else extra <- list(w = w)

    if (!length(etastart)) {
      mu.init <-
        if ( .imethod == 1) {
          if (length( .iprob ))
            matrix( .iprob , n, M, byrow = TRUE) else
          if (length(mustart.orig))
            matrix(rep(mustart.orig, length = n * M), n, M) else
            mustart  # Already n x M
        } else {
          matrix(runif(n * M), n, M)
        }
      etastart <- theta2eta(mu.init, .link , earg = .earg )  # n x M
    }
    mustart <- NULL
  }), list( .link = link, .earg = earg,
            .type.fitted = type.fitted,
            .p.small    = p.small,
            .no.warning = no.warning,
            .iprob = iprob,
            .imethod = imethod ))),

  linkinv = eval(substitute(function(eta, extra = NULL) {
    tau <- extra$ncoly
    taup1 <- tau + 1
    probs <- eta2theta(eta, .link , earg = .earg )
    prc <- probs[, 1:tau]
    prr <- cbind(0,  # == pr1.ignored
                 probs[, taup1:ncol(probs)])  # 1st coln ignored

    logQQQ <- rowSums(log1p(-prc))
    QQQ <- exp(logQQQ)
    AAA <- exp(log1p(-QQQ))  # 1 - QQQ

    type.fitted <- if (length(extra$type.fitted)) extra$type.fitted else {
                     warning("cannot find 'type.fitted'. ",
                             "Returning 'likelihood.cond'.")
                     "likelihood.cond"
                   }

  type.fitted <- match.arg(type.fitted,
                           c("likelihood.cond", "mean.uncond"))[1]



    if ( type.fitted == "likelihood.cond") {
      probs.numer <- prr 
      mat.index <- cbind(1:nrow(prc), extra$cap1)
      probs.numer[mat.index] <- prc[mat.index]
      probs.numer[extra$cap.hist1 == 0] <- prc[extra$cap.hist1 == 0]
      fv <- probs.numer / AAA
    } else {
      fv <- matrix(prc[, 1] / AAA, nrow(prc), ncol(prc))

      fv[, 2] <- (prc[, 2] + prc[, 1] * (prr[, 2] - prc[, 2])) / AAA

      if (tau >= 3) {
        QQQcummat <- exp(t( apply(log1p(-prc), 1, cumsum)))
        for (jay in 3:tau) {
          sum1 <- prc[, 1]
          for (kay in 2:(jay-1))
            sum1 <- sum1 + prc[, kay] * QQQcummat[, kay-1]
          fv[, jay] <- prc[, jay] * QQQcummat[, jay-1] +
                       prr[, jay] * sum1
        }
        fv[, 3:tau] <- fv[, 3:tau] / AAA
      }
    }



    ans <- fv
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
  }, list( .link = link,
           .earg = earg ))),
  last = eval(substitute(expression({
    extra$w   <- NULL   # Kill it off 


    misc$link <- rep( .link , length = M)
    names(misc$link) <- c(dn2.cap, dn2.recap)

    misc$earg <- vector("list", M)
    names(misc$earg) <- names(misc$link)
    for (ii in 1:M)
      misc$earg[[ii]] <- .earg


    misc$multiple.responses <- TRUE
    misc$iprob              <- .iprob



    R <- tfit$qr$qr[1:ncol.X.vlm, 1:ncol.X.vlm, drop = FALSE]
    R[lower.tri(R)] <- 0
    tmp6 <- N.hat.posbernoulli(eta = eta, link = .link , earg = .earg ,
                               R = R, w = w,
                               X.vlm = X.vlm.save,
                               Hlist = Hlist,  # 20150428; bug fixed here
                               extra = extra, model.type = "tb")
    extra$N.hat    <- tmp6$N.hat
    extra$SE.N.hat <- tmp6$SE.N.hat


    misc$drop.b             <- .drop.b
    misc$parallel.t         <- .parallel.t
    misc$parallel.b         <- .parallel.b
    misc$apply.parint.b     <- .apply.parint.b
    misc$apply.parint.t     <- .apply.parint.t
    misc$ridge.constant <- .ridge.constant
    misc$ridge.power    <- .ridge.power

  }), list( .link = link, .earg = earg,
            .apply.parint.b = apply.parint.b,
            .apply.parint.t = apply.parint.t,
            .parallel.t = parallel.t,
            .parallel.b = parallel.b,
            .drop.b     = drop.b,
            .ridge.constant = ridge.constant,
            .ridge.power = ridge.power,
            .iprob = iprob ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta,
             extra = NULL,
             summation = TRUE) {

    tau <- extra$ncoly
    taup1 <- tau + 1
    ycounts <- y
    use.orig.w <- if (length(extra$orig.w)) extra$orig.w else 1

    probs <- eta2theta(eta, .link , earg = .earg )
    prc <- probs[, 1:tau]

    prr <- cbind(0,  # pr1.ignored
                 probs[, taup1:ncol(probs)])  # 1st coln ignored


    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    } else {
      probs.numer <- prr 
      mat.index <- cbind(1:nrow(prc), extra$cap1)
      probs.numer[mat.index] <- prc[mat.index]
      probs.numer[extra$cap.hist1 == 0] <- prc[extra$cap.hist1 == 0]

      ll.elts <-
        c(use.orig.w) *
            dposbern(x = ycounts,  # size = 1,  # Bernoulli trials
                     prob = probs.numer, prob0 = prc, log = TRUE)
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),
  vfamily = c("posbernoulli.tb"),
  deriv = eval(substitute(expression({
    tau <- extra$ncoly
    taup1 <- tau + 1
    probs <- eta2theta(eta, .link , earg = .earg )

    prc <- probs[, 1:tau]
    prr <- cbind(pr1.ignored = 0,
                 probs[, taup1:ncol(probs)])  # 1st coln ignored
    logQQQ <- rowSums(log1p(-prc))
    QQQ <- exp(logQQQ)
 
    
    
    dprobs.deta <- dtheta.deta(probs, .link , earg = .earg )
    dQ.dprc   <- -QQQ / (1 - prc)
    d2Q.dprc <- array(0, c(n, tau, tau))
    for (jay in 1:(tau-1))
      for (kay in (jay+1):tau)
        d2Q.dprc[, jay, kay] <-
        d2Q.dprc[, kay, jay] <-  QQQ / ((1 - prc[, jay]) *
                                        (1 - prc[, kay]))

    dl.dpc <- dl.dpr <- matrix(0, n, tau)  # First coln of dl.dpr is ignored
    for (jay in 1:tau) {
      dl.dpc[, jay] <- (1 - extra$cap.hist1[, jay]) *
        (    y[, jay]  /      prc[, jay]   -
        (1 - y[, jay]) / (1 - prc[, jay])) +
            dQ.dprc[, jay] / (1 - QQQ)
    }
    for (jay in 2:tau) {
      dl.dpr[, jay] <- extra$cap.hist1[, jay] *
        (    y[, jay]  /      prr[, jay] -
        (1 - y[, jay]) / (1 - prr[, jay]))
    }

    deriv.ans <- c(w) * cbind(dl.dpc, dl.dpr[, -1]) * dprobs.deta
    deriv.ans
  }), list( .link = link,
            .earg = earg ))),

  weight = eval(substitute(expression({
    wz <- matrix(0, n, sum(M:(M - (tau - 1))))



    QQQcummat <- exp(t( apply(log1p(-prc), 1, cumsum)))
    wz.pc <- (QQQcummat / prc - QQQ / (1 - QQQ)) / ((1 - QQQ) *
              (1 - prc)^2)
    wz[, 1:tau] <- wz.pc
  

    wz.pr <- as.matrix((1 - QQQcummat / (1 - prc)) / (
                        prr * (1 - prr) * (1 - QQQ)))
    wz[, taup1:M] <- wz.pr[, -1]
  

    for (jay in 1:(tau-1))
      for (kay in (jay+1):tau)
        wz[, iam(jay, kay, M = M)] <-
          -(d2Q.dprc[, jay, kay] +
             dQ.dprc[, jay] *
             dQ.dprc[, kay] / (1 - QQQ)) / (1 - QQQ)


    cindex <- iam(NA, NA, M = M, both = TRUE)
    cindex$row.index <- cindex$row.index[1:ncol(wz)]
    cindex$col.index <- cindex$col.index[1:ncol(wz)]

    wz <- wz * dprobs.deta[, cindex$row.index] *
               dprobs.deta[, cindex$col.index]


      wz.mean <- mean(wz[, 1:tau])
      wz.adjustment <- wz.mean * .ridge.constant * iter^( .ridge.power )
      wz[, 1:tau] <- wz[, 1:tau] + wz.adjustment

    c(w) * wz
  }), list( .link = link, .earg = earg,
            .ridge.constant = ridge.constant,
            .ridge.power = ridge.power ))))
}







setClass("posbernoulli.tb",     contains = "vglmff")
setClass("posbernoulli.t",      contains = "posbernoulli.tb")
setClass("posbernoulli.b",      contains = "posbernoulli.tb")

 setClass("posbinomial",        contains = "posbernoulli.b")



setMethod("summaryvglmS4VGAM",  signature(VGAMff = "posbernoulli.tb"),
  function(object,
           VGAMff,
           ...) {
  object@post
})



setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "posbernoulli.tb"),
  function(object,
           VGAMff,
           ...) {
 if (length(object@extra$N.hat) == 1 &&
      is.numeric(object@extra$N.hat)) {
    cat("\nEstimate of N: ", round(object@extra$N.hat, digits = 3), "\n")
    cat("\nStd. Error of N: ", round(object@extra$SE.N.hat, digits = 3), "\n")

    confint.N <- object@extra$N.hat + c(Lower = -1, Upper = 1) *
                                      qnorm(0.975) * object@extra$SE.N.hat
    cat("\nApproximate 95 percent confidence interval for N:\n")
    print(round(confint.N, digits = 2))
  }
})



setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "posbernoulli.b"),
  function(object,
           VGAMff,
           ...) {
  callNextMethod(VGAMff = VGAMff, object = object, ...)
})



setMethod("showsummaryvglmS4VGAM",  signature(VGAMff = "posbernoulli.t"),
  function(object,
           VGAMff,
           ...) {
  callNextMethod(VGAMff = VGAMff, object = object, ...)
})






