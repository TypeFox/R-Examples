# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.








rrar.Ci <- function(i, coeffs, aa, Ranks., MM) {
  index <- cumsum(c(aa, MM*Ranks.))
  ans <- matrix(coeffs[(index[i]+1):index[i+1]],
                Ranks.[i], MM, byrow = TRUE)
  t(ans)
}


rrar.Ak1 <- function(MM, coeffs, Ranks., aa) {
  ptr <- 0
  Ak1 <- diag(MM)
  for (jay in 1:MM) {
    for (i in 1:MM) {
      if (i > jay && (MM+1)-(Ranks.[jay]-1) <= i) {
        ptr <- ptr + 1
        Ak1[i,jay] <- coeffs[ptr]
      }
    }
  }
  if (aa > 0 && ptr != aa) stop("something wrong")
  Ak1
}


rrar.Di <- function(i, Ranks.) {
  if (Ranks.[1] == Ranks.[i]) diag(Ranks.[i]) else 
  rbind(diag(Ranks.[i]),
        matrix(0, Ranks.[1] - Ranks.[i], Ranks.[i]))
}


rrar.Mi <- function(i, MM, Ranks., ki) {
  if (Ranks.[ki[i]] == MM)
    return(NULL)
  hi <- Ranks.[ki[i]] - Ranks.[ki[i+1]]
  Ji <- matrix(0, hi, Ranks.[1])
  for (j in 1:hi) {
    Ji[j,j+Ranks.[ki[i+1]]] <- 1
  }
  Mi <- matrix(0, MM-Ranks.[ki[i]], MM)  # dim(Oi) == dim(Ji)
  for (j in 1:(MM-Ranks.[ki[i]])) {
    Mi[j,j+Ranks.[ki[i  ]]] <- 1
  }
  kronecker(Mi, Ji)
}


rrar.Mmat <- function(MM, uu, Ranks., ki) {
  Mmat <- NULL
  for (ii in uu:1) {
    Mmat <- rbind(Mmat, rrar.Mi(ii, MM, Ranks., ki))
  }
  Mmat
}


block.diag <- function(A, B) {
  if (is.null(A) && is.null(B))
    return(NULL)
  if (!is.null(A) && is.null(B))
    return(A)
  if (is.null(A) && !is.null(B))
    return(B)

  A <- as.matrix(A)
  B <- as.matrix(B)
  temp <-  cbind(A, matrix(0, nrow(A), ncol(B)))
  rbind(temp, cbind(matrix(0, nrow(B), ncol(A)), B))
}


rrar.Ht <- function(plag, MM, Ranks., coeffs, aa, uu, ki) {
  Htop <- Hbot <- NULL
  Mmat <- rrar.Mmat(MM, uu, Ranks., ki)   # NULL if full rank
  Ak1 <- rrar.Ak1(MM, coeffs, Ranks., aa)

  if (!is.null(Mmat))
  for (i in 1:plag) {
    Di <- rrar.Di(i, Ranks.)
    Ci <- rrar.Ci(i, coeffs, aa, Ranks., MM)
    temp <- Di %*% t(Ci)
    Htop <- cbind(Htop, Mmat %*% kronecker(diag(MM), temp))
  }

  for (i in 1:plag) {
    Di <- rrar.Di(i, Ranks.)
    temp <- kronecker(t(Di) %*% t(Ak1), diag(MM))
    Hbot <- block.diag(Hbot, temp)
  }
  rbind(Htop, Hbot)
}


rrar.Ut <- function(y, tt, plag, MM) {
  Ut <- NULL
  if (plag>1)
  for (i in 1:plag) {
    Ut <- rbind(Ut, kronecker(diag(MM), cbind(y[tt-i,])))
  }
  Ut
}


rrar.UU <- function(y, plag, MM, n) {
  UU <- NULL
  for (i in (plag+1):n) {
    UU <- rbind(UU, t(rrar.Ut(y, i, plag, MM)))
  }
  UU
}


rrar.Wmat <- function(y, Ranks., MM, ki, plag, aa, uu, n, coeffs) {
  temp1 <- rrar.UU(y, plag, MM, n)
  temp2 <- t(rrar.Ht(plag, MM, Ranks., coeffs, aa, uu, ki))
  list(UU = temp1, Ht = temp2)
}



rrar.control <- function(stepsize = 0.5, save.weights = TRUE, ...) {

  if (stepsize <= 0 || stepsize > 1) {
    warning("bad value of stepsize; using 0.5 instead")
    stepsize <- 0.5
  }
  list(stepsize = stepsize,
       save.weights = as.logical(save.weights)[1])
}



 rrar <- function(Ranks = 1, coefstart = NULL) {
  lag.p <- length(Ranks)

  new("vglmff",
  blurb = c("Nested reduced-rank vector autoregressive model AR(",
            lag.p, ")\n\n",
            "Link:     ",
            namesof("mu_t", "identitylink"),
            ", t = ", paste(paste(1:lag.p, coll = ",", sep = ""))),
  initialize = eval(substitute(expression({
      Ranks. <- .Ranks
      plag <- length(Ranks.)
      nn <- nrow(x)  # original n
      indices <- 1:plag

      copy.X.vlm <- TRUE  # X.vlm.save matrix changes at each iteration 

      dsrank <- -sort(-Ranks.)  # ==rev(sort(Ranks.))
      if (any(dsrank != Ranks.))
          stop("Ranks must be a non-increasing sequence")
      if (!is.matrix(y) || ncol(y) == 1) {
          stop("response must be a matrix with more than one column")
      } else {
          MM <- ncol(y)
          ki <- udsrank <- unique(dsrank)
          uu <- length(udsrank)
          for (i in 1:uu)
             ki[i] <- max((1:plag)[dsrank == udsrank[i]])
          ki <- c(ki, plag+1)  # For computing a
          Ranks. <- c(Ranks., 0)  # For computing a
          aa <- sum( (MM-Ranks.[ki[1:uu]]) * (Ranks.[ki[1:uu]]-Ranks.[ki[-1]]) )
      }
      if (!intercept.only)
        warning("ignoring explanatory variables")

      if (any(MM < Ranks.))
        stop("'max(Ranks)' can only be ", MM, " or less")
      y.save <- y  # Save the original
      if (any(w != 1))
        stop("all weights should be 1")

      new.coeffs <- .coefstart  # Needed for iter = 1 of $weight
      new.coeffs <- if (length(new.coeffs))
                        rep(new.coeffs, len = aa+sum(Ranks.)*MM) else
                        runif(aa+sum(Ranks.)*MM)
      temp8 <- rrar.Wmat(y.save, Ranks., MM, ki, plag,
                         aa, uu, nn, new.coeffs)
      X.vlm.save <- temp8$UU %*% temp8$Ht 

      if (!length(etastart)) {
        etastart <- X.vlm.save %*% new.coeffs
        etastart <- matrix(etastart, ncol = ncol(y), byrow = TRUE)
      }

      extra$Ranks. <- Ranks.; extra$aa <- aa
      extra$plag <- plag; extra$nn <- nn
      extra$MM <- MM; extra$coeffs <- new.coeffs;
      extra$y.save <- y.save

      keep.assign <- attr(x, "assign")
      x <- x[-indices, , drop = FALSE]
      if (is.R())
          attr(x, "assign") <- keep.assign
      y <- y[-indices, , drop = FALSE]
      w <- w[-indices]
      n.save <- n <- nn - plag
  }), list( .Ranks = Ranks, .coefstart = coefstart ))), 

  linkinv = function(eta, extra = NULL) {
    aa <- extra$aa
    coeffs <- extra$coeffs
    MM <- extra$MM
    nn <- extra$nn
    plag <- extra$plag
    Ranks. <- extra$Ranks.
    y.save <- extra$y.save

    tt <- (1+plag):nn
    mu <- matrix(0, nn-plag, MM)
    Ak1 <- rrar.Ak1(MM, coeffs, Ranks., aa)
    for (i in 1:plag) {
      Di <- rrar.Di(i, Ranks.)
      Ci <- rrar.Ci(i, coeffs, aa, Ranks., MM)
      mu <- mu + y.save[tt-i, , drop = FALSE] %*%
                 t(Ak1 %*% Di %*% t(Ci))
    }
    mu
  },
  last = expression({
    misc$plag <- plag
    misc$Ranks <- Ranks.
    misc$Ak1 <- Ak1
    misc$omegahat <- omegahat
    misc$Cmatrices <- Cmatrices
    misc$Dmatrices <- Dmatrices
    misc$Hmatrix <- temp8$Ht
    misc$Phimatrices <- vector("list", plag)
    for (ii in 1:plag) {
      misc$Phimatrices[[ii]] <- Ak1 %*% Dmatrices[[ii]] %*%
                                t(Cmatrices[[ii]])
    }
    misc$Z <- y.save %*% t(solve(Ak1)) 
  }),
  vfamily = "rrar",
  deriv = expression({
    temp8 <- rrar.Wmat(y.save, Ranks., MM, ki, plag,
                       aa, uu, nn, new.coeffs)
    X.vlm.save <- temp8$UU %*% temp8$Ht 

    extra$coeffs <- new.coeffs

    resmat <- y
    tt <- (1+plag):nn
    Ak1 <- rrar.Ak1(MM, new.coeffs, Ranks., aa)
    Cmatrices <- Dmatrices <- vector("list", plag)
    for (ii in 1:plag) {
      Dmatrices[[ii]] <- Di <- rrar.Di(ii, Ranks.)
      Cmatrices[[ii]] <- Ci <- rrar.Ci(ii, new.coeffs, aa, Ranks., MM)
      resmat <- resmat - y.save[tt - ii, , drop = FALSE] %*%
                         t(Ak1 %*% Di %*% t(Ci))
    }
    omegahat <- (t(resmat) %*% resmat) / n  # MM x MM
    omegainv <- solve(omegahat)

    omegainv <- solve(omegahat)
    ind1 <- iam(NA, NA, MM, both = TRUE)

    wz <- matrix(omegainv[cbind(ind1$row, ind1$col)],
                 nn-plag, length(ind1$row), byrow = TRUE)
    mux22(t(wz), y-mu, M = extra$MM, as.matrix = TRUE)
  }),
  weight = expression({
    wz
  }))
}








vglm.garma.control <- function(save.weights = TRUE, ...) {
    list(save.weights = as.logical(save.weights)[1])
}


 garma <- function(link = "identitylink",
                   p.ar.lag = 1,
                   q.ma.lag = 0,
                   coefstart = NULL,
                   step = 1.0) {

  if (!is.Numeric(p.ar.lag, integer.valued = TRUE, length.arg = 1))
    stop("bad input for argument 'p.ar.lag'")
  if (!is.Numeric(q.ma.lag, integer.valued = TRUE, length.arg = 1))
    stop("bad input for argument 'q.ma.lag'")
  if (q.ma.lag != 0)
    stop("sorry, only q.ma.lag = 0 is currently implemented")


  link <- as.list(substitute(link))
  earg <- link2list(link)
  link <- attr(earg, "function.name")


  new("vglmff",
  blurb = c("GARMA(", p.ar.lag, ",", q.ma.lag, ")\n\n",
            "Link:     ",
            namesof("mu_t", link, earg = earg),
            ", t = ", paste(paste(1:p.ar.lag, coll = ",", sep = ""))),
  initialize = eval(substitute(expression({
    plag <- .p.ar.lag
    predictors.names <- namesof("mu", .link , earg = .earg , tag = FALSE)

    indices <- 1:plag
    tt.index <- (1 + plag):nrow(x)
    p.lm <- ncol(x)

    copy.X.vlm <- TRUE  # x matrix changes at each iteration 

    if ( .link == "logit"   || .link == "probit" ||
         .link == "cloglog" || .link == "cauchit") {
        delete.zero.colns <- TRUE
        eval(process.categorical.data.VGAM)
        mustart <- mustart[tt.index, 2]
        y <- y[, 2]
    } else {
    }


    x.save <- x  # Save the original
    y.save <- y  # Save the original
    w.save <- w  # Save the original

    new.coeffs <- .coefstart  # Needed for iter = 1 of @weight
    new.coeffs <- if (length(new.coeffs))
                    rep(new.coeffs, len = p.lm + plag) else
                    c(rnorm(p.lm, sd = 0.1), rep(0, plag)) 

    if (!length(etastart)) {
      etastart <- x[-indices, , drop = FALSE] %*% new.coeffs[1:p.lm]
    }

    x <- cbind(x, matrix(NA_real_, n, plag))  # Right size now
    dx <- dimnames(x.save)
    morenames <- paste("(lag", 1:plag, ")", sep = "") 
    dimnames(x) <- list(dx[[1]], c(dx[[2]], morenames)) 

    x <- x[-indices, , drop = FALSE]
    class(x) <- "matrix"
    y <- y[-indices]
    w <- w[-indices]
    n.save <- n <- n - plag

    more <- vector("list", plag)
    names(more) <- morenames
    for (ii in 1:plag)
      more[[ii]] <- ii + max(unlist(attr(x.save, "assign")))
    attr(x, "assign") <- c(attr(x.save, "assign"), more)
  }), list( .link = link, .p.ar.lag = p.ar.lag,
            .coefstart = coefstart, .earg = earg ))), 
  linkinv = eval(substitute(function(eta, extra = NULL) {
    eta2theta(eta, link = .link , earg = .earg)
  }, list( .link = link, .earg = earg ))),
  last = eval(substitute(expression({
    misc$link <-    c(mu = .link )

    misc$earg <- list(mu = .earg )

    misc$plag <- plag
  }), list( .link = link, .earg = earg ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL,
             summation = TRUE) {
    if (residuals) {
      switch( .link ,
              identitylink   = y - mu,
              loge       = w * (y / mu - 1),
              reciprocal = w * (y / mu - 1),
              inverse    = w * (y / mu - 1),
              w * (y / mu - (1-y) / (1 - mu)))
    } else {
      ll.elts <-
      switch( .link ,
              identitylink   = c(w) * (y - mu)^2,
              loge       = c(w) * (-mu + y * log(mu)),
              reciprocal = c(w) * (-mu + y * log(mu)),
              inverse    = c(w) * (-mu + y * log(mu)),
                           c(w) * (y * log(mu) + (1-y) * log1p(-mu)))
      if (summation) {
        sum(ll.elts)
      } else {
        ll.elts
      }
    }
  }, list( .link = link, .earg = earg ))),

  middle2 = eval(substitute(expression({
    realfv <- fv
    for (ii in 1:plag) {
      realfv <- realfv + old.coeffs[ii + p.lm] *
        (x.save[tt.index-ii, 1:p.lm, drop = FALSE] %*%
         new.coeffs[1:p.lm])  # +
    }

    true.eta <- realfv + offset  
    mu <- family@linkinv(true.eta, extra)  # overwrite mu with correct one
  }), list( .link = link, .earg = earg ))),
  vfamily = c("garma", "vglmgam"),
  deriv = eval(substitute(expression({
    dl.dmu <- switch( .link ,
                  identitylink = y-mu,
                  loge       = (y - mu) / mu,
                  reciprocal = (y - mu) / mu,
                  inverse    = (y - mu) / mu,
                  (y - mu) / (mu * (1 - mu)))
    dmu.deta <- dtheta.deta(mu, .link , earg = .earg)
    Step <- .step # This is another method of adjusting step lengths
    Step * c(w) * dl.dmu * dmu.deta
  }), list( .link = link,
            .step = step,
            .earg = earg ))),

  weight = eval(substitute(expression({
    x[, 1:p.lm] <- x.save[tt.index, 1:p.lm]  # Reinstate 

    for (ii in 1:plag) {
        temp <- theta2eta(y.save[tt.index-ii], .link , earg = .earg )


        x[, 1:p.lm] <- x[, 1:p.lm] -
                       x.save[tt.index-ii, 1:p.lm] * new.coeffs[ii + p.lm]
        x[, p.lm+ii] <- temp - x.save[tt.index-ii, 1:p.lm, drop = FALSE] %*%
                            new.coeffs[1:p.lm]
    }
    class(x) <- "matrix" # Added 27/2/02; 26/2/04

    if (iter == 1)
      old.coeffs <- new.coeffs 

    X.vlm.save <- lm2vlm.model.matrix(x, Hlist, xij = control$xij)

    vary <- switch( .link ,
                   identitylink = 1,
                   loge       = mu,
                   reciprocal = mu^2,
                   inverse    = mu^2,
                   mu * (1 - mu))
    c(w) * dtheta.deta(mu, link = .link , earg = .earg )^2 / vary
  }), list( .link = link,
            .earg = earg ))))
}






 if (FALSE) {
setClass(Class = "Coef.rrar", representation(
         "plag"    = "integer",
         "Ranks"   = "integer",
         "omega"   = "integer",
         "C"       = "matrix",
         "D"       = "matrix",
         "H"       = "matrix",
         "Z"       = "matrix",
         "Phi"     = "list",  # list of matrices
         "Ak1"     = "matrix"))



Coef.rrar <- function(object, ...) {
    result = new(Class = "Coef.rrar",
         "plag"     = object@misc$plag,
         "Ranks"    = object@misc$Ranks,
         "omega"    = object@misc$omega,
         "C"        = object@misc$C,
         "D"        = object@misc$D,
         "H"        = object@misc$H,
         "Z"        = object@misc$Z,
         "Phi"      = object@misc$Phi,
         "Ak1"      = object@misc$Ak1)
}





show.Coef.rrar <- function(object) {
  cat(object@plag)
}


setMethod("Coef", "rrar",
          function(object, ...)
          Coef(object, ...))




setMethod("show", "Coef.rrar",
          function(object)
          show.Coef.rrar(object))

}








if (FALSE)
 AR1.control <- function(criterion = "loglikelihood",
                         stepsize = 1,
                         maxit = 30, ...) {
  list(criterion = criterion,
       stepsize  = stepsize,
       maxit     = maxit)
}



if (TRUE)
 AR1.control <-
    function(half.stepsizing = FALSE,  # Avoids jittering very near the solution
                         ...) {
  list(half.stepsizing = half.stepsizing
       )
}



 AR1 <-
  function(ldrift = "identitylink",
           lsd  = "loge",
           lvar = "loge",
           lrho = "rhobit",
           idrift  = NULL,
           isd  = NULL,
           ivar = NULL,
           irho = NULL,
           imethod = 1,
           ishrinkage = 1,  # 0.90; unity means a constant
           type.likelihood = c("exact", "conditional"),
           var.arg = FALSE,  # TRUE,
           nodrift = FALSE,  # TRUE,
           almost1 = 0.99,
           zero = c(if (var.arg) "var" else "sd", "rho")  # "ARcoef1"
          ) {
  type.likelihood <- match.arg(type.likelihood,
                               c("exact", "conditional"))[1]

  if (!is.Numeric(almost1, length.arg = 1) || almost1 < 0.9 ||
      almost1 >= 1)
    stop("Bad input for argument 'almost1'")
  
  if (length(isd) && !is.Numeric(isd, positive = TRUE))
    stop("Bad input for argument 'isd'")
  
  if (length(ivar) && !is.Numeric(ivar, positive = TRUE))
    stop("Bad input for argument 'ivar'")
  
  if (length(irho) &&
     (!is.Numeric(irho) || any(abs(irho) > 1.0)))
    stop("Bad input for argument 'irho'")
  


  if (!is.logical(nodrift) ||
      length(nodrift) != 1)
    stop("argument 'nodrift' must be a single logical")

  if (!is.logical(var.arg) ||
      length(var.arg) != 1)
    stop("argument 'var.arg' must be a single logical")

  ismn <- idrift
  lsmn <- as.list(substitute(ldrift))
  esmn <- link2list(lsmn)
  lsmn <- attr(esmn, "function.name")     

  lsdv <- as.list(substitute(lsd))
  esdv <- link2list(lsdv)
  lsdv <- attr(esdv, "function.name")

  lvar  <- as.list(substitute(lvar))
  evar  <- link2list(lvar)
  lvar  <- attr(evar, "function.name")

  lrho <- as.list(substitute(lrho))
  erho <- link2list(lrho)
  lrho <- attr(erho, "function.name")     

  n.sc <- if (var.arg) "var" else "sd"
  l.sc <- if (var.arg) lvar else lsdv
  e.sc <- if (var.arg) evar else esdv


  new("vglmff", 
  blurb = c(ifelse(nodrift, "Two", "Three"),
            "-parameter autoregressive process of order-1\n\n",
            "Links:       ",
            if (nodrift) "" else
              paste(namesof("drift", lsmn, earg = esmn), ", ", sep = ""),
            namesof(n.sc , l.sc, earg = e.sc), ", ",
            namesof("rho", lrho, earg = erho), "\n",
            "Model:       Y_t = drift + rho * Y_{t-1} + error_{t},", "\n",
            "             where 'error_{2:n}' ~ N(0, sigma^2) independently",
            if (nodrift) ", and drift = 0" else "",
                          "\n",
            "Mean:        drift / (1 - rho)", "\n",
            "Correlation: rho = ARcoef1", "\n",
            "Variance:    sd^2 / (1 - rho^2)"),
  constraints = eval(substitute(expression({

    M1 <- 3 - .nodrift
    dotzero <- .zero
    eval(negzero.expression.VGAM)
  }), list( .zero = zero,
            .nodrift = nodrift ))),
  infos = eval(substitute(function(...) {
    list(M1 = 3 - nodrift, 
         Q1 = 1, 
         expected = TRUE, 
         multipleResponse = TRUE,
         type.likelihood = .type.likelihood ,
         ldrift = if ( .nodrift ) NULL else .lsmn ,
         edrift = if ( .nodrift ) NULL else .esmn ,
         lvar = .lvar ,
         lsd  = .lsdv ,
         evar = .evar ,
         esd  = .esdv ,
         lrho = .lrho ,
         erho = .erho ,
         almost1 = .almost1 ,
         zero = .zero )
  }, list( .lsmn = lsmn, .lvar = lvar, .lsdv = lsdv, .lrho = lrho,
           .esmn = esmn, .evar = evar, .esdv = esdv, .erho = erho,
           .type.likelihood = type.likelihood,
           .nodrift = nodrift,
           .almost1 = almost1, .zero = zero))),
  initialize = eval(substitute(expression({
    extra$M1 <- M1 <- 3 - .nodrift
    check <- w.y.check(w = w, y = y,
                       Is.positive.y = FALSE,
                       ncol.w.max = Inf,
                       ncol.y.max = Inf,
                       out.wy = TRUE,
                       colsyperw = 1, 
                       maximize = TRUE)
    w <- check$w
    y <- check$y
    if ( .type.likelihood == "conditional")
      w[1, ] <- 1.0e-8  # 1.0e-6

    
    NOS <- ncoly <- ncol(y)
    M <- M1*NOS
    var.names <- param.names("var",     NOS)
    sdv.names <- param.names("sd",      NOS)
    smn.names <- if ( .nodrift ) NULL else
                 param.names("drift",   NOS)
    rho.names <- param.names("rho",     NOS)

    mynames1 <- smn.names
    mynames2 <- if ( .var.arg ) var.names else sdv.names 
    mynames3 <- rho.names

    predictors.names <-
      c(if ( .nodrift ) NULL else
        namesof(smn.names, .lsmn , earg = .esmn , tag = FALSE),
        if ( .var.arg ) 
        namesof(var.names, .lvar , earg = .evar , tag = FALSE) else
        namesof(sdv.names, .lsdv , earg = .esdv , tag = FALSE),
        namesof(rho.names, .lrho , earg = .erho , tag = FALSE))
    predictors.names <- predictors.names[interleave.VGAM(M, M1 = M1)]


    if (!length(etastart)) {
      init.smn <- Init.mu(y = y, w = w, imethod = .imethod ,  # x = x,
                          imu = .ismn , ishrinkage = .ishrinkage ,
                          pos.only = FALSE)




      init.rho <- matrix(if (length( .irho )) .irho else 0.1,  # Dummy value
                         n, NOS, byrow = TRUE)
      init.sdv <- matrix(if (length( .isdv )) .isdv else 1.0,  # Dummy value
                         n, NOS, byrow = TRUE)
      init.var <- matrix(if (length( .ivar )) .ivar else 1.0,  # Dummy value
                         n, NOS, byrow = TRUE)
      for (jay in 1: NOS) {
        mycor <- cor(y[-1, jay], y[-n, jay])
        init.smn[-1, jay] <- init.smn[-1, jay] * (1 - mycor)
        if (!length( .irho ))
          init.rho[, jay] <- sign(mycor) * min(0.95, abs(mycor))
        if (!length( .ivar ))
          init.var[, jay] <- var(y[, jay]) * (1 - mycor^2)
        if (!length( .isdv ))
          init.sdv[, jay] <- sqrt(init.var[, jay])
      }  # for

      etastart <-
        cbind(if ( .nodrift ) NULL else
              theta2eta(init.smn, .lsmn , earg = .esmn ),
              if ( .var.arg ) 
              theta2eta(init.var, .lvar , earg = .evar ) else
              theta2eta(init.sdv, .lsdv , earg = .esdv ),
              theta2eta(init.rho, .lrho , earg = .erho ))
      etastart <- etastart[, interleave.VGAM(M, M1 = M1), drop = FALSE]
    }  # end of etastart
  }), list( .lsmn = lsmn, .lrho = lrho, .lsdv = lsdv, .lvar = lvar,
            .esmn = esmn, .erho = erho, .esdv = esdv, .evar = evar,
            .ismn = ismn, .irho = irho, .isdv = isd , .ivar = ivar,
            .type.likelihood = type.likelihood, .ishrinkage = ishrinkage,
            .var.arg = var.arg,
            .nodrift = nodrift,
            .imethod = imethod ))),
  linkinv = eval(substitute(function(eta, extra = NULL) {
    M1  <- 3 - .nodrift
    NOS <- ncol(eta)/M1
    ar.smn <- if ( .nodrift ) 0 else
      eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                .lsmn , earg = .esmn )
    ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                        .lrho , earg = .erho )
    ar.smn / (1 - ar.rho)
  }, list ( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
            .var.arg = var.arg, .type.likelihood = type.likelihood,
            .nodrift = nodrift,
            .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ))),
  last = eval(substitute(expression({
    if (any(abs(ar.rho) > 1)) 
      warning("Regularity conditions are violated at the final",
              "IRLS iteration, since 'abs(rho) > 1")

    M1 <- extra$M1

    temp.names <- c(mynames1, mynames2, mynames3)
    temp.names <- temp.names[interleave.VGAM(M1 * ncoly, M1 = M1)]

    misc$link <- rep( .lrho , length = M1 * ncoly)
    misc$earg <- vector("list", M1 * ncoly)
    names(misc$link) <-
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      if ( !( .nodrift ))
        misc$link[ M1*ii-2 ] <- .lsmn
      misc$link[ M1*ii-1 ] <- if ( .var.arg ) .lvar else .lsdv
      misc$link[ M1*ii   ] <- .lrho
      if ( !( .nodrift ))
        misc$earg[[M1*ii-2]] <- .esmn
      misc$earg[[M1*ii-1]] <- if ( .var.arg ) .evar else .esdv
      misc$earg[[M1*ii  ]] <- .erho
    }

    misc$type.likelihood <- .type.likelihood
    misc$var.arg <- .var.arg
    misc$M1 <- M1
    misc$expected <- TRUE
    misc$imethod <- .imethod
    misc$multipleResponses <- TRUE
    misc$nodrift <- .nodrift


  }), list( .lsmn = lsmn, .lrho = lrho, .lsdv = lsdv, .lvar = lvar,
            .esmn = esmn, .erho = erho, .esdv = esdv, .evar = evar,
            .irho = irho, .isdv = isd , .ivar = ivar,
            .nodrift = nodrift,
            .var.arg = var.arg, .type.likelihood = type.likelihood,
            .imethod = imethod ))),
  loglikelihood = eval(substitute(
    function(mu, y, w, residuals= FALSE, eta, 
             extra = NULL, summation = TRUE) {
      M1  <- 3 - .nodrift
      NOS <- ncol(eta)/M1

      if ( .var.arg ) {
        ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                            .lvar , earg = .evar )
      } else {
        ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                            .lsdv , earg = .esdv )
        ar.var <- ar.sdv^2
      }  
      ar.smn <- if ( .nodrift ) 0 else
        eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                  .lsmn , earg = .esmn )
      ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                          .lrho , earg = .erho )
    
      if (residuals) {
        stop("Loglikelihood not implemented yet to handle",
             "residuals.")
      } else {
        loglik.terms <- c(w) * dAR1(x = y,
                                    drift = ar.smn,
                                    var.error = ar.var,
                                    type.likelihood = .type.likelihood ,
                                    ARcoef1 = ar.rho, log = TRUE)
        loglik.terms <- as.matrix(loglik.terms)
        if (summation) {
          sum(if ( .type.likelihood == "exact") loglik.terms else
                   loglik.terms[-1, ])
        } else {
          loglik.terms
        }
      }
          
    }, list( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
             .var.arg = var.arg, .type.likelihood = type.likelihood,
             .nodrift = nodrift,
             .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ))),

  vfamily = c("AR1"),


  simslot = eval(substitute(
  function(object, nsim) {

      pwts <- if (length(pwts <- object@prior.weights) > 0)
                pwts else weights(object, type = "prior")
      if (any(pwts != 1))
        warning("ignoring prior weights")
      eta <- predict(object)
      fva <- fitted(object)      
      M1  <- 3 - .nodrift
      NOS <- ncol(eta)/M1

      if ( .var.arg ) {
        ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                            .lvar , earg = .evar )
      } else {
        ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                            .lsdv , earg = .esdv )
        ar.var <- ar.sdv^2
      }  
      ar.smn <- if ( .nodrift ) matrix(0, n, NOS) else
        eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                  .lsmn , earg = .esmn )
      ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                          .lrho , earg = .erho )

      ans <- array(0, c(nrow(eta), NOS, nsim))
    for (jay in 1:NOS) {
      ans[1, jay, ] <- rnorm(nsim, m = fva[1, jay],  # zz
                             sd = sqrt(ar.var[1, jay]))
      for (ii in 2:nrow(eta))
        ans[ii, jay, ] <- ar.smn[ii, jay] +
                          ar.rho[ii, jay] * ans[ii-1, jay, ] +
                          rnorm(nsim, sd = sqrt(ar.var[ii, jay]))
    }
    ans <- matrix(c(ans), c(nrow(eta) * NOS, nsim))
    ans
  }, list( .lsmn = lsmn, .lrho = lrho , .lsdv = lsdv, .lvar = lvar ,
           .var.arg = var.arg, .type.likelihood = type.likelihood,
           .nodrift = nodrift,
           .esmn = esmn, .erho = erho , .esdv = esdv, .evar = evar ))),




  deriv = eval(substitute(expression({
    M1  <- 3 - .nodrift
    NOS <- ncol(eta)/M1
    ncoly <- ncol(as.matrix(y))
    
    if ( .var.arg ) {
      ar.var <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                          .lvar , earg = .evar )
    } else {
      ar.sdv <- eta2theta(eta[, M1*(1:NOS) - 1, drop = FALSE],
                          .lsdv , earg = .esdv )
      ar.var <- ar.sdv^2
    }  
    
    ar.smn <- if ( .nodrift ) matrix(0, n, NOS) else
      eta2theta(eta[, M1*(1:NOS) - 2, drop = FALSE],
                .lsmn , earg = .esmn )
    ar.rho <- eta2theta(eta[, M1*(1:NOS)    , drop = FALSE],
                        .lrho , earg = .erho )


    temp1 <- y[-1, , drop = FALSE] -
             ar.smn[-1, , drop = FALSE] -
             y[-n, , drop = FALSE] * ar.rho[-1, , drop = FALSE]
    temp5 <- 1 - ar.rho^2

    dl.dsmn <- rbind(rep(0, length = ncoly),
                     temp1 / ar.var[-1, , drop = FALSE])
    if ( .var.arg ) {
      dl.dvar <- rbind(rep(0, length = ncoly),
                       0.5 * (temp1 / ar.var[-1, , drop = FALSE])^2 -
                       0.5 / ar.var[-1, , drop = FALSE])
    } else {
      dl.dsdv <- rbind(rep(0, length = ncoly),
                       temp1^2 / (ar.sdv[-1, , drop = FALSE])^3 -
                       1 / ar.sdv[-1, , drop = FALSE])
    }
    dl.drho <- rbind(rep(0, length = ncoly),
                     y[-n, , drop = FALSE] *
                     temp1 / ar.var[-1, , drop = FALSE])
    dl.dsmn[1, ] <- (y[1, ] - mu[1, ]) * (1 + ar.rho[1, ]) / ar.var[1, ]
    dl.drho[1, ] <- ar.rho[1, ] * (y[1, ] - mu[1, ])^2 / ar.var[1, ] -
                    ar.rho[1, ] / temp5[1, ]
    if ( .var.arg ) {
      dl.dvar[1, ] <- -0.5 / ar.var[1, ] +
                       0.5 *  temp5[1, ] *
                       ((y[1, ] - mu[1, ]) / ar.var[1, ])^2
    } else {
      dl.dsdv[1, ] <- -1 / ar.sdv[1, ] +
                       temp5[1, ] * (y[1, ] - mu[1, ])^2 / (ar.sdv[1, ])^3
    }



        
    dsmn.deta <- dtheta.deta(ar.smn, .lsmn , earg = .esmn )
    drho.deta <- dtheta.deta(ar.rho, .lrho , earg = .erho )
    if ( .var.arg ) {
      dvar.deta <- dtheta.deta(ar.var, .lvar , earg = .evar )
    } else {
      dsdv.deta <- dtheta.deta(ar.sdv, .lsdv , earg = .esdv )
    }
    myderiv <-
      c(w) * cbind(if ( .nodrift ) NULL else dl.dsmn * dsmn.deta,
                   if ( .var.arg ) dl.dvar * dvar.deta else
                                   dl.dsdv * dsdv.deta,
                   dl.drho * drho.deta)
                   
    myderiv[, interleave.VGAM(M, M1 = M1)]
  }), list( .lsmn = lsmn, .lrho = lrho, .lsdv = lsdv, .lvar = lvar,
            .esmn = esmn, .erho = erho, .esdv = esdv, .evar = evar,
            .nodrift = nodrift,
            .var.arg = var.arg, .type.likelihood = type.likelihood ))),
  weight = eval(substitute(expression({
    if ( .var.arg ) {
      ned2l.dvar    <- 0.5 / ar.var^2
      ned2l.drhovar <-  matrix(0, n, ncoly)
      ned2l.drhovar[1, ] <- .almost1 * ar.rho[1, ] / (ar.var[1, ] *
                                                      temp5[1, ])
    } else {
      ned2l.dsdv    <- 2 / ar.var
      ned2l.drhosdv <-  matrix(0, n, ncoly)
      ned2l.drhosdv[1, ] <- 2 *
                            .almost1 * ar.rho[1, ] / (ar.sdv[1, ] *
                                                      temp5[1, ])
    }
    
    if ( !( .nodrift )) {
      ned2l.dsmn <- 1 / ar.var
      ned2l.dsmn[1, ] <- (1 + ar.rho[1, ]) / ((1 - ar.rho[1, ]) *
                                              ar.var[1, ])
      ned2l.dsmnrho <- mu / ar.var
      ned2l.dsmnrho[1, ] <- 0
    }

    ned2l.drho <- ((    mu[-n, , drop = FALSE])^2 +
                    ar.var[-n, , drop = FALSE] /
                     temp5[-n, , drop = FALSE]) / ar.var[-1, , drop = FALSE]
    ned2l.drho <- rbind(0, ned2l.drho)
    ned2l.drho[1, ]  <- 2 * (ar.rho[1, ] / temp5[1, ])^2
    
    ncol.wz <- M + (M - 1) + ifelse( .nodrift , 0, M - 2)
    wz <- matrix(0, n, ncol.wz)

    if ( !( .nodrift ))
      wz[, M1*(1:NOS) - 2] <- ned2l.dsmn * dsmn.deta^2
    wz[, M1*(1:NOS) - 1] <-
      if ( .var.arg ) ned2l.dvar * dvar.deta^2 else
                      ned2l.dsdv * dsdv.deta^2
    wz[, M1*(1:NOS)    ] <- ned2l.drho * drho.deta^2

    if ( !( .nodrift ))
      wz[, M1*(1:NOS) + M + (M - 1) - M1 + 1] <- ned2l.dsmnrho *
                                                 dsmn.deta * drho.deta
    wz[, M1*(1:NOS) + M - 1] <-
      if ( .var.arg ) ned2l.drhovar * drho.deta * dvar.deta else
                      ned2l.drhosdv * drho.deta * dsdv.deta

    wz <- w.wz.merge(w = w, wz = wz, n = n, M = ncol.wz, ndepy = NOS)

    wz
  }), list( .var.arg = var.arg, .type.likelihood = type.likelihood,
            .nodrift = nodrift,
            .almost1 = almost1)))
  )
}  # End of function 'AR1'










dAR1 <- function(x,
                 drift = 0,  # Stationarity is the default
                 var.error = 1, ARcoef1 = 0.0,
                 type.likelihood = c("exact", "conditional"),
                 log = FALSE) {

  type.likelihood <- match.arg(type.likelihood,
                               c("exact", "conditional"))[1]

  is.vector.x <- is.vector(x)

  x <- as.matrix(x)
  drift <- as.matrix(drift)
  var.error <- as.matrix(var.error)
  ARcoef1 <- as.matrix(ARcoef1)
  LLL <- max(nrow(x), nrow(drift), nrow(var.error), nrow(ARcoef1))
  UUU <- max(ncol(x), ncol(drift), ncol(var.error), ncol(ARcoef1))
  x          <- matrix(x,         LLL, UUU)
  drift    <- matrix(drift,   LLL, UUU)
  var.error  <- matrix(var.error, LLL, UUU)
  rho        <- matrix(ARcoef1,   LLL, UUU)

  if (any(abs(rho) > 1))
    warning("Values of argument 'ARcoef1' are greater ",
            "than 1 in absolute value")

  if (!is.logical(log.arg <- log) || length(log) != 1) 
    stop("Bad input for argument 'log'.")
  rm(log)

  ans <- matrix(0.0, LLL, UUU)
  var.noise <- var.error / (1 - rho^2)

  ans[ 1, ] <- dnorm(x    = x[1, ],
                     mean = drift[ 1, ] / (1 - rho[1, ]), 
                     sd   = sqrt(var.noise[1, ]), log = log.arg)
  ans[-1, ] <- dnorm(x    = x[-1, ],
                     mean = drift[-1, ] + rho[-1, ] * x[-nrow(x), ],
                     sd   = sqrt(var.noise[-1, ]), log = log.arg)

  if (type.likelihood == "conditional")
    ans[1, ] <- NA
  
  if (is.vector.x) as.vector(ans) else ans
}






