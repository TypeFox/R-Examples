# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.











 SURff <-
  function(mle.normal = FALSE,
           divisor = c("n", "n-max(pj,pk)", "sqrt((n-pj)*(n-pk))"),
           parallel = FALSE, 
           Varcov = NULL,
           matrix.arg = FALSE) {








  apply.parint <- TRUE


  lmean <- "identitylink"
  lsdev <- "loge"
  emean <- list()
  esdev <- list()


  if (!is.logical(mle.normal) ||
      length(mle.normal) != 1)
    stop("argument 'mle.normal' must be a single logical")

  if (!is.logical(apply.parint) ||
      length(apply.parint) != 1)
    stop("argument 'apply.parint' must be a single logical")




  divisor <- match.arg(divisor,
      c("n", "n-max(pj,pk)", "sqrt((n-pj)*(n-pk))"))[1]

  if (mle.normal && divisor != "n")
    warning("MLE requires 'n' as the value of argument 'divisor'. ",
            "The solution will probably not be the MLE")


  ret.ff <-
  new("vglmff",
  blurb = c("Seemingly unrelated regressions"),
  constraints = eval(substitute(expression({
    constraints <- cm.VGAM(matrix(1, M, 1), x = x,
                           bool = .parallel , 
                           constraints = constraints,
                           apply.int = .apply.parint )
  }), list( .parallel = parallel,
            .apply.parint = apply.parint ))),




  infos = eval(substitute(function(...) {
    list(M1 = 1,  # zz???
         Q1 = 1,
         parallel = .parallel ,
         expected = TRUE,
         multipleResponses = TRUE,
         parameters.names = as.character(NA))
  }, list( .parallel = parallel ))),

  initialize = eval(substitute(expression({

    if (!is.matrix(y) || ncol(y) == 1)
      stop("response must be a matrix with at least 2 columns")
    ncoly <- ncol(y)

   if (is.logical( .parallel ) &&
       .parallel &&
       !all(as.logical(trivial.constraints(constraints))))
     warning("setting 'parallel = TRUE' with nontrivial constraints may not ",
             "make sense")

   temp5 <-
    w.y.check(w = w, y = y,
              ncol.w.min = 1,
              ncol.w.max = 1,
              ncol.y.max = Inf,
              Is.integer.y = FALSE,
              Is.positive.y = FALSE,
              out.wy = TRUE,
              colsyperw = ncoly,
              maximize = TRUE)
    w <- temp5$w
    y <- temp5$y

    if (!all(w[1, 1] == w))
      stop("all prior 'weights' must currently have equal values")


    ncoly <- ncol(y)
    M1 <- 1
    extra$ncoly <- ncoly
    extra$M1 <- M1
    M <- M1 * ncoly


    predictors.names <- if (!length(ddd <- dimnames(y)[[2]]))
        paste("Y", 1:M, sep = "") else ddd






    extra$wz <- matrix(1, nrow(x), M)


    if (!length(etastart)) {
      etastart <- matrix(0, n, M)


      Hlist.early <- process.constraints(constraints, x, M,
                                         specialCM = specialCM)
      X.vlm.early  <- lm2vlm.model.matrix(x, Hlist.early,
                                          xij = control$xij,
                                          Xm2 = Xm2)

      Hmatrices <- matrix(c(unlist(Hlist.early)), nrow = M)
      jay.index <- 1:ncol(Hmatrices)


      extra$ncols.X.lm <- numeric(ncoly)
      for (jay in 1:ncoly) {

        X.lm.jay <- vlm2lm.model.matrix(x.vlm = X.vlm.early,
                                        Hlist = Hlist.early,
                                        which.linpred = jay, M = M)

        extra$ncols.X.lm[jay] <- ncol(X.lm.jay)

        etastart[, jay] <- y[, jay] -
                           lsfit(x = X.lm.jay, y = y[, jay],
                                 wt = c(w), intercept = FALSE)$residuals
      }  # jay
    }  # !length(etastart)
  }), list(
            .parallel = parallel 
          ))),
  linkinv = function(eta, extra = NULL) eta, 
  last = eval(substitute(expression({

    M1 <- extra$M1
    misc$link <- c(rep( .lmean , length = ncoly))
    temp.names <- predictors.names
    names(misc$link) <- temp.names

    misc$earg <- vector("list", M1 * ncoly)
    names(misc$earg) <- temp.names
    for (ii in 1:ncoly) {
      misc$earg[[M1*ii]] <- .emean
    }
    names(misc$earg) <- temp.names

    misc$M1 <- M1
    misc$expected <- TRUE
    misc$divisor <- .divisor
    misc$values.divisor <- round(n / ratio.df)

  }), list( .lmean = lmean, .lsdev = lsdev,
            .emean = emean, .esdev = esdev,
            .divisor = divisor
          ))),

  vfamily = "SURff",


  deriv = eval(substitute(expression({
    mymu <- eta
    iam.indices <- iam(NA, NA, M = M, both = TRUE)
    resmat <- y - mymu
    Sigma.elts <- colMeans(resmat[, iam.indices$row.index] *
                           resmat[, iam.indices$col.index])

    if ( .divisor != "n") {
      ratio.df <- n / switch( .divisor ,
        "n-max(pj,pk)" = n - pmax(extra$ncols.X.lm[iam.indices$row.index],
                                  extra$ncols.X.lm[iam.indices$col.index]),
        "sqrt((n-pj)*(n-pk))" =
        sqrt((n - extra$ncols.X.lm[iam.indices$row.index]) *
             (n - extra$ncols.X.lm[iam.indices$col.index])),
        stop("argument 'divisor' unmatched"))
      Sigma.elts <- Sigma.elts * ratio.df
    } else {
      ratio.df <- rep(1, length = M*(M+1)/2)
    }

    Sigma.mat <- matrix(0, M, M)
    Sigma.mat[cbind(iam.indices$row.index,
                    iam.indices$col.index)] <- Sigma.elts
    Sigma.mat[cbind(iam.indices$col.index,
                    iam.indices$row.index)] <- Sigma.elts

    invSigma.mat <- chol2inv(chol(Sigma.mat))


    temp3 <- matrix(invSigma.mat[cbind(iam.indices$row.index,
                                       iam.indices$col.index)],
                    M*(M+1)/2, n)
    dl.dmu <- mux22(temp3, y - mymu, M = M,
                    upper = FALSE, as.matrix = TRUE)
    dmu.deta <- dtheta.deta(mymu,   .lmean , earg = .emean )

    c(w) * dl.dmu * dmu.deta
  }), list( .lmean = lmean,
            .emean = emean,
            .divisor = divisor ))),


  weight = eval(substitute(expression({


    if (length( .Varcov )) {
      Sigma.mat <- if ( .matrix.arg ) .Varcov else {
                     temp.vec <- rep( .Varcov , len = M*(M+1)/2)
                     temp.mat <- matrix(0, M, M)
                     temp.mat[cbind(iam.indices$col.index,
                                    iam.indices$row.index)] <- temp.vec
                     temp.mat[cbind(iam.indices$row.index,
                                    iam.indices$col.index)] <- temp.vec
                     temp.mat
                   }
      invSigma.mat <- chol2inv(chol(Sigma.mat))
    }


    wz <-
    extra$wz <- c(w) * matrix(invSigma.mat[cbind(iam.indices$col.index,
                                                 iam.indices$row.index)],
                              n, M*(M+1)/2, byrow = TRUE)
    extra$Sigma.mat <- Sigma.mat
    extra$invSigma.mat <- invSigma.mat

    wz
  }), list( .divisor = divisor,
            .Varcov = Varcov,
            .matrix.arg = matrix.arg ))))



  if (mle.normal) {


    ret.ff@loglikelihood <-
      function(mu, y, w, residuals = FALSE, eta, extra = NULL,
               summation = TRUE) {

      if (!summation)
        stop("cannot handle 'summation = FALSE' yet")


      M <- if (is.matrix(y)) ncol(y) else 1
      n <- if (is.matrix(y)) nrow(y) else length(y)

      wz <- extra$wz

      temp1 <- ResSS.vgam(y-mu, wz = wz, M = M)
      onewz <- if (length(extra$invSigma.mat))
                 extra$invSigma.mat else
                 (m2a(wz[1, , drop = FALSE], M = M))[,, 1]  # M x M


      logdet <- determinant(onewz)$modulus
      logretval <- -0.5 * temp1 + 0.5 * n * logdet -
                   n * (M / 2) * log(2*pi)
      logretval
    }
  }

  ret.ff
}












