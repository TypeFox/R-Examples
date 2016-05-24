# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.







replace.constraints <- function(Hlist, cm, index) {

  for (iii in index)
    Hlist[[iii]] <- cm
  Hlist
}




 valt.control <- function(
                 Alphavec = c(2, 4, 6, 9, 12, 16, 20, 25, 30, 40, 50,
                              60, 80, 100, 125, 2^(8:12)),
                 Criterion = c("ResSS", "coefficients"),
                 Linesearch = FALSE,
                 Maxit = 7,
                 Suppress.warning = TRUE,
                 Tolerance = 1e-7, ...) {

  if (mode(Criterion) != "character" && mode(Criterion) != "name")
    Criterion <- as.character(substitute(Criterion))
  Criterion <- match.arg(Criterion, c("ResSS", "coefficients"))[1]

  list(Alphavec = Alphavec,
       Criterion = Criterion, 
       Linesearch = Linesearch,
       Maxit = Maxit,
       Suppress.warning = Suppress.warning,
       Tolerance = Tolerance)
}




qrrvglm.xprod <- function(numat, Aoffset, Quadratic, I.tolerances) {
  Rank <- ncol(numat)
  moff <- NULL
  ans <- if (Quadratic) {
           index <- iam(NA, NA, M = Rank, diag = TRUE, both = TRUE) 
           temp1 <- cbind(numat[, index$row] * numat[, index$col])
           if (I.tolerances) {
             moff <- 0
             for (ii in 1:Rank)
               moff <- moff - 0.5 * temp1[, ii]
           }
           cbind(numat, if (I.tolerances) NULL else temp1)
  } else {
    as.matrix(numat)
  }
  list(matrix = if (Aoffset > 0) ans else ans[, -(1:Rank), drop = FALSE],
       offset = moff)
}



 valt <- function(x, z, U, Rank = 1,
                  Hlist = NULL, 
                  Cinit = NULL,
                  Alphavec = c(2, 4, 6, 9, 12, 16, 20, 25, 30, 40, 50,
                               60, 80, 100, 125, 2^(8:12)),
                  Criterion = c("ResSS", "coefficients"),
                  Crow1positive = rep(TRUE, length.out = Rank),
                  colx1.index,
                  Linesearch = FALSE,
                  Maxit = 20, 
                  str0 = NULL,
                  sd.Cinit = 0.02,
                  Suppress.warning = FALSE,
                  Tolerance = 1e-6, 
                  trace = FALSE,
                  xij = NULL) {




                 


  if (mode(Criterion) != "character" && mode(Criterion) != "name")
    Criterion <- as.character(substitute(Criterion))
  Criterion <- match.arg(Criterion, c("ResSS", "coefficients"))[1]

  if (any(diff(Alphavec) <= 0))
    stop("'Alphavec' must be an increasing sequence") 

  if (!is.matrix(z))
    z <- as.matrix(z)
  n <- nrow(z)
  M <- ncol(z)
  if (!is.matrix(x))
    x <- as.matrix(x)

  colx2.index <- if (is.null(colx1.index)) 1:ncol(x) else
                 (1:ncol(x))[-colx1.index]

  p1 <- length(colx1.index)
  p2 <- length(colx2.index)
  p  <- p1 + p2
  if (!p2)
    stop("'p2', the number of variables for the ",
         "reduced-rank regression, must be > 0")

  if (!length(Hlist)) {
    Hlist <- replace.constraints(vector("list", p), diag(M), 1:p)
  }

  dU <- dim(U)
  if (dU[2] != n)
    stop("input unconformable")

  clist2 <- replace.constraints(vector("list", Rank+p1),
                                if (length(str0))
                                diag(M)[, -str0, drop = FALSE] else
                                diag(M), 1:Rank)
  if (p1) {
    for (kk in 1:p1)
      clist2[[Rank+kk]] <- Hlist[[colx1.index[kk]]]
  }

  if (is.null(Cinit))
    Cinit <- matrix(rnorm(p2*Rank, sd = sd.Cinit), p2, Rank)

  fit <- list(ResSS = 0)  # Only for initial old.crit below

  C <- Cinit  # This is input for the main iter loop
  old.crit <- switch(Criterion, coefficients = C, ResSS = fit$ResSS)

  recover <- 0  # Allow a few iterations between different line searches 
  for (iter in 1:Maxit) {
    iter.save <- iter

      latvar.mat <- x[, colx2.index, drop = FALSE] %*% C
      new.latvar.model.matrix <- cbind(latvar.mat,
                                       if (p1) x[, colx1.index] else NULL)
      fit <- vlm.wfit(xmat = new.latvar.model.matrix, z, Hlist = clist2,
                      U = U, matrix.out = TRUE, is.vlmX = FALSE,
                      ResSS = FALSE, qr = FALSE, xij = xij)
      A <- t(fit$mat.coef[1:Rank, , drop = FALSE])

      clist1 <- replace.constraints(Hlist, A, colx2.index)
      fit <- vlm.wfit(xmat = x, z, Hlist = clist1, U = U,
                      matrix.out = TRUE, is.vlmX = FALSE,
                      ResSS = TRUE, qr = FALSE, xij = xij)
      C <- fit$mat.coef[colx2.index, , drop = FALSE] %*% A %*%
           solve(t(A) %*% A)

      numat <- x[, colx2.index, drop = FALSE] %*% C
      evnu <- eigen(var(numat))
      temp7 <- if (Rank > 1) evnu$vector %*% diag(evnu$value^(-0.5)) else
                 evnu$vector %*% evnu$value^(-0.5)
      C <- C %*% temp7
      A <- A %*% t(solve(temp7))
      temp8 <- crow1C(cmat = C, Crow1positive, amat = A)
      C <- temp8$cmat
      A <- temp8$amat


      ratio <-
          switch(Criterion,
                 coefficients = max(abs(C - old.crit) / (
                                    Tolerance + abs(C))),
                 ResSS = max(abs(fit$ResSS - old.crit) / (
                              Tolerance + fit$ResSS)))

        if (trace) {
          cat("   Alternating iteration", iter,
              ",   Convergence criterion  = ", format(ratio), "\n")
          if (!is.null(fit$ResSS))
              cat("    ResSS  = ", fit$ResSS, "\n")
          flush.console()
      }

      if (ratio < Tolerance) {
        if (!Linesearch || (Linesearch && iter >= 3))
          break
      } else if (iter == Maxit && !Suppress.warning) {
        warning("did not converge")
      }

      fini.linesearch <- FALSE
      if (Linesearch && iter - recover >= 2) {
          xnew <- C

          direction1 <- (xnew - xold)  # / sqrt(1 + sum((xnew-xold)^2))
          ftemp <- fit$ResSS  # Most recent objective function 
          use.alpha <- 0  # The current step relative to (xold, yold)
          for (itter in 1:length(Alphavec)) {
            CC <- xold + Alphavec[itter] * direction1

            try.latvar.mat <- x[, colx2.index, drop = FALSE] %*% CC
            try.new.latvar.model.matrix <-
              cbind(try.latvar.mat,
                    if (p1) x[, colx1.index] else NULL)

            try <- vlm.wfit(xmat = try.new.latvar.model.matrix, z,
                            Hlist = clist2, U = U, matrix.out = TRUE,
                            is.vlmX = FALSE, ResSS = TRUE, qr = FALSE,
                            xij = xij)
            if (try$ResSS < ftemp) {
              use.alpha <- Alphavec[itter]
              fit <- try 
              ftemp <- try$ResSS
              C <- CC 
              A <- t(fit$mat.coef[1:Rank, , drop = FALSE])
              latvar.mat <- x[, colx2.index, drop = FALSE] %*% C
              recover <- iter  # Give it some altg iters to recover
            } else {
              if (trace && use.alpha > 0) {
                cat("    Finished line search using Alpha  = ",
                    use.alpha, "\n")
                flush.console()
              }
              fini.linesearch <- TRUE
            }
          if (fini.linesearch) break 
        }  # End of itter loop 
    }

    xold <- C # Do not take care of drift
    old.crit <- switch(Criterion,
                       coefficients = C,
                       ResSS = fit$ResSS)
  }  # End of iter loop

  list(A = A,
       C = C,
       fitted = fit$fitted,
       new.coeffs = fit$coef,
       ResSS = fit$ResSS)
}




 lm2qrrvlm.model.matrix <-
  function(x, Hlist, C, control, assign = TRUE,
           no.thrills = FALSE) {

    Rank <- control$Rank
    colx1.index <- control$colx1.index
    Quadratic <- control$Quadratic
    Dzero <- control$Dzero
    Corner <- control$Corner
    I.tolerances <- control$I.tolerances

    M <- nrow(Hlist[[1]])
    p1 <- length(colx1.index)
    combine2 <- c(control$str0,
                  if (Corner) control$Index.corner else NULL)

    Qoffset <- if (Quadratic) ifelse(I.tolerances, 0, sum(1:Rank)) else 0
    NoA <- length(combine2) == M    # No unknown parameters in A
    clist2 <- if (NoA) {
        Aoffset <- 0
        vector("list", Aoffset+Qoffset+p1)
    } else {
      Aoffset <- Rank
      replace.constraints(vector("list", Aoffset+Qoffset+p1),
           if (length(combine2)) diag(M)[, -combine2, drop = FALSE] else diag(M),
           1:Rank)  # If Corner then does not contain \bI_{Rank}
    }
    if (Quadratic && !I.tolerances)
      clist2 <- replace.constraints(clist2,
          if (control$eq.tolerances)
              matrix(1, M, 1) - eijfun(Dzero, M) else {
          if (length(Dzero)) diag(M)[,-Dzero, drop = FALSE] else diag(M)},
          Aoffset + (1:Qoffset))
    if (p1)
      for (kk in 1:p1)
        clist2[[Aoffset+Qoffset+kk]] <- Hlist[[colx1.index[kk]]]
    if (!no.thrills) {
      i63 <- iam(NA, NA, M=Rank, both = TRUE)
      names(clist2) <- c(
             if (NoA) NULL else paste("(latvar", 1:Rank, ")", sep = ""), 
             if (Quadratic && Rank == 1 && !I.tolerances)
                 "(latvar^2)" else 
             if (Quadratic && Rank>1 && !I.tolerances)
                 paste("(latvar", i63$row, ifelse(i63$row == i63$col, "^2",
                 paste("*latvar", i63$col, sep = "")), ")", sep = "") else NULL,
             if (p1) names(colx1.index) else NULL)
    }

    latvar.mat <- x[, control$colx2.index, drop = FALSE] %*% C


    tmp900 <- qrrvglm.xprod(latvar.mat, Aoffset, Quadratic, I.tolerances)
    new.latvar.model.matrix <- cbind(tmp900$matrix,
                                     if (p1) x[,colx1.index] else NULL)
    if (!no.thrills)
      dimnames(new.latvar.model.matrix) <- list(dimnames(x)[[1]],
                                                names(clist2))

    if (assign) {
      asx <- attr(x, "assign")
      asx <- vector("list", ncol(new.latvar.model.matrix))
      names(asx) <- names(clist2)
      for (ii in 1:length(names(asx))) {
        asx[[ii]] <- ii
      }
      attr(new.latvar.model.matrix, "assign") <- asx
    }

    if (no.thrills)
      list(new.latvar.model.matrix = new.latvar.model.matrix,
           constraints = clist2,
           offset = tmp900$offset) else
      list(new.latvar.model.matrix = new.latvar.model.matrix,
           constraints = clist2,
           NoA = NoA,
           Aoffset = Aoffset,
           latvar.mat = latvar.mat,
           offset = tmp900$offset)
}



valt.2iter <- function(x, z, U, Hlist, A, control) {


  clist1 <- replace.constraints(Hlist, A, control$colx2.index)
  fit <- vlm.wfit(xmat = x, z, Hlist = clist1, U = U, matrix.out = TRUE, 
                  is.vlmX = FALSE, ResSS = TRUE, qr = FALSE, xij = control$xij)
  C <- fit$mat.coef[control$colx2.index, , drop = FALSE] %*%
       A %*% solve(t(A) %*% A)

  list(A = A, C = C,
       fitted = fit$fitted, new.coeffs = fit$coef,
       Hlist = clist1, ResSS = fit$ResSS)
}



valt.1iter <- function(x, z, U, Hlist, C, control,
                      lp.names = NULL, nice31 = FALSE,
                      MSratio = 1) {

    Rank <- control$Rank
    Quadratic <- control$Quadratic
    Index.corner <- control$Index.corner
    p1 <- length(control$colx1.index)
    M <- ncol(zedd <- as.matrix(z))
    NOS <- M / MSratio
    Corner <- control$Corner
    I.tolerances <- control$I.tolerances

    Qoffset <- if (Quadratic) ifelse(I.tolerances, 0, sum(1:Rank)) else 0
    tmp833 <- lm2qrrvlm.model.matrix(x = x, Hlist = Hlist, C = C,
                                     control = control)
    new.latvar.model.matrix <- tmp833$new.latvar.model.matrix 
    clist2 <- tmp833$constraints # Does not contain \bI_{Rank}
    latvar.mat <- tmp833$latvar.mat
    if (Corner)
        zedd[,Index.corner] <- zedd[,Index.corner] - latvar.mat

    if (nice31 && MSratio == 1) {
      fit <- list(mat.coef = NULL, fitted.values = NULL, ResSS = 0)

      clist2 <- NULL # for vlm.wfit

      i5 <- rep(0, length.out = MSratio)
      for (ii in 1:NOS) {
        i5 <- i5 + 1:MSratio

        tmp100 <- vlm.wfit(xmat = new.latvar.model.matrix,
                           zedd[, i5, drop = FALSE],
                           Hlist = clist2,
                           U = U[i5,, drop = FALSE],
                           matrix.out = TRUE,
                           is.vlmX = FALSE, ResSS = TRUE,
                           qr = FALSE,
                           Eta.range = control$Eta.range,
                           xij = control$xij,
                           lp.names = lp.names[i5])
        fit$ResSS <- fit$ResSS + tmp100$ResSS
        fit$mat.coef <- cbind(fit$mat.coef, tmp100$mat.coef)
        fit$fitted.values <- cbind(fit$fitted.values,
                                   tmp100$fitted.values)
      }
    } else {
      fit <- vlm.wfit(xmat = new.latvar.model.matrix,
                      zedd, Hlist = clist2, U = U,
                      matrix.out = TRUE,
                      is.vlmX = FALSE, ResSS = TRUE, qr = FALSE,
                      Eta.range = control$Eta.range,
                      xij = control$xij, lp.names = lp.names)
    }
    A <- if (tmp833$NoA) matrix(0, M, Rank) else
        t(fit$mat.coef[1:Rank,, drop = FALSE])
    if (Corner)
        A[Index.corner,] <- diag(Rank)     

    B1 <- if (p1)
      fit$mat.coef[-(1:(tmp833$Aoffset+Qoffset)),, drop = FALSE] else
      NULL
    fv <- as.matrix(fit$fitted.values)
    if (Corner)
        fv[,Index.corner] <- fv[,Index.corner] + latvar.mat
    Dmat <- if (Quadratic) {
            if (I.tolerances) {
                tmp800 <- matrix(0, M, Rank*(Rank+1)/2)
                tmp800[if (MSratio == 2) c(TRUE, FALSE) else
                       TRUE, 1:Rank] <- -0.5
                tmp800
            } else 
                t(fit$mat.coef[(tmp833$Aoffset+1):
                  (tmp833$Aoffset+Qoffset),, drop = FALSE])
    } else
        NULL

    list(Amat = A, B1 = B1, Cmat = C, Dmat = Dmat,
         fitted = if (M == 1) c(fv) else fv,
         new.coeffs = fit$coef, constraints = clist2, ResSS = fit$ResSS,
         offset = if (length(tmp833$offset)) tmp833$offset else NULL)
}





rrr.init.expression <- expression({
    if (length(control$Quadratic) && control$Quadratic)
      copy.X.vlm <- TRUE




  if (function.name %in% c("cqo", "cao")) {

    modelno <- switch(family@vfamily[1], "poissonff" = 2,
              "quasipoissonff" = 2, "quasipoisson" = 2,
              "binomialff" = 1, "quasibinomialff" = 1,
              "quasibinomial" = 1, "negbinomial" = 3,
              "gamma2" = 5, "gaussianff" = 8,
              0)  # stop("cannot fit this model using fast algorithm")
    if (modelno == 1) modelno = get("modelno", envir = VGAMenv)
    rrcontrol$modelno = control$modelno = modelno
    if (modelno == 3 || modelno == 5) {


      M <- 2 * ifelse(is.matrix(y), ncol(y), 1)
        control$str0 <-
      rrcontrol$str0 <- seq(from = 2, to = M, by = 2)  # Handles A
        control$Dzero <-
      rrcontrol$Dzero <- seq(from = 2, to = M, by = 2)  # Handles D


    }
  } else {
    modelno <- 0  # Any value will do as the variable is unused.
  }


})



rrr.alternating.expression <- expression({

    alt <- valt(x, z, U, Rank = Rank,
                Hlist = Hlist,
                Cinit = rrcontrol$Cinit,
                Criterion = rrcontrol$Criterion,
                colx1.index = rrcontrol$colx1.index,
                Linesearch = rrcontrol$Linesearch,
                Maxit = rrcontrol$Maxit,
                str0 = rrcontrol$str0,
                sd.Cinit = rrcontrol$sd.Cinit,
                Suppress.warning = rrcontrol$Suppress.warning,
                Tolerance = rrcontrol$Tolerance,
                trace = trace,
                xij = control$xij)  # This is subject to drift in A and C

    ans2 <- rrr.normalize(rrcontrol = rrcontrol, A=alt$A, C=alt$C, x = x)

    Amat <- ans2$A           # Fed into Hlist below (in rrr.end.expression)
    tmp.fitted <- alt$fitted # Also fed; was alt2$fitted 

    rrcontrol$Cinit <- ans2$C   # For next valt() call

    eval(rrr.end.expression)    # Put Amat into Hlist, and create new z
})



  adjust.Dmat.expression <- function(Mmat, Rank, Dmat, M) {

    if (length(Dmat)) {
      ind0 <- iam(NA, NA, both = TRUE, M = Rank)
      for (kay in 1:M) {
        elts <- Dmat[kay, , drop = FALSE]  # Manual recycling
        if (length(elts) < Rank)
          elts <- matrix(elts, 1, Rank)
        Dk <- m2a(elts, M = Rank)[, , 1]
        Dk <- matrix(Dk, Rank, Rank)
        Dk <- t(Mmat) %*% Dk  %*% Mmat  # 20030822; Not diagonal in general.
        Dmat[kay, ] <- Dk[cbind(ind0$row.index[1:ncol(Dmat)],
                                ind0$col.index[1:ncol(Dmat)])] 
      }
    }
    Dmat
  }



rrr.normalize <- function(rrcontrol, A, C, x, Dmat = NULL) {



    colx2.index <- rrcontrol$colx2.index
    Rank <- rrcontrol$Rank
    Index.corner <- rrcontrol$Index.corner
    M <- nrow(A)
    C.old <- C

    if (rrcontrol$Corner) {
      tmp87 <- A[Index.corner,, drop = FALSE]
      Mmat <- solve(tmp87)  # The normalizing matrix
      C <- C %*% t(tmp87)
      A <- A %*% Mmat
      A[Index.corner,] <- diag(Rank)  # Make sure 

      Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                     Dmat = Dmat, M = M)
    }

    if (rrcontrol$Svd.arg) {
      temp <- svd(C %*% t(A))
      if (!is.matrix(temp$v))
        temp$v <- as.matrix(temp$v) 
      C <- temp$u[, 1:Rank, drop = FALSE] %*%
           diag(temp$d[1:Rank]^(1-rrcontrol$Alpha), nrow = Rank)
      A <- diag(temp$d[1:Rank]^(  rrcontrol$Alpha), nrow = Rank) %*%
           t(temp$v[, 1:Rank, drop = FALSE])
      A <- t(A)
      Mmat <- t(C.old)  %*% C.old %*% solve(t(C) %*% C.old)


      Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                     Dmat = Dmat, M = M)
    }

    if (rrcontrol$Uncorrelated.latvar) {
        latvar.mat <- x[, colx2.index, drop = FALSE] %*% C
        var.latvar.mat <- var(latvar.mat)
        UU <- chol(var.latvar.mat)
        Ut <- solve(UU)
        Mmat <- t(UU)
        C <- C %*% Ut
        A <- A %*% t(UU)



      Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                     Dmat = Dmat, M = M)
    }


    if (rrcontrol$Quadratic) {
        Mmat <- diag(Rank)
        for (LV in 1:Rank)
            if (( rrcontrol$Crow1positive[LV] && C[1,LV] < 0) ||
               (!rrcontrol$Crow1positive[LV] && C[1,LV] > 0)) {
                C[,LV] <- -C[,LV]
                A[,LV] <- -A[,LV]
                Mmat[LV,LV] <- -1
            }



      Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                     Dmat = Dmat, M = M)
    }


    list(Amat = A, Cmat = C, Dmat = Dmat)
}





rrr.end.expression <- expression({

  if (exists(".VGAM.etamat", envir = VGAMenv))
    rm(".VGAM.etamat", envir = VGAMenv)


  if (control$Quadratic) {
    if (!length(extra))
      extra <- list()
    extra$Cmat <- Cmat      # Saves the latest iteration 
    extra$Dmat <- Dmat      # Not the latest iteration
    extra$B1   <- B1.save   # Not the latest iteration (not good)
  } else {
    Hlist <- replace.constraints(Hlist.save, Amat, colx2.index)
  }

    X.vlm.save <- if (control$Quadratic) {
        tmp300 <- lm2qrrvlm.model.matrix(x = x, Hlist = Hlist.save,
                                         C = Cmat, control = control)
        latvar.mat <- tmp300$latvar.mat  # Needed at the top of new.s.call

        lm2vlm.model.matrix(tmp300$new.latvar.model.matrix,
                            H.list,
                            xij = control$xij)
    } else {
        lm2vlm.model.matrix(x, Hlist, xij = control$xij)
    }


    fv <- tmp.fitted  # Contains \bI \bnu
    eta <- fv + offset
    if (FALSE && control$Rank == 1) {
      ooo <- order(latvar.mat[, 1])
    }
    mu <- family@linkinv(eta, extra)

    if (any(is.na(mu)))
        warning("there are NAs in mu") 

    deriv.mu <- eval(family@deriv)
    wz <- eval(family@weight)
    if (control$checkwz)
      wz <- checkwz(wz, M = M, trace = trace,
                    wzepsilon = control$wzepsilon)
    U <- vchol(wz, M = M, n = n, silent=!trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
    z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset  # Contains \bI \bnu



})



rrr.derivative.expression <- expression({






    which.optimizer <- if (control$Quadratic && control$FastAlgorithm) {
      "BFGS" 
    } else {
      if (iter <= rrcontrol$Switch.optimizer) "Nelder-Mead" else "BFGS"
    }
    if (trace && control$OptimizeWrtC) {
      cat("\n\n")
      cat("Using", which.optimizer, "\n")
      flush.console()
    } 

    constraints <- replace.constraints(constraints, diag(M),
                                       rrcontrol$colx2.index)
    nice31 <- (!control$eq.tol || control$I.tolerances) &&
              all(trivial.constraints(constraints) == 1)

    theta0 <- c(Cmat)
    assign(".VGAM.dot.counter", 0, envir = VGAMenv)
    if (control$OptimizeWrtC) {
      if (control$Quadratic && control$FastAlgorithm) {
        if (iter == 2) {
          if (exists(".VGAM.etamat", envir = VGAMenv))
            rm(".VGAM.etamat", envir = VGAMenv)
        }
        if (iter > 2 && !quasi.newton$convergence) {
          if (zthere <- exists(".VGAM.z", envir = VGAMenv)) {
              ..VGAM.z <- get(".VGAM.z", envir = VGAMenv)
              ..VGAM.U <- get(".VGAM.U", envir = VGAMenv)
              ..VGAM.beta <- get(".VGAM.beta", envir = VGAMenv)
                }
                if (zthere) {
                    z <- matrix(..VGAM.z, n, M)  # minus any offset
                    U <- matrix(..VGAM.U, M, n)
                }

          }
    
          if (iter == 2 || quasi.newton$convergence) {
              NOS <- ifelse(modelno == 3 || modelno == 5, M/2, M)

              canfitok <-
                (exists("CQO.FastAlgorithm", envir=VGAMenv) &&
                get("CQO.FastAlgorithm", envir = VGAMenv))
              if (!canfitok)
                stop("cannot fit this model using fast algorithm")
              p2star <- if (nice31) 
                        ifelse(control$I.toleran,
                               Rank,
                               Rank+0.5*Rank*(Rank+1)) else
                        (NOS*Rank +
                         Rank*(Rank+1)/2 * ifelse(control$eq.tol, 1,NOS))
              p1star <- if (nice31) p1 *
                        ifelse(modelno == 3 || modelno == 5, 2, 1) else
                        (ncol(X.vlm.save) - p2star)
              X.vlm.1save <- if (p1star > 0)
                             X.vlm.save[,-(1:p2star)] else NULL
              quasi.newton <-
                optim(par = Cmat, fn = callcqof, 
                      gr <- if (control$GradientFunction) calldcqo else NULL,
                      method = which.optimizer,
                      control = list(fnscale = 1,
                                     trace = as.integer(control$trace),
                                     parscale = rep(control$Parscale,
                                                    length.out=length(Cmat)),
                                     maxit = 250),
                      etamat = eta, xmat = x, ymat = y, wvec = w,
                      X.vlm.1save = if (nice31) NULL else X.vlm.1save,
                      modelno = modelno, Control = control,
                      n = n, M = M, p1star = p1star,
                      p2star = p2star, nice31 = nice31)


                if (zthere <- exists(".VGAM.z", envir = VGAMenv)) {
                  ..VGAM.z <- get(".VGAM.z", envir = VGAMenv)
                  ..VGAM.U <- get(".VGAM.U", envir = VGAMenv)
                  ..VGAM.beta <- get(".VGAM.beta", envir = VGAMenv)
                }
                if (zthere) {
                  z <- matrix(..VGAM.z, n, M)  # minus any offset
                  U <- matrix(..VGAM.U, M, n)
                }
          } else {
            if (exists(".VGAM.offset", envir = VGAMenv))
              rm(".VGAM.offset", envir = VGAMenv)
          }
        } else {
          use.reltol <- if (length(rrcontrol$Reltol) >= iter) 
              rrcontrol$Reltol[iter] else rev(rrcontrol$Reltol)[1]
          quasi.newton <-
            optim(par = theta0,
                  fn = rrr.derivC.ResSS, 
                  method = which.optimizer,
                  control = list(fnscale = rrcontrol$Fnscale, 
                                 maxit = rrcontrol$Maxit,
                                 abstol = rrcontrol$Abstol,
                                 reltol = use.reltol),
                  U = U, z = if (control$I.tolerances) z + offset else z,
                  M = M, xmat = x,  # varbix2 = varbix2,
                  Hlist = Hlist, rrcontrol = rrcontrol)
        }




      Cmat <- matrix(quasi.newton$par, p2, Rank, byrow = FALSE)

      if (Rank > 1 && rrcontrol$I.tolerances) {
        numat <- x[, rrcontrol$colx2.index, drop = FALSE] %*% Cmat
        evnu <- eigen(var(numat))
        Cmat <- Cmat %*% evnu$vector
        numat <- x[, rrcontrol$colx2.index, drop = FALSE] %*% Cmat
        offset <- if (Rank > 1) -0.5*rowSums(numat^2) else -0.5*numat^2
      }
    }


    alt <- valt.1iter(x = x, z = z, U = U, Hlist = Hlist,
                     C = Cmat, nice31 = nice31,
                     control = rrcontrol,
                     lp.names = predictors.names)


    if (length(alt$offset))
        offset <- alt$offset

    B1.save <- alt$B1 # Put later into extra  
    tmp.fitted <- alt$fitted  # contains \bI_{Rank} \bnu if Corner

    if (modelno != 33 && control$OptimizeWrtC)
        alt <- rrr.normalize(rrc = rrcontrol, A = alt$Amat, C = alt$Cmat,
                             x = x, Dmat = alt$Dmat)

    if (trace && control$OptimizeWrtC) {
      cat("\n")
      cat(which.optimizer, "using optim():\n")
      cat("Objective  = ", quasi.newton$value, "\n")
      cat("Parameters (= c(C)) = ", if (length(quasi.newton$par) < 5)
          "" else "\n")
      cat(alt$Cmat, fill = TRUE)
      cat("\n")
      cat("Number of function evaluations  = ", quasi.newton$count[1], "\n")
      if (length(quasi.newton$message))
        cat("Message  = ", quasi.newton$message, "\n")
      cat("\n")
      flush.console()
    }



    Amat <- alt$Amat  # Needed in rrr.end.expression 
    Cmat <- alt$Cmat  # Needed in rrr.end.expression if Quadratic 
    Dmat <- alt$Dmat  # Put later into extra  

    eval(rrr.end.expression)    # Put Amat into Hlist, and create new z
})



rrr.derivC.ResSS <- function(theta, U, z, M, xmat, Hlist, rrcontrol,
                          omit.these = NULL) {

  if (rrcontrol$trace) {
      cat(".")
      flush.console()
  }
  alreadyThere <- exists(".VGAM.dot.counter", envir = VGAMenv)
  if (alreadyThere) {
    VGAM.dot.counter <- get(".VGAM.dot.counter", envir = VGAMenv)
    VGAM.dot.counter <- VGAM.dot.counter + 1 
    assign(".VGAM.dot.counter", VGAM.dot.counter, envir = VGAMenv)
    if (VGAM.dot.counter > max(50, options()$width - 5)) {
      if (rrcontrol$trace) {
          cat("\n")
          flush.console()
      }
      assign(".VGAM.dot.counter", 0, envir = VGAMenv)
    }
  }

  Cmat <- matrix(theta, length(rrcontrol$colx2.index), rrcontrol$Rank)


    tmp700 <-
      lm2qrrvlm.model.matrix(x = xmat, Hlist = Hlist,
                             no.thrills = !rrcontrol$Corner,
                             C = Cmat, control = rrcontrol, assign = FALSE)
    Hlist <- tmp700$constraints  # Does not contain \bI_{Rank} \bnu

    if (rrcontrol$Corner) {
      z <- as.matrix(z)  # should actually call this zedd
      z[, rrcontrol$Index.corner] <-
      z[, rrcontrol$Index.corner] - tmp700$latvar.mat
    }

    if (length(tmp700$offset)) z <- z - tmp700$offset


    vlm.wfit(xmat = tmp700$new.latvar.model.matrix, zmat = z,
             Hlist = Hlist, ncolx = ncol(xmat), U = U, only.ResSS = TRUE,
             matrix.out = FALSE, is.vlmX = FALSE, ResSS = TRUE,
             qr = FALSE, Eta.range = rrcontrol$Eta.range,
             xij = rrcontrol$xij)$ResSS
}




rrvglm.optim.control <- function(Fnscale = 1,
                                 Maxit = 100, 
                                 Switch.optimizer = 3,
                                 Abstol = -Inf, 
                                 Reltol = sqrt(.Machine$double.eps),
                                 ...) {




    list(Fnscale = Fnscale, 
         Maxit = Maxit,
         Switch.optimizer = Switch.optimizer,
         Abstol = Abstol,
         Reltol = Reltol)
}



nlminbcontrol <- function(Abs.tol = 10^(-6),
                          Eval.max = 91,
                          Iter.max = 91,
                          Rel.err = 10^(-6),
                          Rel.tol = 10^(-6),
                          Step.min = 10^(-6),
                          X.tol = 10^(-6),
                          ...) {


  list(Abs.tol = Abs.tol,
       Eval.max = Eval.max,
       Iter.max = Iter.max,
       Rel.err = Rel.err,
       Rel.tol = Rel.tol,
       Step.min = Step.min,
       X.tol = X.tol)
}




Coef.qrrvglm <-
  function(object,
           varI.latvar = FALSE,
           refResponse = NULL, ...) {




  if (length(varI.latvar) != 1 || !is.logical(varI.latvar)) 
    stop("'varI.latvar' must be TRUE or FALSE")
  if (length(refResponse) > 1)
    stop("argument 'refResponse' must be of length 0 or 1")
  if (length(refResponse) &&
      is.Numeric(refResponse))
      if (!is.Numeric(refResponse, length.arg = 1,
                      integer.valued = TRUE))
        stop("bad input for argument 'refResponse'")
  if (!is.logical(ConstrainedQO <- object@control$ConstrainedQO))
    stop("cannot determine whether the model is constrained or not")

  ocontrol <- object@control
  coef.object <- object@coefficients 
  Rank <- ocontrol$Rank 
  M <- object@misc$M
  NOS <- if (length(object@y)) ncol(object@y) else M
  MSratio <- M / NOS  # First value is g(mean) = quadratic form in latvar
  Quadratic <- if (ConstrainedQO) ocontrol$Quadratic else TRUE
  if (!Quadratic) stop("object is not a quadratic ordination object")
  p1 <- length(ocontrol$colx1.index)
  p2 <- length(ocontrol$colx2.index)
  Index.corner <- ocontrol$Index.corner
  str0 <- ocontrol$str0
  eq.tolerances <- ocontrol$eq.tolerances
  Dzero <- ocontrol$Dzero
  Corner <- if (ConstrainedQO) ocontrol$Corner else FALSE

  estI.tol <- if (ConstrainedQO) object@control$I.tolerances else FALSE
  modelno <- object@control$modelno  # 1, 2, 3, 4, 5, 6, 7 or 0
  combine2 <- c(str0, if (Corner) Index.corner else NULL)
  NoA <- length(combine2) == M  # A is fully known.

  Qoffset <- if (Quadratic) ifelse(estI.tol, 0, sum(1:Rank)) else 0

  ynames <- object@misc$ynames
  if (!length(ynames)) ynames <- object@misc$predictors.names
  if (!length(ynames)) ynames <- object@misc$ynames
  if (!length(ynames)) ynames <- paste("Y", 1:NOS, sep = "")
  lp.names <- object@misc$predictors.names
  if (!length(lp.names)) lp.names <- NULL 

  dzero.vector <- rep(FALSE, length = M)
  if (length(Dzero))
    dzero.vector[Dzero] <- TRUE
  names(dzero.vector) <- ynames 
  latvar.names <- if (Rank == 1)
    "latvar" else
    paste("latvar", 1:Rank, sep = "")



  td.expression <- function(Dmat, Amat, M, Dzero, Rank, bellshaped) {


    Tolerance <- Darray <- m2a(Dmat, M = Rank)
    for (ii in 1:M)
      if (length(Dzero) && any(Dzero == ii)) {
        Tolerance[, , ii] <- NA   # Darray[,,ii] == O 
        bellshaped[ii] <- FALSE 
      } else {
        Tolerance[, , ii] <- -0.5 * solve(Darray[, , ii])
        bellshaped[ii] <- all(eigen(Tolerance[, , ii])$values > 0)
      }
    optimum <- matrix(NA_real_, Rank, M)
    for (ii in 1:M)
      if (bellshaped[ii])
        optimum[, ii] <- Tolerance[, , ii] %*% cbind(Amat[ii, ])

    list(optimum    = optimum,
         Tolerance  = Tolerance,
         Darray     = Darray,
         bellshaped = bellshaped)
  }





  Amat <- object@extra$Amat  # M  x Rank
  Cmat <- object@extra$Cmat  # p2 x Rank
  Dmat <- object@extra$Dmat  #
  B1   <- object@extra$B1    #
  bellshaped <- rep(FALSE, length = M)

  if (is.character(refResponse)) {
      refResponse <- (1:NOS)[refResponse == ynames]
      if (length(refResponse) != 1)
         stop("could not match argument 'refResponse' with any response")
  }
  ptr1 <- 1
  candidates <- if (length(refResponse)) refResponse else {
      if (length(ocontrol$Dzero)) (1:M)[-ocontrol$Dzero] else (1:M)}
  repeat {
    if (ptr1 > 0) {
      this.spp <- candidates[ptr1]
    }
  elts <- Dmat[this.spp,, drop = FALSE]
      if (length(elts) < Rank)
        elts <- matrix(elts, 1, Rank)
      Dk <- m2a(elts, M = Rank)[, , 1]  # Hopefully negative-def 
      temp400 <- eigen(Dk)
      ptr1 <- ptr1 + 1 
      if (all(temp400$value < 0))
        break
      if (ptr1 > length(candidates))
        break
  }
  if (all(temp400$value < 0)) {
    temp1tol <- -0.5 * solve(Dk)
    dim(temp1tol) <- c(Rank,Rank)
    Mmat <- t(chol(temp1tol))
    if (ConstrainedQO) {
      temp900 <- solve(t(Mmat))
      Cmat <- Cmat %*% temp900
      Amat <- Amat %*% Mmat
    }
    if (length(Cmat)) {
      temp800 <- crow1C(Cmat, ocontrol$Crow1positive, amat = Amat)
      Cmat <- temp800$cmat
      Amat <- temp800$amat
    }



    Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                   Dmat = Dmat, M = M)



    retlist <- td.expression(Dmat = Dmat, Amat = Amat, M = M,
                             Dzero = Dzero, Rank = Rank,
                             bellshaped = bellshaped)
    optimum    <- retlist$optimum
    Tolerance  <- retlist$Tolerance
    Darray     <- retlist$Darray
    bellshaped <- retlist$bellshaped




    } else {
      if (length(refResponse) == 1) 
        stop("tolerance matrix specified by 'refResponse' ",
             "is not positive-definite") else
        warning("could not find any positive-definite ",
                "tolerance matrix")
    }


  if (ConstrainedQO)
    if (Rank > 1) {
      if (!length(xmat <- object@x))
        stop("cannot obtain the model matrix")
      numat <- xmat[,ocontrol$colx2.index, drop = FALSE] %*% Cmat
      evnu <- eigen(var(numat))
      Mmat <- solve(t(evnu$vector))
      Cmat <- Cmat %*% evnu$vector  # == Cmat %*% solve(t(Mmat))
      Amat <- Amat %*% Mmat
      temp800 <- crow1C(Cmat, ocontrol$Crow1positive, amat = Amat)
      Cmat <- temp800$cmat
      Amat <- temp800$amat


      Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                     Dmat = Dmat, M = M)



      retlist <- td.expression(Dmat = Dmat, Amat = Amat, M = M,
                               Dzero = Dzero, Rank = Rank,
                               bellshaped = bellshaped)
      optimum    <- retlist$optimum
      Tolerance  <- retlist$Tolerance
      Darray     <- retlist$Darray
      bellshaped <- retlist$bellshaped
  }


  if (ConstrainedQO)
    if (varI.latvar) {
      if (!length(xmat <- object@x))
        stop("cannot obtain the model matrix")
      numat <- xmat[,ocontrol$colx2.index, drop = FALSE] %*% Cmat
      sdnumat <- apply(cbind(numat), 2, sd)
      Mmat <- if (Rank > 1) diag(sdnumat) else matrix(sdnumat, 1, 1)
      Cmat <- Cmat %*% solve(t(Mmat))
      Amat <- Amat %*% Mmat
      temp800 <- crow1C(Cmat, ocontrol$Crow1positive, amat = Amat)
      Cmat <- temp800$cmat
      Amat <- temp800$amat
               Cmat # Not needed



      Dmat <- adjust.Dmat.expression(Mmat = Mmat, Rank = Rank,
                                     Dmat = Dmat, M = M)



      retlist <- td.expression(Dmat = Dmat, Amat = Amat, M = M,
                               Dzero = Dzero, Rank = Rank,
                               bellshaped = bellshaped)
      optimum    <- retlist$optimum
      Tolerance  <- retlist$Tolerance
      Darray     <- retlist$Darray
      bellshaped <- retlist$bellshaped
    }


  cx1i <- ocontrol$colx1.index
  maximum <- if (length(cx1i) == 1 && names(cx1i) == "(Intercept)") {
      eta.temp <- B1
      for (ii in 1:M)
        eta.temp[ii] <- eta.temp[ii] + 
            Amat[ii, , drop = FALSE] %*% optimum[, ii, drop = FALSE] +
            t(optimum[, ii, drop = FALSE]) %*%
            Darray[,, ii, drop = TRUE] %*% optimum[, ii, drop = FALSE]
      mymax <- object@family@linkinv(rbind(eta.temp), extra = object@extra)  
      c(mymax)  # Convert from matrix to vector 
    } else {
      5 * rep(NA_real_, length.out = M)  # Make "numeric"
  }
  names(maximum) <- ynames
    
  latvar.mat <- if (ConstrainedQO) {
    object@x[, ocontrol$colx2.index, drop = FALSE] %*% Cmat 
  } else {
    object@latvar
  }

  dimnames(Amat) <- list(lp.names, latvar.names)
  if (ConstrainedQO)
    dimnames(Cmat) <- list(names(ocontrol$colx2.index), latvar.names)
  if (!length(xmat <- object@x)) stop("cannot obtain the model matrix")
  dimnames(latvar.mat) <- list(dimnames(xmat)[[1]], latvar.names)

  ans <- 
  new(Class <- if (ConstrainedQO) "Coef.qrrvglm" else "Coef.uqo",
       A = Amat, B1 = B1, Constrained = ConstrainedQO, D = Darray,
       NOS = NOS, Rank = Rank,
       latvar = latvar.mat,
       latvar.order = latvar.mat,
       Optimum = optimum, 
       Optimum.order = optimum, 
       bellshaped = bellshaped,
       Dzero = dzero.vector,
       Maximum = maximum,
       Tolerance = Tolerance)
  if (ConstrainedQO) {ans@C <- Cmat} else {Cmat <- NULL}

  for (rrr in 1:Rank)
    ans@Optimum.order[rrr, ] <- order(ans@Optimum[rrr, ])
  for (rrr in 1:Rank)
    ans@latvar.order[, rrr] <- order(ans@latvar[, rrr])

  if (length(object@misc$estimated.dispersion) &&
      object@misc$estimated.dispersion) {
    p <- length(object@coefficients)
    n <- object@misc$n
    M <- object@misc$M
    NOS <- if (length(object@y)) ncol(object@y) else M
    pstar <- if (ConstrainedQO) (p + length(Cmat)) else
             p + n*Rank # Adjustment; not sure about UQO 
    adjusted.dispersion <- object@misc$dispersion * (n*M - p) /
            (n*M - pstar)
    ans@dispersion <- adjusted.dispersion 
  }

  if (MSratio > 1) {
    keepIndex <- seq(from = 1, to = M, by = MSratio)
    ans@Dzero <- ans@Dzero[keepIndex]
    ans@Optimum <- ans@Optimum[,keepIndex, drop = FALSE]
    ans@Tolerance <- ans@Tolerance[,,keepIndex, drop = FALSE]
    ans@bellshaped <- ans@bellshaped[keepIndex]
    names(ans@Dzero) <- ynames
  } else {
    dimnames(ans@D) <- list(latvar.names, latvar.names, ynames)
  }
  names(ans@bellshaped) <- ynames 
  dimnames(ans@Optimum) <- list(latvar.names, ynames)
  dimnames(ans@Tolerance) <- list(latvar.names, latvar.names, ynames)
  ans 
}  # End of Coef.qrrvglm


setClass(Class = "Coef.rrvglm", representation(
      "A"             = "matrix",
      "B1"            = "matrix",  # This may be unassigned if p1 = 0.
      "C"             = "matrix",
      "Rank"          = "numeric",
      "colx1.index"   = "numeric",
      "colx2.index"   = "numeric",
      "Atilde"        = "matrix"))


setClass(Class = "Coef.uqo", representation(
      "A"             = "matrix",
      "B1"            = "matrix",
      "Constrained"   = "logical",
      "D"             = "array",
      "NOS"           = "numeric",
      "Rank"          = "numeric",
      "latvar"        = "matrix",
      "latvar.order"  = "matrix",
      "Maximum"       = "numeric",
      "Optimum"       = "matrix",
      "Optimum.order" = "matrix",
      "bellshaped"    = "logical",
      "dispersion"    = "numeric",
      "Dzero"         = "logical",
      "Tolerance"     = "array"))


setClass(Class = "Coef.qrrvglm", representation(
      "C"            = "matrix"),
    contains = "Coef.uqo")




show.Coef.qrrvglm <- function(x, ...) {

  object <- x 
  Rank <- object@Rank
  M <- nrow(object@A)
  NOS <- object@NOS
  mymat <- matrix(NA_real_, NOS, Rank)
  if (Rank == 1) {  # || object@Diagonal
    for (ii in 1:NOS) {
      fred <- if (Rank > 1)
                diag(object@Tolerance[, , ii, drop = FALSE]) else
                object@Tolerance[, , ii]
      if (all(fred > 0))
        mymat[ii,] <- sqrt(fred)
    }
    dimnames(mymat) <- list(dimnames(object@Tolerance)[[3]],
                            if (Rank == 1) "latvar" else
                            paste("Tolerance", dimnames(mymat)[[2]],
                                  sep = ""))
    } else {
      for (ii in 1:NOS) {
        fred <- eigen(object@Tolerance[, , ii])
          if (all(fred$value > 0))
              mymat[ii, ] <- sqrt(fred$value)
      }
      dimnames(mymat) <- list(dimnames(object@Tolerance)[[3]],
                              paste("tol", 1:Rank, sep = ""))
    }

    dimnames(object@A) <- list(dimnames(object@A)[[1]],
      if (Rank > 1) paste("A", dimnames(object@A)[[2]], sep = ".") else
                          "A")

    Maximum <- if (length(object@Maximum))
               cbind(Maximum = object@Maximum) else NULL
    if (length(Maximum) && length(mymat) && Rank == 1)
      Maximum[is.na(mymat),] <- NA

   optmat <- cbind(t(object@Optimum))
    dimnames(optmat) <- list(dimnames(optmat)[[1]],
        if (Rank > 1)
          paste("Optimum", dimnames(optmat)[[2]], sep = ".") else
          "Optimum")
    if (length(optmat) && length(mymat) && Rank == 1)
        optmat[is.na(mymat), ] <- NA

    if ( object@Constrained ) {
      cat("\nC matrix (constrained/canonical coefficients)\n")
      print(object@C, ...)
    }
    cat("\nB1 and A matrices\n")
    print(cbind(t(object@B1),
                A = object@A), ...)
    cat("\nOptimums and maximums\n")
    print(cbind(Optimum = optmat,
                Maximum), ...)
    if (Rank > 1) {  # !object@Diagonal && Rank > 1
      cat("\nTolerances\n")
    } else {
      cat("\nTolerance\n")
    }
    print(mymat, ...)

    cat("\nStandard deviation of the latent variables (site scores)\n")
    print(apply(cbind(object@latvar), 2, sd))
    invisible(object)
}





setMethod("show", "Coef.qrrvglm", function(object)
    show.Coef.qrrvglm(object))








setMethod("summary", "qrrvglm", function(object, ...)
    summary.qrrvglm(object, ...))



predictqrrvglm <-
  function(object,
           newdata = NULL,
           type = c("link", "response", "latvar", "terms"),
           se.fit = FALSE,
           deriv = 0,
           dispersion = NULL,
           extra = object@extra, 
           varI.latvar = FALSE, refResponse = NULL, ...) {
  if (se.fit)
    stop("cannot handle se.fit == TRUE yet")
  if (deriv != 0)
    stop("derivative is not equal to 0")

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("link", "response", "latvar", "terms"))[1]
  if (type == "latvar")
    stop("cannot handle type='latvar' yet")
  if (type == "terms")
    stop("cannot handle type='terms' yet")

  M <- object@misc$M
  Rank  <- object@control$Rank

  na.act <- object@na.action
  object@na.action <- list()

  if (!length(newdata) &&
      type == "response" &&
      length(object@fitted.values)) {
    if (length(na.act)) {
      return(napredict(na.act[[1]], object@fitted.values))
    } else {
      return(object@fitted.values)
    }
  }

  if (!length(newdata)) {
    X <- model.matrixvlm(object, type = "lm", ...)
    offset <- object@offset
    tt <- object@terms$terms   # terms(object)
    if (!length(object@x))
      attr(X, "assign") <- attrassignlm(X, tt)
  } else {
    if (is.smart(object) && length(object@smart.prediction)) {
      setup.smart("read", smart.prediction = object@smart.prediction)
    }

    tt <- object@terms$terms  # terms(object)  # 20030811; object@terms$terms
    X <- model.matrix(delete.response(tt), newdata,
                      contrasts = if (length(object@contrasts))
                                  object@contrasts else NULL,
                      xlev = object@xlevels)

    if (nrow(X) != nrow(newdata)) {
      as.save <- attr(X, "assign")
      X <- X[rep(1, nrow(newdata)),, drop = FALSE]
      dimnames(X) <- list(dimnames(newdata)[[1]], "(Intercept)")
      attr(X, "assign") <- as.save  # Restored
    }

    offset <- if (!is.null(off.num<-attr(tt,"offset"))) {
      eval(attr(tt,"variables")[[off.num+1]], newdata)
    } else if (!is.null(object@offset))
      eval(object@call$offset, newdata)

      if (any(c(offset) != 0))
        stop("currently cannot handle nonzero offsets")

      if (is.smart(object) && length(object@smart.prediction)) {
        wrapup.smart()
    }

    attr(X, "assign") <- attrassigndefault(X, tt)
  }

  ocontrol <- object@control

    Rank <- ocontrol$Rank
    NOS <- ncol(object@y)
    sppnames <- dimnames(object@y)[[2]]
    modelno <- ocontrol$modelno  # 1, 2, 3, 5 or 0
    M <- if (any(slotNames(object) == "predictors") &&
             is.matrix(object@predictors))
           ncol(object@predictors) else
           object@misc$M
    MSratio <- M / NOS  # First value is g(mean) = quadratic form in latvar
    if (MSratio != 1) stop("can only handle MSratio == 1 for now")


    if (length(newdata)) {
      Coefs <- Coef(object, varI.latvar = varI.latvar, refResponse = refResponse)
      X1mat <- X[, ocontrol$colx1.index, drop = FALSE]
      X2mat <- X[, ocontrol$colx2.index, drop = FALSE]
      latvarmat <- as.matrix(X2mat %*% Coefs@C)  # n x Rank

      etamat <- as.matrix(X1mat %*% Coefs@B1 + latvarmat %*% t(Coefs@A))
      which.species <- 1:NOS  # Do it all for all species
      for (sppno in 1:length(which.species)) {
        thisSpecies <- which.species[sppno]
        Dmat <- matrix(Coefs@D[,,thisSpecies], Rank, Rank)
        etamat[, thisSpecies] <- etamat[, thisSpecies] +
                                 mux34(latvarmat, Dmat, symmetric = TRUE)
      }
    } else {
      etamat <-  object@predictors
  }

  pred <-
    switch(type,
           response = {
           fv <- if (length(newdata))
                  object@family@linkinv(etamat, extra) else
                  fitted(object)
            if (M > 1 && is.matrix(fv)) {
      dimnames(fv) <- list(dimnames(fv)[[1]],
                           dimnames(object@fitted.values)[[2]])
    }
    fv
  },
           link = etamat,
           latvar = stop("failure here"),
           terms  = stop("failure here"))

  if (!length(newdata) && length(na.act)) {
    if (se.fit) {
      pred$fitted.values <- napredict(na.act[[1]], pred$fitted.values)
      pred$se.fit <- napredict(na.act[[1]], pred$se.fit)
    } else {
        pred <- napredict(na.act[[1]], pred)
    }
  }
  pred
}


setMethod("predict", "qrrvglm", function(object, ...)
  predictqrrvglm(object, ...))


coefqrrvglm <- function(object, matrix.out = FALSE,
                       label = TRUE) {
  if (matrix.out)
    stop("currently cannot handle matrix.out = TRUE")
  coefvlm(object, matrix.out = matrix.out, label = label)
}



residualsqrrvglm  <-
  function(object,
           type = c("deviance", "pearson", "working", "response", "ldot"),
           matrix.arg = TRUE) {
  stop("this function has not been written yet")
}


setMethod("residuals",  "qrrvglm",
  function(object, ...)
    residualsqrrvglm(object, ...))





show.rrvglm <- function(x, ...) {
  if (!is.null(cl <- x@call)) {
    cat("Call:\n")
    dput(cl)
  }
  vecOfBetas <- x@coefficients
  if (any(nas <- is.na(vecOfBetas))) {
    if (is.null(names(vecOfBetas)))
      names(vecOfBetas) <- paste("b",
            1:length(vecOfBetas), sep = "")
    cat("\nCoefficients: (", sum(nas),
        " not defined because of singularities)\n", sep = "")
  } else 
      cat("\nCoefficients:\n")
  print.default(vecOfBetas, ...)    # used to be print()

  if (FALSE) {
    Rank <- x@Rank
    if (!length(Rank))
      Rank <- sum(!nas)
  }

  if (FALSE) {
    nobs <- if (length(x@df.total)) x@df.total else length(x@residuals)
    rdf <- x@df.residual
    if (!length(rdf))
      rdf <- nobs - Rank
  }
  cat("\n")

  if (length(deviance(x)))
    cat("Residual deviance:", format(deviance(x)), "\n")
  if (length(vll <- logLik.vlm(x)))
    cat("Log-likelihood:", format(vll), "\n")

  if (length(x@criterion)) {
    ncrit <- names(x@criterion)
    for (iii in ncrit)
      if (iii != "loglikelihood" &&
          iii != "deviance")
        cat(paste(iii, ":", sep = ""),
            format(x@criterion[[iii]]), "\n")
  }

  invisible(x)
}






setMethod("show", "rrvglm", function(object) show.rrvglm(object))








summary.rrvglm <- function(object, correlation = FALSE,
                           dispersion = NULL, digits = NULL, 
                           numerical = TRUE,
                           h.step = 0.0001, 
                           kill.all = FALSE, omit13 = FALSE,
                           fixA = FALSE, 
                           presid = TRUE, 
                           nopredictors = FALSE, ...) {









    if (!is.Numeric(h.step, length.arg = 1) ||
        abs(h.step) > 1)
      stop("bad input for 'h.step'")

    if (!object@control$Corner)
      stop("this function works with corner constraints only")

    if (is.null(dispersion))
      dispersion <- object@misc$dispersion

    newobject <- as(object, "vglm")


    stuff <- summaryvglm(newobject,
                         correlation = correlation,
                         dispersion = dispersion,
                         presid = presid)

    answer <-
    new(Class = "summary.rrvglm",
        object,
        call = stuff@call,
        coef3 = stuff@coef3,
        cov.unscaled = stuff@cov.unscaled,
        correlation = stuff@correlation,
        df = stuff@df,
        sigma = stuff@sigma)


    if (is.numeric(stuff@dispersion))
      slot(answer, "dispersion") <- stuff@dispersion

    if (presid && length(stuff@pearson.resid))
      slot(answer, "pearson.resid") <- stuff@pearson.resid



    tmp5 <- get.rrvglm.se1(object, omit13 = omit13,
                           numerical = numerical, h.step = h.step,
                           kill.all = kill.all, fixA = fixA, ...) 
    if (any(diag(tmp5$cov.unscaled) <= 0) ||
       any(eigen(tmp5$cov.unscaled)$value <= 0)) {
        warning("cov.unscaled is not positive definite") 
    }

    answer@cov.unscaled <- tmp5$cov.unscaled 

    od <- if (is.numeric(object@misc$disper))
        object@misc$disper else
        object@misc$default.disper
    if (is.numeric(dispersion)) {
      if (is.numeric(od) && dispersion != od)
          warning("dispersion != object@misc$dispersion; ",
                  "using the former")
    } else {
      dispersion <- if (is.numeric(od)) od else 1
    }

    tmp8 <- object@misc$M - object@control$Rank - 
            length(object@control$str0)
    answer@df[1] <- answer@df[1] + tmp8 * object@control$Rank
    answer@df[2] <- answer@df[2] - tmp8 * object@control$Rank
    if (dispersion == 0) {
      dispersion <- tmp5$ResSS / answer@df[2]  # Estimate 
    }

    answer@coef3 <- get.rrvglm.se2(answer@cov.unscaled,
                                   dispersion = dispersion,
                                   coefficients = tmp5$coefficients)

    answer@dispersion <- dispersion
    answer@sigma <- dispersion^0.5


    answer@misc$nopredictors <- nopredictors  # 20150925

    answer
}






get.rrvglm.se1 <- function(fit, omit13 = FALSE, kill.all = FALSE,
                           numerical = TRUE,
                           fixA = FALSE, h.step = 0.0001,
                           trace.arg = FALSE, ...) {




  if (length(fit@control$Nested) && fit@control$Nested)
    stop("sorry, cannot handle nested models yet")

  str0 <- fit@control$str0


  if (!length(fit@x))
    stop("fix@x is empty. Run rrvglm(... , x = TRUE)")

  colx1.index <- fit@control$colx1.index  # May be NULL
  colx2.index <- fit@control$colx2.index 
  Hlist <- fit@constraints
  ncolHlist <- unlist(lapply(Hlist, ncol))

  p1 <- length(colx1.index)  # May be 0
  p2 <- length(colx2.index)

  Rank <- fit@control$Rank  # fit@misc$Nested.Rank   

  Amat <- fit@constraints[[colx2.index[1]]]
  B1mat <- if (p1)
    coefvlm(fit, matrix.out = TRUE)[colx1.index, , drop = FALSE] else
    NULL
  C.try <- coefvlm(fit, matrix.out= TRUE)[colx2.index, , drop = FALSE]
  Cmat <- C.try %*% Amat %*% solve(t(Amat) %*% Amat)

  x1mat <- if (p1) fit@x[, colx1.index, drop = FALSE] else NULL
  x2mat <- fit@x[, colx2.index, drop = FALSE]

  wz <- weights(fit, type = "work")  # old: wweights(fit)  #fit@weights
  if (!length(wz))
    stop("cannot get fit@weights")

  M <- fit@misc$M
  n <- fit@misc$n
  Index.corner <- fit@control$Index.corner   # used to be (1:Rank);
  zmat <- fit@predictors + fit@residuals
  theta <- c(Amat[-c(Index.corner,str0), ])
  if (fit@control$checkwz)
    wz <- checkwz(wz, M = M, trace = trace,
                  wzepsilon = fit@control$wzepsilon)
   U <- vchol(wz, M = M, n = n, silent= TRUE)

  delct.da <- if (numerical) {
    num.deriv.rrr(fit, M = M, r = Rank,
                  x1mat = x1mat, x2mat = x2mat, p2 = p2, 
                  Index.corner, Aimat = Amat,
                  B1mat = B1mat, Cimat = Cmat,
                  h.step = h.step,
                  colx2.index = colx2.index,
                  xij = fit@control$xij,
                  str0 = str0)
  } else {
    dctda.fast.only(theta = theta, wz = wz,
                    U = U, zmat,
                    M = M, r = Rank, x1mat = x1mat,
                    x2mat = x2mat, p2 = p2,
                    Index.corner, Aimat = Amat,
                    B1mat = B1mat, Cimat = Cmat,
                    xij = fit@control$xij,
                    str0 = str0)
  }




  newobject <- as(fit, "vglm")




  sfit2233 <- summaryvglm(newobject) 
  d8 <-  dimnames(sfit2233@cov.unscaled)[[1]]
  cov2233 <- solve(sfit2233@cov.unscaled)  # Includes any intercepts
  dimnames(cov2233) <- list(d8, d8)

  log.vec33 <- NULL 
  nassign <- names(fit@constraints) 
  choose.from <-  varassign(fit@constraints, nassign)
  for (ii in nassign)
    if (any(ii == names(colx2.index))) {
      log.vec33 <- c(log.vec33, choose.from[[ii]])
    }
    cov33 <- cov2233[ log.vec33, log.vec33, drop = FALSE]  # r*p2 by r*p2
    cov23 <- cov2233[-log.vec33, log.vec33, drop = FALSE]
    cov22 <- cov2233[-log.vec33,-log.vec33, drop = FALSE]


    latvar.mat <- x2mat %*% Cmat
    offs <- matrix(0, n, M)  # The "0" handles str0's 
    offs[, Index.corner] <- latvar.mat
    if (M == (Rank + length(str0)))
      stop("cannot handle full-rank models yet")
    cm <- matrix(0, M, M - Rank - length(str0))
    cm[-c(Index.corner, str0), ] <- diag(M - Rank - length(str0))

    Hlist <- vector("list", length(colx1.index)+1) 
    names(Hlist) <- c(names(colx1.index), "I(latvar.mat)")
    for (ii in names(colx1.index))
      Hlist[[ii]] <- fit@constraints[[ii]]
    Hlist[["I(latvar.mat)"]] <- cm


    if (p1) {
      ooo <- fit@assign
      bb <- NULL 
      for (ii in 1:length(ooo)) {
        if (any(ooo[[ii]][1] == colx1.index))
          bb <- c(bb, names(ooo)[ii])
      }

      has.intercept <- any(bb == "(Intercept)")
      bb[bb == "(Intercept)"] <- "1"
      if (p1 > 1)
        bb <- paste(bb, collapse = "+")
      if (has.intercept) {
        bb <- paste("zmat - offs ~ ", bb, " + I(latvar.mat)", collapse = " ")
      } else {
        bb <- paste("zmat - offs ~ -1 + ", bb, " + I(latvar.mat)", collapse = " ")
      }
      bb <- as.formula(bb)
    } else {
      bb <- as.formula("zmat - offs ~ -1 + I(latvar.mat)")
    }


    if (fit@misc$dataname == "list") {
      dspec <- FALSE
    } else {
      mytext1 <- "exists(x = fit@misc$dataname, envir = VGAMenv)"
      myexp1 <- parse(text = mytext1)
      is.there <- eval(myexp1)
      bbdata <- if (is.there)
                get(fit@misc$dataname, envir = VGAMenv) else
                get(fit@misc$dataname)
      dspec <- TRUE
    }


    fit1122 <- if (dspec)
               vlm(bb,
                   constraints = Hlist, criterion = "d", weights = wz,
                   data = bbdata,
                   save.weights = TRUE, smart = FALSE, trace = trace.arg,
                   x.arg = TRUE) else
               vlm(bb,
                   constraints = Hlist, criterion = "d", weights = wz,
                   save.weights = TRUE, smart = FALSE, trace = trace.arg,
                   x.arg = TRUE)



    sfit1122 <- summaryvlm(fit1122)
    d8 <-  dimnames(sfit1122@cov.unscaled)[[1]]
    cov1122 <- solve(sfit1122@cov.unscaled)
    dimnames(cov1122) <- list(d8, d8)

    lcs <- length(coefvlm(sfit1122))
    log.vec11 <- (lcs-(M-Rank-length(str0))*Rank+1):lcs
    cov11 <- cov1122[log.vec11,  log.vec11, drop = FALSE]
    cov12 <- cov1122[ log.vec11, -log.vec11, drop = FALSE]
    cov22 <- cov1122[-log.vec11, -log.vec11, drop = FALSE]
    cov13 <- delct.da %*% cov33


    if (omit13) 
      cov13 <- cov13 * 0   # zero it

    if (kill.all) {
      cov13 <- cov13 * 0   # zero it
      if (fixA) {
        cov12 <- cov12 * 0   # zero it
      } else {
        cov23 <- cov23 * 0   # zero it
      }
    }

   cov13 <- -cov13  # Richards (1961)

    if (fixA) {
      cov.unscaled <- rbind(cbind(cov1122, rbind(cov13, cov23)),
                            cbind(t(cov13), t(cov23), cov33))
    } else {
      cov.unscaled <- rbind(cbind(cov11, cov12, cov13),
                            cbind(rbind(t(cov12), t(cov13)), cov2233))
    }

    ans <- solve(cov.unscaled)

    acoefs <- c(fit1122@coefficients[log.vec11], fit@coefficients)
    dimnames(ans) <- list(names(acoefs), names(acoefs))
    list(cov.unscaled = ans,
         coefficients = acoefs,
         ResSS       = sfit1122@ResSS)
}



get.rrvglm.se2 <- function(cov.unscaled, dispersion = 1, coefficients) {

  d8 <-  dimnames(cov.unscaled)[[1]]
  ans <- matrix(coefficients, length(coefficients), 3) 
  ans[, 2] <- sqrt(dispersion) * sqrt(diag(cov.unscaled))
  ans[, 3] <- ans[, 1] / ans[, 2]
  dimnames(ans) <- list(d8, c("Estimate", "Std. Error", "z value"))
  ans
}



num.deriv.rrr <- function(fit, M, r, x1mat, x2mat,
                          p2, Index.corner, Aimat, B1mat, Cimat, 
                          h.step = 0.0001, colx2.index,
                          xij = NULL, str0 = NULL) {


  nn <- nrow(x2mat)
  if (nrow(Cimat) != p2 || ncol(Cimat) != r)
    stop("'Cimat' wrong shape")

  dct.da <- matrix(NA_real_, (M-r-length(str0))*r, r*p2)

  if ((length(Index.corner) + length(str0)) == M)
    stop("cannot handle full rank models yet")
  cbindex <- (1:M)[-c(Index.corner, str0)]

  ptr <- 1
  for (sss in 1:r)
    for (tt in cbindex) {
      small.Hlist <- vector("list", p2)
      pAmat <- Aimat
      pAmat[tt,sss] <- pAmat[tt,sss] + h.step   # Perturb it
      for (ii in 1:p2)
        small.Hlist[[ii]] <- pAmat

      offset <- if (length(fit@offset)) fit@offset else 0
      if (all(offset == 0))
        offset <- 0
      neweta <- x2mat %*% Cimat %*% t(pAmat)
      if (is.numeric(x1mat))
        neweta <- neweta + x1mat %*% B1mat
      fit@predictors <- neweta


      newmu <- fit@family@linkinv(neweta, fit@extra) 
      fit@fitted.values <- as.matrix(newmu)  # 20100909

      fred <- weights(fit, type = "w", deriv = TRUE, ignore.slot = TRUE)
      if (!length(fred))
        stop("cannot get @weights and @deriv from object")
      wz <- fred$weights
      deriv.mu <- fred$deriv

      U <- vchol(wz, M = M, n = nn, silent = TRUE)
      tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = nn)
      newzmat <- neweta + vbacksub(U, tvfor, M = M, n = nn) - offset
      if (is.numeric(x1mat))
        newzmat <- newzmat - x1mat %*% B1mat

      newfit <- vlm.wfit(xmat = x2mat, zmat = newzmat,
                               Hlist = small.Hlist, U = U,
                               matrix.out = FALSE, is.vlmX = FALSE,
                               ResSS = TRUE, qr = FALSE, x.ret = FALSE,
                               offset = NULL, xij = xij)
      dct.da[ptr, ] <- (newfit$coef - t(Cimat)) / h.step
      ptr <- ptr + 1
    }

    dct.da
}




dctda.fast.only <- function(theta, wz, U, zmat, M, r, x1mat, x2mat,
                            p2, Index.corner, Aimat, B1mat, Cimat,
                            xij = NULL,
                            str0 = NULL) {


  if (length(str0))
    stop("cannot handle 'str0' in dctda.fast.only()")

  nn <- nrow(x2mat)
  if (nrow(Cimat) != p2 || ncol(Cimat) != r)
    stop("Cimat wrong shape")

  fred <- kronecker(matrix(1, 1,r), x2mat)
  fred <- kronecker(fred, matrix(1,M, 1))
  barney <- kronecker(Aimat, matrix(1, 1,p2))
  barney <- kronecker(matrix(1, nn, 1), barney)

  temp <- array(t(barney*fred), c(p2*r, M, nn))
  temp <- aperm(temp, c(2, 1, 3))  # M by p2*r by nn
  temp <- mux5(wz, temp, M = M, matrix.arg= TRUE)
  temp <- m2a(temp, M = p2 * r)  # Note M != M here!
  G <- solve(rowSums(temp, dims = 2))  # p2*r by p2*r 

  dc.da <- array(NA, c(p2, r, M, r))  # different from other functions
  if (length(Index.corner) == M)
      stop("cannot handle full rank models yet")
  cbindex <- (1:M)[-Index.corner]  # complement of Index.corner 
  resid2 <- if (length(x1mat))
    mux22(t(wz), zmat - x1mat %*% B1mat, M = M,
          upper = FALSE, as.matrix = TRUE) else
    mux22(t(wz), zmat                  , M = M,
          upper = FALSE, as.matrix = TRUE)

  for (sss in 1:r)
    for (ttt in cbindex) {
      fred <- t(x2mat) *
              matrix(resid2[, ttt], p2, nn, byrow = TRUE)  # p2 * nn
      temp2 <- kronecker(I.col(sss, r), rowSums(fred))
      for (kkk in 1:r) {
        Wiak <- mux22(t(wz), matrix(Aimat[,kkk], nn, M, byrow = TRUE),
                      M = M, upper = FALSE,
                      as.matrix = TRUE)  # nn * M
        wxx <- Wiak[,ttt] * x2mat
        blocki <- t(x2mat) %*% wxx 
        temp4a <- blocki %*% Cimat[,kkk]
        if (kkk == 1) {
            temp4b <- blocki %*% Cimat[,sss]
        }
        temp2 <- temp2 - kronecker(I.col(sss, r), temp4a) -
                         kronecker(I.col(kkk, r), temp4b)
      }
      dc.da[,,ttt,sss] <- G %*% temp2 
    }
  ans1 <- dc.da[,,cbindex,, drop = FALSE]  # p2 x r x (M-r) x r 
  ans1 <- aperm(ans1, c(2, 1, 3, 4))  # r x p2 x (M-r) x r 

  ans1 <- matrix(c(ans1), r*p2, (M-r)*r)
  ans1 <- t(ans1)
  ans1
}



dcda.fast <- function(theta, wz, U, z, M, r, xmat, pp, Index.corner,
                      intercept = TRUE, xij = NULL) {



  nn <- nrow(xmat)

  Aimat <- matrix(NA_real_, M, r)
  Aimat[Index.corner,] <- diag(r)
  Aimat[-Index.corner,] <- theta    # [-(1:M)]

  if (intercept) {
    Hlist <- vector("list", pp+1)
    Hlist[[1]] <- diag(M)
    for (ii in 2:(pp+1))
      Hlist[[ii]] <- Aimat
  } else {
    Hlist <- vector("list", pp)
    for (ii in 1:pp)
      Hlist[[ii]] <- Aimat
  }

  coeffs <- vlm.wfit(xmat = xmat, z, Hlist, U = U, matrix.out = TRUE,
                     xij = xij)$mat.coef
  c3 <- coeffs <- t(coeffs)  # transpose to make M x (pp+1)


  int.vec <- if (intercept) c3[, 1] else 0  # \boldeta_0
  Cimat <- if (intercept) t(c3[Index.corner,-1, drop = FALSE]) else 
           t(c3[Index.corner,, drop = FALSE])
  if (nrow(Cimat)!=pp || ncol(Cimat)!=r)
    stop("Cimat wrong shape")

  fred <- kronecker(matrix(1, 1,r),
                    if (intercept) xmat[,-1, drop = FALSE] else xmat)
  fred <- kronecker(fred, matrix(1,M, 1))
  barney <- kronecker(Aimat, matrix(1, 1,pp))
  barney <- kronecker(matrix(1, nn, 1), barney)

  temp <- array(t(barney*fred), c(r*pp,M,nn))
  temp <- aperm(temp, c(2, 1, 3))
  temp <- mux5(wz, temp, M = M, matrix.arg = TRUE)
  temp <- m2a(temp, M = r * pp)     # Note M != M here!
  G <- solve(rowSums(temp, dims = 2))

  dc.da <- array(NA, c(pp,r,M,r))  # different from other functions
  cbindex <- (1:M)[-Index.corner]
  resid2 <- mux22(t(wz),
                  z - matrix(int.vec, nn, M, byrow = TRUE), M = M,
                  upper = FALSE, as.matrix = TRUE)  # mat = TRUE,

  for (s in 1:r)
    for (tt in cbindex) {
      fred <- (if (intercept) t(xmat[, -1, drop = FALSE]) else
               t(xmat)) * matrix(resid2[, tt], pp, nn, byrow = TRUE) 
      temp2 <- kronecker(I.col(s, r), rowSums(fred))

      temp4 <- rep(0,pp)
      for (k in 1:r) {
        Wiak <- mux22(t(wz),
                      matrix(Aimat[, k], nn, M, byrow = TRUE),
                      M = M, upper = FALSE, as.matrix = TRUE)
        wxx <- Wiak[,tt] * (if (intercept)
                            xmat[, -1, drop = FALSE] else
                            xmat)
        blocki <- (if (intercept)
                  t(xmat[, -1, drop = FALSE]) else
                  t(xmat)) %*% wxx
        temp4 <- temp4 + blocki %*% Cimat[, k]
      }
      dc.da[,,tt,s] <- G %*% (temp2 - 2 * kronecker(I.col(s, r), temp4))
    }
  ans1 <- dc.da[,,cbindex,, drop = FALSE]  # pp x r x (M-r) x r 
  ans1 <- aperm(ans1, c(2, 1, 3, 4))   # r x pp x (M-r) x r 

  ans1 <- matrix(c(ans1), (M-r)*r, r*pp, byrow = TRUE)


  detastar.da <- array(0,c(M,r,r,nn))
  for (s in 1:r)
    for (j in 1:r) {
      t1 <- t(dc.da[,j,,s])
      t1 <- matrix(t1, M, pp)
      detastar.da[,j,s,] <- t1 %*% (if (intercept)
                            t(xmat[,-1, drop = FALSE]) else t(xmat))
    }

  etastar <- (if (intercept) xmat[,-1, drop = FALSE] else xmat) %*% Cimat
  eta <- matrix(int.vec, nn, M, byrow = TRUE) + etastar %*% t(Aimat)

  sumWinv <- solve((m2a(t(colSums(wz)), M = M))[, , 1])

  deta0.da <- array(0,c(M,M,r))
  AtWi <- kronecker(matrix(1, nn, 1), Aimat)
  AtWi <- mux111(t(wz), AtWi, M = M, upper= FALSE)  # matrix.arg= TRUE, 
  AtWi <- array(t(AtWi), c(r, M, nn))
  for (ss in 1:r) {
    temp90 <- (m2a(t(colSums(etastar[, ss]*wz)), M = M))[, , 1]  # MxM
    temp92 <- array(detastar.da[,,ss,], c(M, r, nn))
    temp93 <- mux7(temp92, AtWi)
    temp91 <- rowSums(temp93, dims = 2)  # M x M
    deta0.da[,,ss] <- -(temp90 + temp91) %*% sumWinv
  }
  ans2 <- deta0.da[-(1:r), , , drop = FALSE]  # (M-r) x M x r
  ans2 <- aperm(ans2, c(1, 3, 2))  # (M-r) x r x M
  ans2 <- matrix(c(ans2), (M-r)*r, M) 

  list(dc.da = ans1, dint.da = ans2)
}



rrr.deriv.ResSS <- function(theta, wz, U, z, M, r, xmat,
                            pp, Index.corner, intercept = TRUE,
                            xij = NULL) {

  Amat <- matrix(NA_real_, M, r)
  Amat[Index.corner,] <- diag(r)
  Amat[-Index.corner,] <- theta    # [-(1:M)]

  if (intercept) {
    Hlist <- vector("list", pp+1)
    Hlist[[1]] <- diag(M)
    for (ii in 2:(pp+1))
      Hlist[[ii]] <- Amat
  } else {
    Hlist <- vector("list", pp)
    for (ii in 1:pp)
      Hlist[[ii]] <- Amat
  }

  vlm.wfit(xmat = xmat, z, Hlist, U = U, matrix.out = FALSE,
           ResSS = TRUE, xij = xij)$ResSS
}




rrr.deriv.gradient.fast <- function(theta, wz, U, z, M, r, xmat,
                                    pp, Index.corner,
                                    intercept = TRUE) {




  nn <- nrow(xmat)

  Aimat <- matrix(NA_real_, M, r)
  Aimat[Index.corner,] <- diag(r)
  Aimat[-Index.corner,] <- theta    # [-(1:M)]

  if (intercept) {
    Hlist <- vector("list", pp+1)
    Hlist[[1]] <- diag(M)
    for (i in 2:(pp+1))
      Hlist[[i]] <- Aimat
  } else {
    Hlist <- vector("list", pp)
    for (i in 1:(pp))
      Hlist[[i]] <- Aimat
  }

  coeffs <- vlm.wfit(xmat, z, Hlist, U = U, matrix.out= TRUE,
                     xij = NULL)$mat.coef
  c3 <- coeffs <- t(coeffs)  # transpose to make M x (pp+1)


  int.vec <- if (intercept) c3[, 1] else 0  # \boldeta_0
  Cimat <- if (intercept) t(c3[Index.corner, -1, drop = FALSE]) else
           t(c3[Index.corner,, drop = FALSE])
  if (nrow(Cimat) != pp || ncol(Cimat) != r)
      stop("Cimat wrong shape")

  fred <- kronecker(matrix(1, 1,r),
                    if (intercept) xmat[, -1, drop = FALSE] else xmat)
  fred <- kronecker(fred, matrix(1, M, 1))
  barney <- kronecker(Aimat, matrix(1, 1, pp))
  barney <- kronecker(matrix(1, nn, 1), barney)

  temp <- array(t(barney*fred), c(r*pp, M, nn))
  temp <- aperm(temp, c(2, 1, 3))
  temp <- mux5(wz, temp, M = M, matrix.arg = TRUE)
  temp <- m2a(temp, M = r * pp)  # Note M != M here!
  G <- solve(rowSums(temp, dims = 2))

  dc.da <- array(NA,c(pp,r,r,M))
  cbindex <- (1:M)[-Index.corner]
  resid2 <- mux22(t(wz), z - matrix(int.vec, nn, M, byrow = TRUE),
                  M = M,
                  upper = FALSE, as.matrix = TRUE)

  for (s in 1:r)
    for (tt in cbindex) {
      fred <- (if (intercept) t(xmat[, -1, drop = FALSE]) else
               t(xmat)) * matrix(resid2[, tt], pp, nn, byrow = TRUE) 
      temp2 <- kronecker(I.col(s, r), rowSums(fred))

      temp4 <- rep(0,pp)
      for (k in 1:r) {
        Wiak <- mux22(t(wz),
                     matrix(Aimat[, k], nn, M, byrow = TRUE),
                     M = M, upper = FALSE, as.matrix = TRUE)
        wxx <- Wiak[,tt] * (if (intercept)
                            xmat[, -1, drop = FALSE] else xmat)
        blocki <- (if (intercept) t(xmat[, -1, drop = FALSE]) else
                  t(xmat)) %*% wxx 
        temp4 <- temp4 + blocki %*% Cimat[, k]
      }
      dc.da[,,s,tt] <- G %*% (temp2 - 2 * kronecker(I.col(s, r), temp4))
    }

  detastar.da <- array(0,c(M,r,r,nn))
  for (s in 1:r)
    for (j in 1:r) {
      t1 <- t(dc.da[,j,s,])
      t1 <- matrix(t1, M, pp)
      detastar.da[,j,s,] <- t1 %*% (if (intercept)
                            t(xmat[, -1, drop = FALSE]) else t(xmat))
    }

  etastar <- (if (intercept) xmat[, -1, drop = FALSE] else xmat) %*% Cimat
  eta <- matrix(int.vec, nn, M, byrow = TRUE) + etastar %*% t(Aimat)

  sumWinv <- solve((m2a(t(colSums(wz)), M = M))[, , 1])

  deta0.da <- array(0, c(M, M, r))

  AtWi <- kronecker(matrix(1, nn, 1), Aimat)
  AtWi <- mux111(t(wz), AtWi, M = M, upper = FALSE)  # matrix.arg= TRUE, 
  AtWi <- array(t(AtWi), c(r, M, nn))

  for (ss in 1:r) {
    temp90 <- (m2a(t(colSums(etastar[, ss] * wz)), M = M))[, , 1]
    temp92 <- array(detastar.da[, , ss, ], c(M, r, nn))
    temp93 <- mux7(temp92,AtWi)
    temp91 <- apply(temp93, 1:2,sum)  # M x M
    temp91 <- rowSums(temp93, dims = 2)  # M x M
    deta0.da[,,ss] <- -(temp90 + temp91) %*% sumWinv
  }

  ans <- matrix(0,M,r)
  fred <- mux22(t(wz), z - eta, M = M,
                upper = FALSE, as.matrix = TRUE)
  fred.array <- array(t(fred %*% Aimat),c(r, 1, nn))
  for (s in 1:r) {
    a1 <- colSums(fred %*% t(deta0.da[,, s]))
    a2 <- colSums(fred * etastar[, s])
    temp92 <- array(detastar.da[, , s, ],c(M, r, nn))
    temp93 <- mux7(temp92, fred.array)
    a3 <- rowSums(temp93, dims = 2)
    ans[,s] <- a1 + a2 + a3
  }

  ans <- -2 * c(ans[cbindex, ])

  ans
}








vellipse <- function(R, ratio = 1, orientation = 0,
                     center = c(0, 0), N = 300) {
  if (length(center) != 2)
    stop("argument 'center' must be of length 2")
  theta <-       2*pi*(0:N)/N
  x1 <-       R*cos(theta)
  y1 <- ratio*R*sin(theta)
  x <- center[1] + cos(orientation)*x1 - sin(orientation)*y1
  y <- center[2] + sin(orientation)*x1 + cos(orientation)*y1
  cbind(x, y)
}


biplot.qrrvglm <- function(x, ...) {
  stop("biplot.qrrvglm has been replaced by the function lvplot.qrrvglm")
}



 lvplot.qrrvglm <-
  function(object, varI.latvar = FALSE, refResponse = NULL,
           add = FALSE, show.plot = TRUE, rug = TRUE, y = FALSE, 
           type = c("fitted.values", "predictors"),
           xlab = paste("Latent Variable",
                        if (Rank == 1) "" else " 1", sep = ""),
           ylab = if (Rank == 1) switch(type, predictors = "Predictors", 
              fitted.values = "Fitted values") else "Latent Variable 2",
          pcex = par()$cex, pcol = par()$col, pch = par()$pch, 
          llty = par()$lty, lcol = par()$col, llwd = par()$lwd,
          label.arg = FALSE, adj.arg = -0.1, 
          ellipse = 0.95, Absolute = FALSE, 
              elty = par()$lty, ecol = par()$col, elwd = par()$lwd,
              egrid = 200,
          chull.arg = FALSE, clty = 2, ccol = par()$col, clwd = par()$lwd,
              cpch = "   ",
          C = FALSE,
              OriginC = c("origin", "mean"),
              Clty = par()$lty, Ccol = par()$col, Clwd = par()$lwd,
              Ccex = par()$cex, Cadj.arg = -0.1, stretchC = 1, 
          sites = FALSE, spch = NULL, scol = par()$col, scex = par()$cex,
          sfont = par()$font,
          check.ok = TRUE, ...) {
    if (mode(type) != "character" && mode(type) != "name")
      type <- as.character(substitute(type))
    type <- match.arg(type, c("fitted.values", "predictors"))[1]

    if (is.numeric(OriginC))
      OriginC <- rep(OriginC, length.out = 2) else {
      if (mode(OriginC) != "character" && mode(OriginC) != "name")
        OriginC <- as.character(substitute(OriginC))
      OriginC <- match.arg(OriginC, c("origin","mean"))[1]
    }

    if (length(ellipse) > 1)
      stop("ellipse must be of length 1 or 0")
    if (is.logical(ellipse)) {ellipse <- if (ellipse) 0.95 else NULL}

    Rank <- object@control$Rank
    if (Rank > 2)
      stop("can only handle rank 1 or 2 models")
    M <- object@misc$M
    NOS <- ncol(object@y)
    MSratio <- M / NOS  # First value is g(mean) = quadratic form in latvar
    n <- object@misc$n
    colx2.index <- object@control$colx2.index
    cx1i <- object@control$colx1.index  # May be NULL
    if (check.ok)
      if (!(length(cx1i) == 1 && names(cx1i) == "(Intercept)"))
        stop("latent variable plots allowable only for ",
             "noRRR = ~ 1 models")

    Coef.list <- Coef(object, varI.latvar = varI.latvar,
                      refResponse = refResponse)
    if ( C) Cmat <- Coef.list@C
    nustar <- Coef.list@latvar  # n x Rank 

    if (!show.plot) return(nustar)

    r.curves <- slot(object, type)   # n times M (\boldeta or \boldmu) 
    if (!add) {
      if (Rank == 1) {
        matplot(nustar,
                if ( y && type == "fitted.values")
                object@y else r.curves,
                type = "n", xlab = xlab, ylab = ylab, ...)
      } else {  # Rank == 2
        matplot(c(Coef.list@Optimum[1, ], nustar[, 1]),
                c(Coef.list@Optimum[2, ], nustar[, 2]),
                type = "n", xlab = xlab, ylab = ylab, ...)
      }
    }




    pch  <- rep(pch,  length = ncol(r.curves))
    pcol <- rep(pcol, length = ncol(r.curves))
    pcex <- rep(pcex, length = ncol(r.curves))
    llty <- rep(llty, length = ncol(r.curves))
    lcol <- rep(lcol, length = ncol(r.curves))
    llwd <- rep(llwd, length = ncol(r.curves))
    elty <- rep(elty, length = ncol(r.curves))
    ecol <- rep(ecol, length = ncol(r.curves))
    elwd <- rep(elwd, length = ncol(r.curves))
    adj.arg <- rep(adj.arg, length = ncol(r.curves))
    if ( C ) {
      Clwd <- rep(Clwd, length = nrow(Cmat))
      Clty <- rep(Clty, length = nrow(Cmat))
      Ccol <- rep(Ccol, length = nrow(Cmat))
      Ccex <- rep(Ccex, length = nrow(Cmat))
      Cadj.arg <- rep(Cadj.arg, length = nrow(Cmat))
    }

    if (Rank == 1) {
      for (i in 1:ncol(r.curves)) {
        xx <- nustar 
        yy <- r.curves[,i]
        o <- sort.list(xx)
        xx <- xx[o]
        yy <- yy[o]
        lines(xx, yy, col = lcol[i], lwd = llwd[i], lty = llty[i])
        if ( y && type == "fitted.values") {
          ypts <- object@y
          if (ncol(as.matrix(ypts)) == ncol(r.curves))
            points(xx, ypts[o,i], col = pcol[i],
                   cex = pcex[i], pch = pch[i])
        } 
      } 
      if (rug)
        rug(xx) 
    } else {
     for (i in 1:ncol(r.curves))
      points(Coef.list@Optimum[1, i], Coef.list@Optimum[2, i],
             col = pcol[i], cex = pcex[i], pch = pch[i])
     if (label.arg) {
      for (i in 1:ncol(r.curves))
          text(Coef.list@Optimum[1, i], Coef.list@Optimum[2, i],
               labels = (dimnames(Coef.list@Optimum)[[2]])[i], 
               adj = adj.arg[i], col = pcol[i], cex = pcex[i])
    }
    if (chull.arg) {
      hull <- chull(nustar[, 1], nustar[, 2])
      hull <- c(hull, hull[1])
      lines(nustar[hull, 1], nustar[hull, 2], type = "b", pch = cpch,
            lty = clty, col = ccol, lwd = clwd)
    }
    if (length(ellipse)) {
      ellipse.temp <- if (ellipse > 0) ellipse else 0.95
      if (ellipse < 0 && (!object@control$eq.tolerances || varI.latvar))
        stop("an equal-tolerances assumption and 'varI.latvar = FALSE' ",
             "is needed for 'ellipse' < 0")
      if ( check.ok ) {
        colx1.index <- object@control$colx1.index
        if (!(length(colx1.index) == 1 &&
              names(colx1.index) == "(Intercept)"))
          stop("can only plot ellipses for intercept models only")
      }
      for (i in 1:ncol(r.curves)) {
        cutpoint <- object@family@linkfun( if (Absolute) ellipse.temp
                        else Coef.list@Maximum[i] * ellipse.temp,
                        extra = object@extra)
        if (MSratio > 1) 
          cutpoint <- cutpoint[1, 1]

          cutpoint <- object@family@linkfun(Coef.list@Maximum[i],
                      extra = object@extra) - cutpoint
          if (is.finite(cutpoint) && cutpoint > 0) {
            Mmat <- diag(rep(ifelse(object@control$Crow1positive, 1, -1),
                            length.out = Rank))
            etoli <- eigen(t(Mmat) %*% Coef.list@Tolerance[,,i] %*% Mmat)
            A <- ifelse(etoli$val[1]>0, sqrt(2*cutpoint*etoli$val[1]), Inf)
            B <- ifelse(etoli$val[2]>0, sqrt(2*cutpoint*etoli$val[2]), Inf)
            if (ellipse < 0)
              A <- B <- -ellipse / 2

            theta.angle <- asin(etoli$vector[2, 1]) *
                           ifelse(object@control$Crow1positive[2], 1, -1)
            if (object@control$Crow1positive[1])
              theta.angle <- pi - theta.angle
            if (all(is.finite(c(A,B))))
              lines(vellipse(R = 2*A, ratio = B/A,
                             orientation = theta.angle,
                             center = Coef.list@Optimum[, i],
                             N = egrid),
                    lwd = elwd[i], col =ecol[i], lty = elty[i])
            }
        }
      }

    if ( C ) {
      if (is.character(OriginC) && OriginC == "mean")
        OriginC <- c(mean(nustar[, 1]), mean(nustar[, 2]))
      if (is.character(OriginC) && OriginC == "origin")
        OriginC <- c(0,0)
      for (i in 1:nrow(Cmat))
        arrows(x0 = OriginC[1], y0 = OriginC[2],
               x1 = OriginC[1] + stretchC * Cmat[i, 1],
               y1 = OriginC[2] + stretchC * Cmat[i, 2],
               lty = Clty[i], col = Ccol[i], lwd = Clwd[i])
      if (label.arg) {
        temp200 <- dimnames(Cmat)[[1]]
        for (i in 1:nrow(Cmat))
          text(OriginC[1] + stretchC * Cmat[i, 1],
               OriginC[2] + stretchC * Cmat[i, 2], col = Ccol[i],
               labels = temp200[i], adj = Cadj.arg[i],
               cex = Ccex[i])
      }
    }
    if (sites) {
      text(nustar[, 1], nustar[, 2], adj = 0.5,
           labels = if (is.null(spch)) dimnames(nustar)[[1]] else
           rep(spch, length = nrow(nustar)), col = scol,
           cex = scex, font = sfont)
    }
  }
  invisible(nustar)
}



lvplot.rrvglm <- function(object,
                          A = TRUE,
                          C = TRUE,
                          scores = FALSE, show.plot = TRUE,
                          groups = rep(1, n),
                          gapC = sqrt(sum(par()$cxy^2)), scaleA = 1,
                          xlab = "Latent Variable 1",
                          ylab = "Latent Variable 2",
         Alabels= if (length(object@misc$predictors.names))
         object@misc$predictors.names else paste("LP", 1:M, sep = ""),
                          Aadj = par()$adj,
                          Acex = par()$cex,
                          Acol = par()$col,
                          Apch = NULL,
                          Clabels=rownames(Cmat),
                          Cadj = par()$adj,
                          Ccex = par()$cex,
                          Ccol = par()$col, 
                          Clty = par()$lty, 
                          Clwd = par()$lwd, 
                          chull.arg = FALSE,
                          ccex = par()$cex,
                          ccol = par()$col,
                          clty = par()$lty,
                          clwd = par()$lwd,
                          spch = NULL,
                          scex = par()$cex,
                          scol = par()$col,
                          slabels = rownames(x2mat),
                          ...) {


    if (object@control$Rank != 2 && show.plot)
        stop("can only handle rank-2 models")
    M <- object@misc$M
    n <- object@misc$n
    colx2.index <- object@control$colx2.index
    Coef.list <- Coef(object)
    Amat <- Coef.list@A
    Cmat <- Coef.list@C

    Amat <- Amat * scaleA
    dimnames(Amat) <- list(object@misc$predictors.names, NULL) 
    Cmat <- Cmat / scaleA

    if (!length(object@x)) {
        object@x <- model.matrixvlm(object, type = "lm")
    }
    x2mat <- object@x[, colx2.index, drop = FALSE]
    nuhat <- x2mat %*% Cmat
    if (!show.plot) return(as.matrix(nuhat))

    index.nosz <- 1:M
    allmat <- rbind(if (A) Amat else NULL, 
                   if (C) Cmat else NULL, 
                   if (scores) nuhat else NULL)

    plot(allmat[, 1], allmat[, 2], type = "n",
         xlab=xlab, ylab=ylab, ...)  # xlim etc. supplied through ...

    if (A) {
        Aadj <- rep(Aadj, length.out = length(index.nosz))
        Acex <- rep(Acex, length.out = length(index.nosz))
        Acol <- rep(Acol, length.out = length(index.nosz))
        if (length(Alabels) != M)
          stop("'Alabels' must be of length ", M)
        if (length(Apch)) {
            Apch <- rep(Apch, length.out = length(index.nosz))
            for (i in index.nosz)
                points(Amat[i, 1],
                       Amat[i, 2],
                       pch=Apch[i],cex = Acex[i],col=Acol[i])
        } else {
            for (i in index.nosz)
                text(Amat[i, 1], Amat[i, 2],
                     Alabels[i], cex = Acex[i],
                     col =Acol[i], adj=Aadj[i])
        }
    }

    if (C) {
      p2 <- nrow(Cmat)
      gapC <- rep(gapC, length.out = p2)
      Cadj <- rep(Cadj, length.out = p2)
      Ccex <- rep(Ccex, length.out = p2)
      Ccol <- rep(Ccol, length.out = p2)
      Clwd <- rep(Clwd, length.out = p2)
      Clty <- rep(Clty, length.out = p2)
      if (length(Clabels) != p2)
        stop("'length(Clabels)' must be equal to ", p2)
      for (ii in 1:p2) {
        arrows(0, 0, Cmat[ii, 1], Cmat[ii, 2],
               lwd = Clwd[ii], lty = Clty[ii], col = Ccol[ii])
        const <- 1 + gapC[ii] / sqrt(Cmat[ii, 1]^2 + Cmat[ii, 2]^2)
        text(const*Cmat[ii, 1], const*Cmat[ii, 2],
             Clabels[ii], cex = Ccex[ii],
             adj = Cadj[ii], col = Ccol[ii])
      }
    }

    if (scores) {
      ugrp <- unique(groups)
      nlev <- length(ugrp)  # number of groups
      clty <- rep(clty, length.out = nlev)
      clwd <- rep(clwd, length.out = nlev)
      ccol <- rep(ccol, length.out = nlev)
      if (length(spch))
        spch <- rep(spch, length.out = n)
      scol <- rep(scol, length.out = n)
      scex <- rep(scex, length.out = n)
      for (ii in ugrp) {
        gp <- groups == ii
        if (nlev > 1 && (length(unique(spch[gp])) != 1 ||
            length(unique(scol[gp])) != 1 ||
            length(unique(scex[gp])) != 1))
          warning("spch/scol/scex is different for individuals ",
                  "from the same group")

        temp <- nuhat[gp,, drop = FALSE]
        if (length(spch)) {
          points(temp[, 1], temp[, 2], cex = scex[gp], pch = spch[gp],
                 col = scol[gp])
        } else {
            text(temp[, 1], temp[, 2], label = slabels, cex = scex[gp],
                 col = scol[gp])
        }
        if (chull.arg) {
          hull <- chull(temp[, 1], temp[, 2])
          hull <- c(hull, hull[1])
          lines(temp[hull, 1], temp[hull, 2],
                type = "b", lty = clty[ii],
                col = ccol[ii], lwd = clwd[ii], pch = "  ")
         }
      }
    }

    invisible(nuhat)
}






 Coef.rrvglm <- function(object, ...) {
    M <- object@misc$M
    n <- object@misc$n
    colx1.index <- object@control$colx1.index
    colx2.index <- object@control$colx2.index
    p1 <- length(colx1.index)  # May be 0
    Amat <- object@constraints[[colx2.index[1]]]

    B1mat <- if (p1)
      coefvlm(object, matrix.out = TRUE)[colx1.index,, drop = FALSE] else
      NULL


    C.try <- coefvlm(object, matrix.out = TRUE)[colx2.index, , drop = FALSE]


    Cmat <- C.try %*% Amat %*% solve(t(Amat) %*% Amat)


    Rank <- object@control$Rank
    latvar.names <- if (Rank > 1)
                      paste("latvar", 1:Rank, sep = "") else
                            "latvar"
    dimnames(Amat) <- list(object@misc$predictors.names, latvar.names)
    dimnames(Cmat) <- list(dimnames(Cmat)[[1]], latvar.names)

    ans <- new(Class = "Coef.rrvglm",
      A            = Amat,
      C            = Cmat,
      Rank         = Rank,
      colx2.index  = colx2.index)

  if (!is.null(colx1.index)) {
    ans@colx1.index  <- colx1.index
    ans@B1 <- B1mat
  }

  if (object@control$Corner)
    ans@Atilde <- Amat[-c(object@control$Index.corner,
                       object@control$str0),, drop = FALSE]
  ans
}




setMethod("Coef", "rrvglm",
          function(object, ...) Coef.rrvglm(object, ...))





show.Coef.rrvglm <- function(x, ...) {

    object <- x

  cat("A matrix:\n")
  print(object@A, ...)

  cat("\nC matrix:\n")
  print(object@C, ...)

  p1 <- length(object@colx1.index)
  if (p1) {
    cat("\nB1 matrix:\n")
    print(object@B1, ...)
  }

  invisible(object)
} 


 if (!isGeneric("biplot"))
    setGeneric("biplot", function(x, ...) standardGeneric("biplot")) 


setMethod("Coef", "qrrvglm", function(object, ...)
          Coef.qrrvglm(object, ...))



setMethod("biplot", "qrrvglm",
           function(x, ...) {
           biplot.qrrvglm(x, ...)})

setMethod("lvplot", "qrrvglm",
           function(object, ...) {
           invisible(lvplot.qrrvglm(object, ...))})

setMethod("lvplot", "rrvglm",
           function(object, ...) {
           invisible(lvplot.rrvglm(object, ...))})


biplot.rrvglm <- function(x, ...)
    lvplot(object = x, ...)

setMethod("biplot",  "rrvglm", function(x, ...)
           invisible(biplot.rrvglm(x, ...)))




summary.qrrvglm <-
  function(object,
           varI.latvar = FALSE, refResponse = NULL, ...) {
    answer <- object
    answer@post$Coef <- Coef(object,
                             varI.latvar = varI.latvar,
                             refResponse = refResponse, 
                             ...)  # Store it here; non-elegant

  if (length((answer@post$Coef)@dispersion) &&
     length(object@misc$estimated.dispersion) &&
     object@misc$estimated.dispersion)
      answer@dispersion <- 
      answer@misc$dispersion <- (answer@post$Coef)@dispersion

  as(answer, "summary.qrrvglm")
}



show.summary.qrrvglm <- function(x, ...) {



  cat("\nCall:\n")
  dput(x@call)

  print(x@post$Coef, ...)  # non-elegant programming

  if (length(x@dispersion) > 1) {
    cat("\nDispersion parameters:\n")
    if (length(x@misc$ynames)) {
      names(x@dispersion) <- x@misc$ynames 
      print(x@dispersion, ...)
    } else {
      cat(x@dispersion, fill = TRUE)
    }
    cat("\n")
  } else if (length(x@dispersion) == 1) {
    cat("\nDispersion parameter:  ", x@dispersion, "\n")
  }

}


 setClass("summary.qrrvglm", contains = "qrrvglm")







setMethod("summary", "qrrvglm",
          function(object, ...)
          summary.qrrvglm(object, ...))





 setMethod("show", "summary.qrrvglm",
           function(object)
           show.summary.qrrvglm(object))





setMethod("show", "Coef.rrvglm", function(object)
          show.Coef.rrvglm(object))





 grc <- function(y, Rank = 1, Index.corner = 2:(1+Rank),
                 str0 = 1,
                 summary.arg = FALSE, h.step = 0.0001, ...) {
                           


    myrrcontrol <- rrvglm.control(Rank = Rank,
                                  Index.corner = Index.corner,
                                  str0 = str0, ...)
    object.save <- y
    if (is(y, "rrvglm")) {
      y <- object.save@y
    } else {
      y <- as.matrix(y)
      y <- as(y, "matrix")
    }
    if (length(dim(y)) != 2 || nrow(y) < 3 || ncol(y) < 3)
     stop("y must be a matrix with >= 3 rows & columns, ",
          "or a rrvglm() object")

    ei <- function(i, n) diag(n)[, i, drop = FALSE]
    .grc.df <- data.frame(Row.2 = I.col(2, nrow(y)))

    yn1 <- if (length(dimnames(y)[[1]])) dimnames(y)[[1]] else
              paste("X2.", 1:nrow(y), sep = "")
    warn.save <- options()$warn
    options(warn = -3)  # Suppress the warnings (hopefully, temporarily)
    if (any(!is.na(as.numeric(substring(yn1, 1, 1)))))
        yn1 <- paste("X2.", 1:nrow(y), sep = "")
    options(warn = warn.save)

    Row. <- factor(1:nrow(y))
    modmat.row <- model.matrix( ~ Row.)
    Col. <- factor(1:ncol(y))
    modmat.col <- model.matrix( ~ Col.)

    cms <- list("(Intercept)" = matrix(1, ncol(y), 1))
    for (ii in 2:nrow(y)) {
            cms[[paste("Row.", ii, sep = "")]] <- matrix(1, ncol(y), 1)
        .grc.df[[paste("Row.", ii, sep = "")]] <- modmat.row[,ii]
    }
    for (ii in 2:ncol(y)) {
            cms[[paste("Col.", ii, sep = "")]] <-
               modmat.col[,ii, drop = FALSE]
        .grc.df[[paste("Col.", ii, sep = "")]] <- rep(1, nrow(y))
    }
    for (ii in 2:nrow(y)) {
            cms[[yn1[ii]]] <- diag(ncol(y))
        .grc.df[[yn1[ii]]] <- I.col(ii, nrow(y))
    }

    dimnames(.grc.df) <- list(if (length(dimnames(y)[[1]]))
                              dimnames(y)[[1]] else 
                              as.character(1:nrow(y)),
                              dimnames(.grc.df)[[2]])

    str1 <- "~ Row.2"
    if (nrow(y) > 2)
      for (ii in 3:nrow(y))
        str1 <- paste(str1, paste("Row.", ii, sep = ""), sep = " + ")
    for (ii in 2:ncol(y))
      str1 <- paste(str1, paste("Col.", ii, sep = ""), sep = " + ")
    str2 <- paste("y ", str1)
    for (ii in 2:nrow(y))
      str2 <- paste(str2, yn1[ii], sep = " + ")
    myrrcontrol$noRRR <- as.formula(str1)  # Overwrite this

    assign(".grc.df", .grc.df, envir = VGAMenv)

    warn.save <- options()$warn
    options(warn = -3)  # Suppress the warnings (hopefully, temporarily)
    answer <- if (is(object.save, "rrvglm")) object.save else 
              rrvglm(as.formula(str2), family = poissonff,
                     constraints = cms, control = myrrcontrol,
                     data = .grc.df)
    options(warn = warn.save)

    if (summary.arg) {
      answer <- as(answer, "rrvglm")

      answer <- summary.rrvglm(answer, h.step = h.step)
    } else { 
      answer <- as(answer, "grc")
    }

    if (exists(".grc.df", envir = VGAMenv))
      rm(".grc.df", envir = VGAMenv)

    answer
}

summary.grc <- function(object, ...) {
    grc(object, summary.arg= TRUE, ...)
}





trplot.qrrvglm <-
  function(object,
           which.species = NULL,
           add = FALSE, show.plot = TRUE,
           label.sites = FALSE, 
           sitenames = rownames(object@y),
           axes.equal = TRUE,
           cex = par()$cex,
           col = 1:(nos*(nos-1)/2),
           log = "", 
           lty = rep(par()$lty, length.out = nos*(nos-1)/2),
           lwd = rep(par()$lwd, length.out = nos*(nos-1)/2),
           tcol = rep(par()$col, length.out = nos*(nos-1)/2),
           xlab = NULL, ylab = NULL, 
           main = "",   # "Trajectory plot",
           type = "b",
           check.ok = TRUE, ...) {
  coef.obj <- Coef(object)  # use defaults for those two arguments
  if (coef.obj@Rank != 1)
    stop("object must be a rank-1 model")
  fv <- fitted(object)
  modelno <- object@control$modelno  # 1, 2, 3, or 0
  NOS <- ncol(fv)   # Number of species
  M <- object@misc$M  #
  nn <- nrow(fv)  # Number of sites 
  if (length(sitenames))
    sitenames <- rep(sitenames, length.out = nn)
  sppNames <- dimnames(object@y)[[2]]
  if (!length(which.species)) {
    which.species <- sppNames[1:NOS]
    which.species.numer <- 1:NOS
  } else
  if (is.numeric(which.species)) {
    which.species.numer <- which.species
    which.species <- sppNames[which.species.numer]  # Convert to character
  } else {
     which.species.numer <- match(which.species, sppNames)
  }
    nos <- length(which.species)  # nos = number of species to be plotted

  if (length(which.species.numer) <= 1)
    stop("must have at least 2 species to be plotted")
  cx1i <- object@control$colx1.index
  if (check.ok)
  if (!(length(cx1i) == 1 && names(cx1i) == "(Intercept)"))
    stop("trajectory plots allowable only for noRRR = ~ 1 models")

  first.spp  <- iam(1, 1,M = M,both = TRUE,diag = FALSE)$row.index
  second.spp <- iam(1, 1,M = M,both = TRUE,diag = FALSE)$col.index
  myxlab <- if (length(which.species.numer) == 2) {
              paste("Fitted value for",
              if (is.character(which.species.numer))
                  which.species.numer[1] else
                  sppNames[which.species.numer[1]])
               } else "Fitted value for 'first' species"
  myxlab <- if (length(xlab)) xlab else myxlab
  myylab <- if (length(which.species.numer) == 2) {
              paste("Fitted value for",
              if (is.character(which.species.numer))
                  which.species.numer[2] else
                  sppNames[which.species.numer[2]])
               } else "Fitted value for 'second' species"
  myylab <- if (length(ylab)) ylab else myylab
  if (!add) {
    xxx <- if (axes.equal) fv[,which.species.numer] else
           fv[,which.species.numer[first.spp]]
    yyy <- if (axes.equal) fv[,which.species.numer] else
           fv[,which.species.numer[second.spp]]
    matplot(xxx, yyy, type = "n", log = log, xlab = myxlab,
            ylab = myylab, main = main, ...)
  }

  lwd  <- rep(lwd,  length.out = nos*(nos-1)/2)
  col  <- rep(col,  length.out = nos*(nos-1)/2)
  lty  <- rep(lty,  length.out = nos*(nos-1)/2)
  tcol <- rep(tcol, length.out = nos*(nos-1)/2)

  oo <- order(coef.obj@latvar)  # Sort by the latent variable
  ii <- 0
  col <- rep(col, length = nos*(nos-1)/2)
  species.names <- NULL
  if (show.plot)
    for (i1 in seq(which.species.numer)) {
      for (i2 in seq(which.species.numer))
        if (i1 < i2) {
          ii <- ii + 1
          species.names <- rbind(species.names,
                                 cbind(sppNames[i1], sppNames[i2]))
          matplot(fv[oo, which.species.numer[i1]],
                  fv[oo, which.species.numer[i2]],
                  type = type, add = TRUE,
                  lty = lty[ii], lwd = lwd[ii], col = col[ii],
                  pch = if (label.sites) "   " else "*" )
          if (label.sites && length(sitenames))
              text(fv[oo, which.species.numer[i1]],
                   fv[oo, which.species.numer[i2]],
                   labels = sitenames[oo], cex = cex, col = tcol[ii])
        }
    }
  invisible(list(species.names = species.names, 
                 sitenames     = sitenames[oo]))
}

 if (!isGeneric("trplot"))
     setGeneric("trplot",
                function(object, ...) standardGeneric("trplot"))
setMethod("trplot", "qrrvglm",
    function(object, ...) trplot.qrrvglm(object, ...))

setMethod("trplot", "rrvgam",
    function(object, ...) trplot.qrrvglm(object, ...))




vcovrrvglm <- function(object, ...) {
  summary.rrvglm(object, ...)@cov.unscaled 
}



vcovqrrvglm <- function(object,
                       I.tolerances = object@control$eq.tolerances,
                       MaxScale = c("predictors", "response"),
           dispersion = rep(if (length(sobj@dispersion))
           sobj@dispersion else 1,
                            length.out = M), ...) {
  stop("this function is not yet completed")

  if (mode(MaxScale) != "character" && mode(MaxScale) != "name")
    MaxScale <- as.character(substitute(MaxScale))
  MaxScale <- match.arg(MaxScale, c("predictors", "response"))[1]
  if (MaxScale != "predictors")
    stop("can currently only handle MaxScale='predictors'")

  sobj <- summary(object)
  cobj <- Coef(object, I.tolerances = I.tolerances, ...)
  M <- nrow(cobj@A)
  dispersion <- rep(dispersion, length.out = M)
  if (cobj@Rank != 1)
    stop("object must be a rank 1 model")

  dvecMax <- cbind(1, -0.5 * cobj@A / c(cobj@D), (cobj@A / c(2*cobj@D))^2)
  dvecTol <- cbind(0, 0, 1 / c(-2 * cobj@D)^1.5)
  dvecOpt <- cbind(0, -0.5 / c(cobj@D), 0.5 * cobj@A / c(cobj@D^2))

  if ((length(object@control$colx1.index) != 1) ||
     (names(object@control$colx1.index) != "(Intercept)"))
    stop("Can only handle noRRR=~1 models")

  okvals <- c(3*M, 2*M+1)
  if (all(length(coef(object)) != okvals))
    stop("Can only handle intercepts-only model with ",
         "eq.tolerances = FALSE")

  answer <- NULL
  Cov.unscaled <- array(NA, c(3, 3, M), dimnames = list(
      c("(Intercept)", "latvar", "latvar^2"),
      c("(Intercept)", "latvar", "latvar^2"), dimnames(cobj@D)[[3]]))
  for (spp in 1:M) {
    index <- c(M + ifelse(object@control$eq.tolerances, 1, M) + spp,
               spp,
               M + ifelse(object@control$eq.tolerances, 1, spp))
    vcov <- Cov.unscaled[,,spp] <-
        sobj@cov.unscaled[index, index]  # Order is A, D, B1
    se2Max <- dvecMax[spp,, drop = FALSE] %*% vcov %*% cbind(dvecMax[spp,])
    se2Tol <- dvecTol[spp,, drop = FALSE] %*% vcov %*% cbind(dvecTol[spp,])
    se2Opt <- dvecOpt[spp,, drop = FALSE] %*% vcov %*% cbind(dvecOpt[spp,])
    answer <- rbind(answer, dispersion[spp]^0.5 *
                   c(se2Opt = se2Opt, se2Tol = se2Tol, se2Max = se2Max))
  }

  link.function <- if (MaxScale == "predictors")
      remove.arg(object@misc$predictors.names[1]) else ""
  dimnames(answer) <- list(dimnames(cobj@D)[[3]],
                           c("Optimum", "Tolerance",
                             if (nchar(link.function))
                           paste(link.function, "(Maximum)", sep = "") else
                           "Maximum"))
  NAthere <- is.na(answer %*% rep(1, length.out = 3))
  answer[NAthere,] <- NA  # NA in tolerance means NA everywhere else
  new(Class = "vcov.qrrvglm",
      Cov.unscaled = Cov.unscaled,
      dispersion = dispersion,
      se = sqrt(answer))
}


setMethod("vcov", "rrvglm", function(object, ...)
    vcovrrvglm(object, ...))


setMethod("vcov", "qrrvglm", function(object, ...)
    vcovqrrvglm(object, ...))


setClass(Class = "vcov.qrrvglm", representation(
         Cov.unscaled = "array",  # permuted cov.unscaled
         dispersion = "numeric",
         se = "matrix"))



model.matrix.qrrvglm <- function(object,
                                 type = c("latvar", "vlm"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type, c("latvar", "vlm"))[1]

  switch(type,
         latvar  = Coef(object, ...)@latvar,
         vlm = object@x) 
}


setMethod("model.matrix",  "qrrvglm", function(object, ...)
           model.matrix.qrrvglm(object, ...))







perspqrrvglm <-
  function(x, varI.latvar = FALSE, refResponse = NULL,
      show.plot = TRUE,
      xlim = NULL, ylim = NULL,
      zlim = NULL,  # zlim ignored if Rank == 1
      gridlength = if (Rank == 1) 301 else c(51, 51),
      which.species = NULL,
      xlab = if (Rank == 1)
      "Latent Variable" else "Latent Variable 1",
      ylab = if (Rank == 1)
      "Expected Value" else "Latent Variable 2",
      zlab = "Expected value",
      labelSpecies = FALSE,  # For Rank == 1 only
      stretch = 1.05,  # quick and dirty, Rank == 1 only
      main = "",
      ticktype = "detailed", 
      col = if (Rank == 1) par()$col else "white",
      llty = par()$lty, llwd = par()$lwd,
      add1 = FALSE,
      ...) {
  oylim <- ylim
  object <- x  # Do not like x as the primary argument 
  coef.obj <- Coef(object, varI.latvar = varI.latvar,
                   refResponse = refResponse)
  if ((Rank <- coef.obj@Rank) > 2)
    stop("object must be a rank-1 or rank-2 model")
  fv <- fitted(object)
  NOS <- ncol(fv)  # Number of species
  M <- object@misc$M # 

  xlim <- rep(if (length(xlim)) xlim else
              range(coef.obj@latvar[, 1]), length = 2)
  if (!length(oylim)) {
    ylim <- if (Rank == 1)
              c(0, max(fv) * stretch) else
              rep(range(coef.obj@latvar[, 2]), length = 2)
  }
  gridlength <- rep(gridlength, length = Rank)
  latvar1 <- seq(xlim[1], xlim[2], length = gridlength[1])
  if (Rank == 1) {
    m <- cbind(latvar1)
  } else {
    latvar2 <- seq(ylim[1], ylim[2], length = gridlength[2])
    m <- expand.grid(latvar1,latvar2)
  }

  if (dim(coef.obj@B1)[1] != 1 ||
      dimnames(coef.obj@B1)[[1]] != "(Intercept)")
    stop("noRRR = ~ 1 is needed")
  LP <- coef.obj@A %*% t(cbind(m))  # M by n
  LP <- LP + c(coef.obj@B1)  # Assumes \bix_1 = 1 (intercept only)

  mm <- as.matrix(m)
  N <- ncol(LP)
  for (jay in 1:M) {
    for (ii in 1:N) {
      LP[jay, ii] <- LP[jay, ii] +
                     mm[ii, , drop = FALSE] %*%
                     coef.obj@D[,,jay] %*%
                     t(mm[ii, , drop = FALSE])
    }
  }
  LP <- t(LP)  # n by M


    fitvals <- object@family@linkinv(LP)  # n by NOS
    dimnames(fitvals) <- list(NULL, dimnames(fv)[[2]])
    sppNames <- dimnames(object@y)[[2]]
    if (!length(which.species)) {
      which.species <- sppNames[1:NOS]
      which.species.numer <- 1:NOS
    } else
    if (is.numeric(which.species)) {
      which.species.numer <- which.species
      which.species <- sppNames[which.species.numer]  # Convert to character
    } else {
      which.species.numer <- match(which.species, sppNames)
    }
    if (Rank == 1) {
      if (show.plot) {
        if (!length(oylim))
          ylim <- c(0, max(fitvals[,which.species.numer]) *
                       stretch)  # A revision
        col <- rep(col, length.out = length(which.species.numer))
        llty <- rep(llty, leng = length(which.species.numer))
        llwd <- rep(llwd, leng = length(which.species.numer))
        if (!add1)
          matplot(latvar1, fitvals, xlab = xlab, ylab = ylab,
                  type = "n",
                  main = main, xlim = xlim, ylim = ylim, ...) 
        for (jloc in 1:length(which.species.numer)) {
          ptr2 <- which.species.numer[jloc]  # points to species column
          lines(latvar1, fitvals[, ptr2],
                col = col[jloc],
                lty = llty[jloc],
                lwd = llwd[jloc], ...)
          if (labelSpecies) {
            ptr1 <- (1:nrow(fitvals))[max(fitvals[, ptr2]) ==
                                              fitvals[, ptr2]]
            ptr1 <- ptr1[1]
            text(latvar1[ptr1],
                 fitvals[ptr1, ptr2] + (stretch-1) * diff(range(ylim)),
                 label = sppNames[jloc], col = col[jloc], ...)
          }
        }
      }
    } else {
      max.fitted <- matrix(fitvals[, which.species[1]],
                           length(latvar1), length(latvar2))
      if (length(which.species) > 1)
      for (jlocal in which.species[-1]) {
        max.fitted <- pmax(max.fitted,
                           matrix(fitvals[, jlocal],
                                  length(latvar1), length(latvar2)))
      }
      if (!length(zlim))
        zlim <- range(max.fitted, na.rm = TRUE)


    perspdefault <- getS3method("persp", "default")
      if (show.plot)
        perspdefault(latvar1, latvar2, max.fitted,
                     zlim = zlim,
                     xlab = xlab, ylab = ylab, zlab = zlab,
                     ticktype = ticktype, col = col, main = main, ...)
    }

    invisible(list(fitted       = fitvals,
                   latvar1grid  = latvar1,
                   latvar2grid  = if (Rank == 2) latvar2 else NULL,
                   max.fitted   = if (Rank == 2) max.fitted else NULL))
}



 if (!isGeneric("persp"))
   setGeneric("persp", function(x, ...) standardGeneric("persp"),
              package = "VGAM")

setMethod("persp", "qrrvglm",
  function(x, ...) perspqrrvglm(x = x, ...))









Rank.rrvglm <- function(object, ...) {
  object@control$Rank
}


Rank.qrrvglm <- function(object, ...) {
  object@control$Rank
}


Rank.rrvgam <- function(object, ...) {
  object@control$Rank
}




concoef.qrrvglm <- function(object, varI.latvar = FALSE,
                          refResponse = NULL, ...) {
  Coef(object, varI.latvar = varI.latvar, refResponse = refResponse, ...)@C
}


concoef.Coef.qrrvglm <- function(object, ...) {
  if (length(list(...)))
    warning("Too late! Ignoring the extra arguments")
  object@C
}


latvar.rrvglm <- function(object, ...) {
  ans <- lvplot(object, show.plot = FALSE)
  if (ncol(ans) == 1)
    dimnames(ans) <- list(dimnames(ans)[[1]], "lv")
  ans
}



latvar.qrrvglm <- function(object,
                           varI.latvar = FALSE,
                           refResponse = NULL, ...) {
  Coef(object,
       varI.latvar = varI.latvar,
       refResponse = refResponse, ...)@latvar
}


latvar.Coef.qrrvglm <- function(object, ...) {
  if (length(list(...)))
    warning("Too late! Ignoring the extra arguments")
  object@latvar
}


Max.qrrvglm <-
  function(object, varI.latvar = FALSE,
           refResponse = NULL, ...) {
  Coef(object, varI.latvar = varI.latvar,
       refResponse = refResponse, ...)@Maximum
}


Max.Coef.qrrvglm <- function(object, ...) {
  if (length(list(...)))
    warning("Too late! Ignoring the extra arguments")
  if (any(slotNames(object) == "Maximum"))
    object@Maximum else
    Max(object, ...)
}


Opt.qrrvglm <-
  function(object, varI.latvar = FALSE, refResponse = NULL, ...) {
      Coef(object, varI.latvar = varI.latvar,
           refResponse = refResponse, ...)@Optimum
}


Opt.Coef.qrrvglm <- function(object, ...) {
  if (length(list(...)))
    warning("Too late! Ignoring the extra arguments")
  object@Optimum
}


Tol.qrrvglm <-
  function(object, varI.latvar = FALSE, refResponse = NULL, ...) {
      Coef(object, varI.latvar = varI.latvar,
           refResponse = refResponse, ...)@Tolerance
}


Tol.Coef.qrrvglm <- function(object, ...) {
  if (length(list(...)))
    warning("Too late! Ignoring the extra arguments")
  if (any(slotNames(object) == "Tolerance"))
   object@Tolerance else Tol(object, ...)
}



 if (FALSE) {
 if (!isGeneric("ccoef"))
    setGeneric("ccoef", function(object, ...) {
    .Deprecated("concoef")
    
    standardGeneric("ccoef")
    }) 

setMethod("ccoef",  "rrvglm",
  function(object, ...) concoef.qrrvglm(object, ...))
setMethod("ccoef", "qrrvglm",
  function(object, ...) concoef.qrrvglm(object, ...))
setMethod("ccoef",  "Coef.rrvglm",
  function(object, ...) concoef.Coef.qrrvglm(object, ...))
setMethod("ccoef", "Coef.qrrvglm",
  function(object, ...) concoef.Coef.qrrvglm(object, ...))
}


 if (!isGeneric("concoef"))
    setGeneric("concoef", function(object, ...)
    standardGeneric("concoef")) 

setMethod("concoef",  "rrvglm",
  function(object, ...) concoef.qrrvglm(object, ...))
setMethod("concoef", "qrrvglm",
  function(object, ...) concoef.qrrvglm(object, ...))
setMethod("concoef",  "Coef.rrvglm",
  function(object, ...) concoef.Coef.qrrvglm(object, ...))
setMethod("concoef", "Coef.qrrvglm",
  function(object, ...) concoef.Coef.qrrvglm(object, ...))








setMethod("coef", "qrrvglm",
  function(object, ...) Coef.qrrvglm(object, ...))
setMethod("coefficients", "qrrvglm",
  function(object, ...) Coef.qrrvglm(object, ...))


 if (!isGeneric("lv"))
    setGeneric("lv",
  function(object, ...) {
    .Deprecated("latvar")
    
    standardGeneric("lv")
  })


setMethod("lv",  "rrvglm",
  function(object, ...) {
    
    latvar.rrvglm(object, ...)
  })
setMethod("lv", "qrrvglm",
  function(object, ...) {
    
    latvar.qrrvglm(object, ...)
  })
setMethod("lv",  "Coef.rrvglm",
  function(object, ...) {
    
    latvar.Coef.qrrvglm(object, ...)
  })
setMethod("lv", "Coef.qrrvglm",
  function(object, ...) {
    
    latvar.Coef.qrrvglm(object, ...)
  })


 if (!isGeneric("latvar"))
     setGeneric("latvar",
  function(object, ...) standardGeneric("latvar")) 
setMethod("latvar",  "rrvglm",
  function(object, ...) latvar.rrvglm(object, ...))
setMethod("latvar", "qrrvglm",
  function(object, ...) latvar.qrrvglm(object, ...))
setMethod("latvar",  "Coef.rrvglm",
  function(object, ...) latvar.Coef.qrrvglm(object, ...))
setMethod("latvar", "Coef.qrrvglm",
  function(object, ...) latvar.Coef.qrrvglm(object, ...))


 if (!isGeneric("Max"))
    setGeneric("Max",
  function(object, ...) standardGeneric("Max")) 
setMethod("Max", "qrrvglm",
  function(object, ...) Max.qrrvglm(object, ...))
setMethod("Max", "Coef.qrrvglm",
  function(object, ...) Max.Coef.qrrvglm(object, ...))



setMethod("Max", "rrvgam",
  function(object, ...) Coef(object, ...)@Maximum)




 if (!isGeneric("Opt"))
    setGeneric("Opt",
  function(object, ...) standardGeneric("Opt"))
setMethod("Opt", "qrrvglm",
  function(object, ...) Opt.qrrvglm(object, ...))
setMethod("Opt", "Coef.qrrvglm",
  function(object, ...) Opt.Coef.qrrvglm(object, ...))


setMethod("Opt", "rrvgam",
  function(object, ...) Coef(object, ...)@Optimum)




 if (!isGeneric("Tol"))
    setGeneric("Tol",
  function(object, ...) standardGeneric("Tol")) 
setMethod("Tol", "qrrvglm",
  function(object, ...) Tol.qrrvglm(object, ...))
setMethod("Tol", "Coef.qrrvglm",
  function(object, ...) Tol.Coef.qrrvglm(object, ...))



cgo <- function(...) {
  stop("The function 'cgo' has been renamed 'cqo'. ",
       "Ouch! Sorry!")
}

clo <- function(...) {
  stop("Constrained linear ordination is fitted with ",
       "the function 'rrvglm'")
}




is.bell.vlm <-
is.bell.rrvglm <- function(object, ...) {
  M <- object@misc$M
  ynames <- object@misc$ynames
  ans <- rep(FALSE, length.out = M)
  if (length(ynames)) names(ans) <- ynames
  ans
}


is.bell.uqo <-
is.bell.qrrvglm <- function(object, ...) {
  is.finite(Max(object, ...))
}


is.bell.rrvgam <- function(object, ...) {
  NA * Max(object, ...)
}


 if (!isGeneric("is.bell"))
    setGeneric("is.bell",
  function(object, ...) standardGeneric("is.bell"))
setMethod("is.bell","qrrvglm",
  function(object,...) is.bell.qrrvglm(object,...))
setMethod("is.bell","rrvglm",
  function(object, ...) is.bell.rrvglm(object, ...))
setMethod("is.bell","vlm",
  function(object, ...) is.bell.vlm(object, ...))
setMethod("is.bell","rrvgam",
  function(object, ...) is.bell.rrvgam(object, ...))
setMethod("is.bell","Coef.qrrvglm",
  function(object,...) is.bell.qrrvglm(object,...))



 if (!isGeneric("Rank"))
    setGeneric("Rank",
  function(object, ...) standardGeneric("Rank"))

setMethod("Rank",  "rrvglm",
  function(object, ...) Rank.rrvglm(object, ...))
setMethod("Rank", "qrrvglm",
  function(object, ...) Rank.qrrvglm(object, ...))
setMethod("Rank", "rrvgam",
  function(object, ...) Rank.rrvgam(object, ...))







