AI <- function (y, X = NULL, ZETA = NULL, R = NULL, draw = TRUE, REML = TRUE, 
          silent = FALSE, iters = 50, constraint = TRUE, init = NULL, 
          sherman = FALSE, che = TRUE, MTG2 = FALSE, Fishers = FALSE, 
          gss = TRUE, forced = NULL) 
{
  y.or <- y
  x.or <- X
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  zig.zag <- function(zz) {
    prove <- diff(zz)
    res <- vector()
    for (i in 2:length(prove)) {
      if ((prove[i] < 0 & prove[i - 1] > 0) | (prove[i] > 
                                               0 & prove[i - 1] < 0)) {
        res[i] <- TRUE
      }
      else {
        res[i] <- FALSE
      }
    }
    res <- res[-1]
    probab <- length(which(res))/length(res)
    return(probab)
  }
  if (MTG2) {
    cat("MTG2 feature activated. Eigen decomposition of K will be performed")
  }
  if (che) {
    if (is.list(ZETA)) {
      if (is.list(ZETA[[1]])) {
        ZETA = ZETA
      }
      else {
        ZETA = list(ZETA)
      }
    }
    else {
      cat("\nThe random effects need to be provided in a list format, please see examples")
    }
  }
  if (is.null(X) & is.null(ZETA)) {
    tn = length(y)
    xm <- matrix(1, tn, 1)
    yv <- scale(y)
    res <- lm(yv ~ xm - 1)
  }
  else {
    if (is.null(X) & !is.null(ZETA)) {
      tn = length(y)
      xm <- matrix(1, tn, 1)
    }
    if (!is.null(X) & !is.null(ZETA)) {
      if (is.list(X)) {
        if (is.list(X[[1]])) {
          xm = X[[1]][[1]]
        }
        else {
          xm = X[[1]]
        }
      }
      else {
        xm = as.matrix(X)
      }
    }
    if (is.null(R)) {
      R <- diag(length(y))
    }
    if (che) {
      ZETA <- lapply(ZETA, function(x) {
        if (length(x) == 1) {
          provided <- names(x)
          if (provided == "Z") {
            y <- list(Z = x[[1]], K = diag(dim(x[[1]])[2]))
          }
          if (provided == "K") {
            y <- list(Z = diag(length(y)), K = x[[1]])
          }
          else {
            stop()
            cat("Names of matrices provided can only be 'Z' or 'K', the names you provided don't match the arguments required")
          }
        }
        else {
          y <- x
        }
        return(y)
      })
    }
    tokeep <- names(ZETA)
    df <- unlist(lapply(ZETA, function(x) {
      dim(x[[1]])[2]
    }))
    df2 <- sort(df, decreasing = FALSE)
    df.ord <- numeric()
    for (u in 1:length(df)) {
      df.ord[u] <- which(df2 %in% df[u])[1]
      df2[df.ord[u]] <- NA
    }
    ZETA <- ZETA[df.ord]
    names(ZETA) <- tokeep[df.ord]
    if (!is.null(init)) {
      init <- init[c(df.ord, (length(init)))]
    }
    if (!is.null(forced)) {
      forced <- forced[c(df.ord, (length(forced)))]
    }
    x.or <- as.matrix(xm)
    zeta.or <- ZETA
    zeta.or <- lapply(zeta.or, function(x) {
      lapply(x, as.matrix)
    })
    if (length(ZETA) == 1 & (dim(ZETA[[1]][[1]])[2] == dim(ZETA[[1]][[2]])[2])) {
      misso <- which(is.na(y))
      if (length(misso) > 0) {
        y[misso] <- median(y, na.rm = TRUE)
      }
    }
    ZETA2 <- ZETA
    y2 <- y
    good <- which(!is.na(y))
    if (length(ZETA) == 1 & MTG2 == TRUE & (dim(ZETA[[1]][[1]])[2] == 
                                            dim(ZETA[[1]][[2]])[2])) {
      ZETA <- lapply(ZETA2, function(x, good) {
        if (dim(x[[1]])[2] == dim(x[[2]])[2]) {
          x[[1]] <- x[[1]][good, good]
          x[[2]] <- x[[2]][good, good]
          return(x)
        }
        else {
          x[[1]] <- x[[1]][good, ]
          x[[2]] <- x[[2]]
          return(x)
        }
      }, good = good)
    }
    else {
      ZETA <- lapply(ZETA2, function(x, good) {
        x[[1]] <- x[[1]][good, ]
        x[[2]] <- x[[2]]
        return(x)
      }, good = good)
    }
    y <- y[good]
    ZETA <- lapply(ZETA, function(x) {
      lapply(x, as.matrix)
    })
    xm <- as.matrix(xm[good, ])
    txm <- t(xm)
    R <- R[good, good]
    if (length(ZETA) == 1 & MTG2 == TRUE & (dim(ZETA[[1]][[1]])[2] == 
                                            dim(ZETA[[1]][[2]])[2])) {
      EIGENS <- lapply(ZETA, function(x) {
        eigen(x[[2]])
      })
      Us <- lapply(EIGENS, function(x) {
        x$vectors
      })
      Usp <- as(do.call("adiag1", Us), Class = "sparseMatrix")
      Ds <- lapply(EIGENS, function(x) {
        diag(x$values)
      })
      Dsp <- as(do.call("adiag1", Ds), Class = "sparseMatrix")
      ZETA <- lapply(as.list(1:length(ZETA)), function(x, 
                                                       zz, kk) {
        list(Z = zz[[x]][[1]], K = kk[[x]])
      }, zz = ZETA, kk = Ds)
    }
    Zs <- lapply(ZETA, function(x) {
      x[[1]]
    })
    Gs <- lapply(ZETA, function(x) {
      x[[2]]
    })
    Zsp <- as(do.call("cbind", Zs), Class = "sparseMatrix")
    tZsp <- t(Zsp)
    Ksp <- as(do.call("adiag1", Gs), Class = "sparseMatrix")
    fail = FALSE
    ZETA2 <- lapply(ZETA, function(x) {
      y = list(Z = as(x[[1]], Class = "sparseMatrix"), 
               K = as(x[[2]], Class = "sparseMatrix"))
    })
    zvar <- which(unlist(lapply(ZETA, function(x) {
      names(x)[1]
    })) == "Z")
    om <- list()
    for (k in zvar) {
      om[[k]] <- tcrossprod(ZETA2[[k]][[1]], ZETA2[[k]][[1]] %*% 
                              (ZETA2[[k]][[2]]))
    }
    om[[length(om) + 1]] <- as(diag(length(y)), Class = "sparseMatrix")
    if (length(ZETA) == 1 & MTG2 == TRUE & (dim(ZETA[[1]][[1]])[2] == 
                                            dim(ZETA[[1]][[2]])[2])) {
      y <- as.vector(t(Usp) %*% as.matrix(y, ncol = 1))
      xm <- t(Usp) %*% xm
      txm <- t(xm)
    }
    var.y <- var(y, na.rm = TRUE)
    yv <- scale(y)
    nvarcom <- length(ZETA) + 1
    base.var <- var(yv, na.rm = TRUE)/nvarcom
    if (length(ZETA) == 1 & MTG2 == TRUE & (dim(ZETA[[1]][[1]])[2] == 
                                            dim(ZETA[[1]][[2]])[2])) {
      if (is.null(init)) {
        var.com <- c(rep(0.01, nvarcom))
      }
      else {
        var.com <- init/var.y
      }
    }
    else {
      if (is.null(init)) {
        var.com <- c(rep(0.01, nvarcom))
      }
      else {
        var.com <- init/var.y
      }
    }
    weird = FALSE
    tn = length(yv)
    logL2 = -1e+07
    logL2.stored <- round(logL2, 0)
    conv = 0
    wi = 0
    record <- matrix(var.com * var.y, ncol = 1)
    if (is.null(names(ZETA))) {
      varosss <- c(paste("u.", df.ord, sep = ""))
    }
    else {
      varosss <- c(names(ZETA))
    }
    lege2 <- list()
    for (k in 1:length(var.com)) {
      gh1 <- varosss[k]
      if (k == length(var.com)) {
        lege2[[k]] <- paste("Var(Error):")
      }
      else {
        lege2[[k]] <- paste("Var(", gh1, "):", sep = "")
      }
    }
    if (!silent) {
      count <- 0
      tot <- 15
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    ups <- numeric()
    if (is.null(forced)) {
      while (conv == 0) {
        wi = wi + 1
        if (!silent) {
          count <- count + 1
        }
        varo <- as.list(var.com)
        Gspo <- lapply(as.list(c(1:length(ZETA))), function(x, 
                                                            K, v) {
          oo = K[[x]] * as.numeric((v[[x]]))
          return(oo)
        }, K = Gs, v = varo)
        Gsp <- as(do.call("adiag1", Gspo), Class = "sparseMatrix")
        Rsp <- as(R * as.numeric(var.com[length(var.com)]), 
                  Class = "sparseMatrix")
        varo <- NULL
        Gspo <- NULL
        if (sherman) {
          Rinv = solve(Rsp, sparse = TRUE, tol = 1e-19)
          Ginv = solve(Gsp, sparse = TRUE, tol = 1e-19)
          ZRZG = solve(as(tZsp %*% Rinv %*% Zsp + Ginv, 
                          Class = "sparseMatrix"), sparse = TRUE, tol = 1e-19)
          vm <- Zsp %*% (Gsp %*% tZsp) + Rsp
          vmi = Rinv - (Rinv %*% Zsp %*% ZRZG %*% t(Zsp) %*% 
                          Rinv)
        }
        else {
          vm <- Zsp %*% crossprod(Gsp, tZsp) + Rsp
          vmi <- solve(vm, sparse = TRUE, tol = 1e-19)
        }
        vmi.xm <- vmi %*% xm
        xvx = txm %*% vmi.xm
        xvxi = solve(xvx, sparse = TRUE, tol = 1e-19)
        s2 = xvxi %*% txm %*% vmi
        pm = vmi - crossprod(t(vmi.xm), s2)
        vmi = NULL
        if (REML) {
          ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
          if (is.infinite(ddv)) {
            stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.")
          }
          else {
            logL = as.numeric(-0.5 * ((ddv) + determinant(xvx, 
                                                          logarithm = TRUE)$modulus[[1]] + t(yv) %*% 
                                        (pm %*% yv)))
          }
        }
        else {
          ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
          if (is.infinite(ddv)) {
            stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.")
          }
          else {
            logL = as.numeric(-0.5 * ((ddv) + t(yv) %*% 
                                        (pm %*% yv)))
          }
        }
        if (wi > 10 & gss == TRUE) {
          if ((zig.zag(logL2.stored[(length(logL2.stored) - 
                                     4):length(logL2.stored)]) == 1)) {
            wi = iters
            if (!silent) {
              setTxtProgressBar(pb, (tot/tot))
            }
            cat("\nA weird likelihood behavior has been encountered. \nBe careful with variance components returned.\nSystem has singularities, ML estimators returned.")
          }
        }
        if (abs(logL - logL2) < 0.001 | wi == iters) {
          conv = 1
          if (!silent) {
            setTxtProgressBar(pb, (tot/tot))
          }
          if (wi == iters) {
            MLE <- which(logL2.stored == max(logL2.stored))[1]
            last20 <- round(length(logL2.stored) * 0.2)
            best.try <- record[, (dim(record)[2] - last20):(dim(record)[2])]
            logL2.stored[which(logL2.stored == max(logL2.stored))]
            fail <- TRUE
          }
        }
        else {
          if (!silent) {
            setTxtProgressBar(pb, (count/tot))
          }
          logL2 = (logL)
          logL2.stored <- c(logL2.stored, logL2)
          aim = matrix(0, nvarcom, nvarcom)
          py = pm %*% yv
          for (i in 1:nvarcom) {
            for (j in 1:i) {
              aim[i, j] = as.numeric(0.5 * (t(yv) %*% 
                                              om[[i]] %*% pm %*% om[[j]] %*% pm %*% 
                                              py))
              if (i != j) {
                aim[j, i] <- aim[i, j]
              }
            }
          }
          if (rankMatrix(aim)[1] == dim(aim)[2]) {
            aimi = solve(aim)
          }
          else {
            aimi = ginv(aim)
          }
          dldv = matrix(0, nvarcom)
          for (k in 1:nvarcom) {
            prm = pm %*% om[[k]]
            tr1 = sum(diag(prm))
            dldv[k, 1] = -0.5 * tr1 + 0.5 * as.numeric(t(yv) %*% 
                                                         prm %*% py)
          }
          up = aimi %*% dldv
          var.com <- as.matrix(var.com) + up
          if (((abs(logL2.stored[length(logL2.stored)]) - 
                logL2.stored[length(logL2.stored) - 1]) > 
               0) & (wi > 15) & (gss == FALSE)) {
            non.zero <- which(as.vector(up) < -0.25)
            if (length(non.zero) > 0) {
              var.com[non.zero] <- base.var
            }
            LRes <- EM2(y = y, X = X, ETA = ZETA, R = R, 
                        iters = 5, REML = REML, draw = FALSE, silent = silent, 
                        init = var.com * var.y)
            var.com <- as.vector(LRes)/var.y
            weird = TRUE
          }
          ups <- cbind(ups, var.com)
          fail <- which(var.com <= 0)
          if (length(fail) > 0) {
            var.com[fail] <- 0.001
          }
          extreme <- which(var.com > 1.5)
          if (length(extreme) > 0 & wi > 1) {
            var.com[extreme] <- record[extreme, (wi - 
                                                   1)]/var.y
          }
          record <- cbind(record, var.com * as.numeric(var.y))
          if (draw) {
            ylim <- max(unlist(record), na.rm = TRUE)
            my.palette <- brewer.pal(7, "Accent")
            layout(matrix(1:2, 2, 1))
            plot(logL2.stored[-1], type = "l", main = "logLikelihood", 
                 col = my.palette[7], lwd = 3, las = 2, 
                 xaxt = "n", ylab = "logLikelihood value", 
                 xlab = "Iterations processed", cex.axis = 0.5)
            axis(1, las = 1, at = 0:10000, labels = 0:10000, 
                 cex.axis = 0.8)
            legend("bottomleft", legend = round(logL2, 
                                                3), bty = "n", cex = 0.7)
            plot(record[1, ], ylim = c(0, ylim), type = "l", 
                 las = 2, xaxt = "n", main = "Average Information algorithm results", 
                 col = my.palette[1], lwd = 3, ylab = "Value of the variance component", 
                 xlab = "Iterations processed", cex.axis = 0.5)
            axis(1, las = 1, at = 0:10000, labels = 0:10000, 
                 cex.axis = 0.8)
            for (t in 1:(dim(record)[1])) {
              lines(record[t, ], col = my.palette[t], 
                    lwd = 3)
            }
            ww <- dim(record)[1]
            lege <- list()
            for (k in 1:length(var.com)) {
              if (k == length(var.com)) {
                lege[[k]] <- paste("Var(e):", round(record[k, 
                                                           wi + 1], 4), sep = "")
              }
              else {
                lege[[k]] <- paste("Var(", varosss[k], 
                                   "):", round(record[k, wi + 1], 4), 
                                   sep = "")
              }
            }
            legend("topleft", bty = "n", cex = 0.7, col = my.palette, 
                   lty = 1, lwd = 3, legend = unlist(lege))
          }
          fail = FALSE
        }
      }
      if (fail) {
        var.com2 <- as.matrix(record[, MLE])
      }
      else {
        var.com2 <- as.matrix(record[, dim(record)[2]])
      }
      pushes <- apply(ups, 1, function(fg) {
        length(which(fg < 0))/length(fg)
      })
      nnn <- length(which(pushes > 0.5))
      if ((nnn >= 1 & nnn <= 4) & (nnn < (dim(var.com2)[1]) - 
                                   1) & (fail == FALSE) & (constraint == TRUE)) {
        cat("\nOne or more variance components close to zero.\nBoundary constraint applied\n")
        zero <- which(pushes > 0.5)
        nonzero <- (1:dim(var.com2)[1])[-zero]
        boost <- AI2(y = y.or, X = x.or, ZETA = zeta.or[-zero], 
                     R = NULL, REML = REML, draw = draw, silent = silent, 
                     iters = 20, init = as.vector(var.com2)[-zero], 
                     sherman = sherman)
        var.com2[nonzero, ] <- boost
        boost2 <- AI3(y = y.or, X = x.or, ZETA = zeta.or, 
                      R = NULL, REML = REML, draw = draw, forced = nonzero, 
                      silent = silent, iters = 15, init = as.vector(var.com2/var.y), 
                      sherman = sherman)
        var.com2[zero, ] <- boost2[zero, ]
      }
    }
    else {
      setTxtProgressBar(pb, (tot/tot))
      cat("\nVariance components forced\n")
      var.com2 <- as.matrix(forced, ncol = 1)
      logL <- 0
      if (is.null(names(ZETA))) {
        varosss <- c(paste("u.", df.ord, sep = ""))
      }
      else {
        varosss <- c(names(ZETA))
      }
    }
  }
  if (weird) {
    cat("\nWeird likelihood behavior found, Additional EM steps were taken to get better initial values\n")
  }
  AIC = as.vector((-2 * logL) + (2 * dim(xm)[2]))
  BIC = as.vector((-2 * logL) + (log(length(y)) * dim(xm)[2]))
  varo <- as.list(var.com2)
  Gspo <- lapply(as.list(c(1:length(ZETA))), function(x, K, 
                                                      v) {
    oo = K[[x]] * as.numeric((v[[x]]))
    return(oo)
  }, K = Gs, v = varo)
  Gsp <- as(do.call("adiag1", Gspo), Class = "sparseMatrix")
  Rsp <- as(R * as.numeric(var.com2[length(var.com2)]), Class = "sparseMatrix")
  varo = NULL
  Gspo = NULL
  if (sherman) {
    Rinv = solve(Rsp, sparse = TRUE, tol = 1e-19)
    Ginv = solve(Gsp, sparse = TRUE, tol = 1e-19)
    ZRZG = solve(tZsp %*% Rinv %*% Zsp + Ginv, sparse = TRUE, 
                 tol = 1e-19)
    vm <- Zsp %*% (Gsp %*% tZsp) + Rsp
    Vinv = Rinv - (Rinv %*% Zsp %*% ZRZG %*% t(Zsp) %*% Rinv)
  }
  else {
    vm <- Zsp %*% crossprod(Gsp, tZsp) + Rsp
    Vinv <- solve(vm, sparse = TRUE, tol = 1e-19)
  }
  xvx <- crossprod(xm, Vinv %*% xm)
  xvxi <- solve(xvx)
  beta <- xvxi %*% crossprod(xm, Vinv %*% y)
  Var.u <- vector(mode = "list", length = length(zvar))
  PEV.u <- Var.u
  u <- Var.u
  pm = Vinv - Vinv %*% xm %*% (xvxi %*% txm %*% (Vinv))
  ZETA3 <- lapply(ZETA, function(x) {
    y = list(Z = as(x[[1]], Class = "sparseMatrix"), K = as(x[[2]], 
                                                            Class = "sparseMatrix"))
  })
  for (h in zvar) {
    Var.u[[h]] <- (as.numeric(var.com2[h, 1])^2) * (crossprod(ZETA3[[h]][[1]] %*% 
                                                                ZETA3[[h]][[2]], pm) %*% (ZETA3[[h]][[1]] %*% ZETA3[[h]][[2]]))
    PEV.u[[h]] <- as.numeric(var.com2[h, 1]) * ZETA3[[h]][[2]] - 
      Var.u[[h]]
  }
  ee <- (y - (xm %*% beta))
  if (length(ZETA) == 1 & MTG2 == TRUE & (dim(ZETA[[1]][[1]])[2] == 
                                          dim(ZETA[[1]][[2]])[2])) {
    for (k in zvar) {
      u[[k]] <- solve(t(Us[[k]])) %*% (((ZETA3[[k]][[2]] * 
                                           as.numeric(var.com2[k, 1])) %*% t(ZETA3[[k]][[1]]) %*% 
                                          Vinv %*% ee))
    }
  }
  else {
    for (k in zvar) {
      u[[k]] <- ((ZETA3[[k]][[2]] * as.numeric(var.com2[k, 
                                                        1])) %*% t(ZETA3[[k]][[1]]) %*% Vinv %*% ee)
    }
  }
  u <- u[zvar]
  fitted.u <- 0
  for (h in 1:length(zeta.or)) {
    fitted.u <- fitted.u + (zeta.or[[h]][[1]] %*% u[[h]])
  }
  fitted.y <- (x.or %*% beta) + fitted.u
  fitted.y.good <- fitted.y[good]
  residuals3 <- y - fitted.y[good]
  rownames(beta) <- colnames(xm)
  for (i in 1:length(ZETA)) {
    rownames(u[[i]]) <- colnames(ZETA[[i]][[1]])
  }
  if (Fishers) {
    fishers = matrix(0, nvarcom, nvarcom)
    for (i in 1:nvarcom) {
      for (j in 1:i) {
        fishers[i, j] = as.numeric(0.5 * sum(diag((om[[i]] %*% 
                                                     pm %*% om[[j]] %*% pm))))
        if (i != j) {
          fishers[j, i] <- fishers[i, j]
        }
      }
    }
    fishers.inv <- solve(fishers)
  }
  else {
    fishers.inv = NULL
  }
  if (!is.null(names(ZETA))) {
    names(u) <- names(ZETA)
    names(Var.u) <- names(ZETA)
    names(PEV.u) <- names(ZETA)
  }
  logL <- as.vector(logL)
  out1 <- as.matrix(var.com2, ncol = 1)
  colnames(out1) <- "Variance Components"
  rownames(out1) <- c(paste("Var(", varosss, ")", sep = ""), 
                      "Var(Error)")
  res <- list(var.comp = out1, V.inv = Vinv, u.hat = u, Var.u.hat = Var.u, 
              PEV.u.hat = PEV.u, beta.hat = beta, Var.beta.hat = xvxi, 
              LL = logL, AIC = AIC, BIC = BIC, X = xm, fitted.y = fitted.y, 
              fitted.u = fitted.u, residuals = ee, cond.residuals = residuals3, 
              fitted.y.good = fitted.y.good, Z = Zsp, K = Ksp, fish.inv = fishers.inv)
  layout(matrix(1, 1, 1))
  return(res)
}