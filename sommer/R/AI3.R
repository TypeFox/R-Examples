AI3 <-function (y, X = NULL, ZETA = NULL, R = NULL, draw = TRUE, REML = TRUE, 
          silent = FALSE, init = NULL, iters = 50, forced = NULL, sherman = FALSE) 
{
  y.or <- y
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
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
    x.or <- as.matrix(xm)
    zeta.or <- ZETA
    zeta.or <- lapply(zeta.or, function(x) {
      lapply(x, as.matrix)
    })
    ZETA2 <- ZETA
    y2 <- y
    good <- which(!is.na(y))
    ZETA <- lapply(ZETA2, function(x, good) {
      x[[1]] <- x[[1]][good, ]
      x[[2]] <- x[[2]]
      return(x)
    }, good = good)
    y <- y[good]
    ZETA <- lapply(ZETA, function(x) {
      lapply(x, as.matrix)
    })
    xm <- as.matrix(xm[good, ])
    R <- R[good, good]
    var.y <- var(y, na.rm = TRUE)
    yv <- scale(y)
    nvarcom <- length(ZETA) + 1
    base.var <- var(yv, na.rm = TRUE)/nvarcom
    if (is.null(init)) {
      var.com <- rep(base.var, nvarcom)
    }
    else {
      var.com <- init/var.y
    }
    zvar <- which(unlist(lapply(ZETA, function(x) {
      names(x)[1]
    })) == "Z")
    tn = length(yv)
    logL2 = -1e+07
    logL2.stored <- round(logL2, 0)
    conv = 0
    wi = 0
    record <- matrix(var.com * var.y, ncol = 1)
    Zs <- lapply(ZETA, function(x) {
      x[[1]]
    })
    Gs <- lapply(ZETA, function(x) {
      x[[2]]
    })
    Zsp <- as(do.call("cbind", Zs), Class = "sparseMatrix")
    fail = FALSE
    ZETA2 <- lapply(ZETA, function(x) {
      y = list(Z = as(x[[1]], Class = "sparseMatrix"), 
               K = as(x[[2]], Class = "sparseMatrix"))
    })
    om <- list()
    for (k in zvar) {
      om[[k]] <- tcrossprod(ZETA2[[k]][[1]], ZETA2[[k]][[1]] %*% 
                              (ZETA2[[k]][[2]]))
    }
    om[[length(om) + 1]] <- as(diag(length(y)), Class = "sparseMatrix")
    if (!silent) {
      count <- 0
      tot <- 15
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
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
      if (sherman) {
        Rinv = solve(Rsp)
        Ginv = solve(Gsp)
        ZRZG = solve(as(t(Zsp) %*% Rinv %*% Zsp + Ginv, 
                        Class = "sparseMatrix"))
        vm <- Zsp %*% Gsp %*% t(Zsp) + Rsp
        vmi = Rinv - (Rinv %*% Zsp %*% ZRZG %*% t(Zsp) %*% 
                        Rinv)
      }
      else {
        vm <- Zsp %*% Gsp %*% t(Zsp) + Rsp
        vmi <- solve(vm)
      }
      xvx = t(xm) %*% vmi %*% xm
      xvxi = solve(xvx)
      s1 = vmi %*% xm
      s2 = xvxi %*% t(xm) %*% vmi
      pm = vmi - s1 %*% s2
      if (REML == TRUE) {
        ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
        if (is.infinite(ddv)) {
          stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.")
        }
        else {
          logL = as.numeric(-0.5 * ((ddv) + log(det(xvx)) + 
                                      t(yv) %*% pm %*% yv))
        }
      }
      else {
        ddv <- determinant(vm, logarithm = TRUE)$modulus[[1]]
        if (is.infinite(ddv)) {
          stop("Infinite values found in the determinant, please make sure your variance-covariance matrices K's, are scaled matrices as regularly should be.")
        }
        else {
          logL = as.numeric(-0.5 * ((ddv) + t(yv) %*% 
                                      pm %*% yv))
        }
      }
      if (abs(logL - logL2) < 0.001 | wi == iters) {
        conv = 1
        if (!silent) {
          setTxtProgressBar(pb, (tot/tot))
        }
        if (wi == iters) {
          cat("\nMaximum number of iterations reached with no convergence. Changing to EM algorithm\n")
          LRes <- EM(y = y, X = X, ETA = ZETA, R = R, 
                     iters = iters, REML = REML, draw = draw, 
                     silent = silent)
          logL <- LRes$LL
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
            aim[i, j] = as.numeric(0.5 * (t(yv) %*% om[[i]] %*% 
                                            pm %*% om[[j]] %*% pm %*% py))
            if (i != j) {
              aim[j, i] <- aim[i, j]
            }
          }
        }
        aimi = solve(aim)
        dldv = matrix(0, nvarcom)
        for (k in 1:nvarcom) {
          prm = pm %*% om[[k]]
          tr1 = sum(diag(prm))
          dldv[k, 1] = -0.5 * tr1 + 0.5 * as.numeric(t(yv) %*% 
                                                       prm %*% py)
        }
        up = aimi %*% dldv
        var.com <- as.matrix(var.com) + up
        if (!is.null(forced)) {
          var.com[forced, ] <- init[forced]
        }
        fail <- which(var.com <= 0)
        if (length(fail) > 0) {
          var.com[fail] <- 0
        }
        extreme <- which(var.com > 1)
        if (length(extreme) > 0) {
          var.com[extreme] <- 0.002
        }
        record <- cbind(record, var.com * as.numeric(var.y))
        lege2 <- list()
        for (k in 1:length(var.com)) {
          if (k == length(var.com)) {
            lege2[[k]] <- paste("Var(e):")
          }
          else {
            lege2[[k]] <- paste("Var(u", k, "):", sep = "")
          }
        }
        if (draw) {
          ylim <- max(unlist(record), na.rm = TRUE)
          my.palette <- brewer.pal(7, "Accent")
          layout(matrix(1:2, 2, 1))
          plot(logL2.stored[-1], type = "l", main = "logLikelihood", 
               col = my.palette[7], lwd = 3, las = 2, xaxt = "n", 
               ylab = "logLikelihood value", xlab = "Iterations processed", 
               cex.axis = 0.5)
          axis(1, las = 1, at = 0:10000, labels = 0:10000, 
               cex.axis = 0.8)
          plot(record[1, ], ylim = c(0, ylim), type = "l", 
               las = 2, xaxt = "n", main = "Average Information algorithm results", 
               col = my.palette[1], lwd = 3, ylab = "Value of the variance component", 
               xlab = "Iterations processed", cex.axis = 0.5)
          axis(1, las = 1, at = 0:10000, labels = 0:10000, 
               cex.axis = 0.8)
          for (t in 1:(dim(record)[1])) {
            lines(record[t, ], col = my.palette[t], lwd = 3)
          }
          ww <- dim(record)[1]
          lege <- list()
          for (k in 1:length(var.com)) {
            if (k == length(var.com)) {
              lege[[k]] <- paste("Var(e):", round(record[k, 
                                                         wi + 1], 4), sep = "")
            }
            else {
              lege[[k]] <- paste("Var(u", k, "):", round(record[k, 
                                                                wi + 1], 4), sep = "")
            }
          }
          legend("topleft", bty = "n", cex = 0.7, col = my.palette, 
                 lty = 1, lwd = 3, legend = unlist(lege))
        }
        fail = FALSE
      }
    }
    if (fail) {
      var.com2 <- as.matrix(LRes$var.comp)
    }
    else {
      var.com2 <- as.matrix(record[, dim(record)[2]])
    }
  }
  return(var.com2)
}