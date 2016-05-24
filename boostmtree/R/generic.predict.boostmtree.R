####**********************************************************************
####**********************************************************************
####
####  BOOSTED MULTIVARIATE TREES FOR LONGITUDINAL DATA (BOOSTMTREE)
####  Version 1.1.0 (_PROJECT_BUILD_ID_)
####
####  Copyright 2016, University of Miami
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 3
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Project Partially Funded By:
####  ----------------------------------------------------------------
####  Dr. Ishwaran's work was funded in part by grant R01 CA163739 from
####  the National Cancer Institute.
####
####  Dr. Kogalur's work was funded in part by grant R01 CA163739 from 
####  the National Cancer Institute.
####  ----------------------------------------------------------------
####  Written by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Amol Pande
####    Division of Biostatistics
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  amoljpande@gmail.com
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Consultant Staff
####    Deptartment of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************


generic.predict.boostmtree <- function(object,
                                       x,
                                       tm,
                                       id,
                                       y,
                                       M,
                                       importance = TRUE,
                                       verbose = TRUE,
                                       eps = 1e-5,
                                       ...)
{
  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }
  if (verbose) cat("  running prediction mode for multivariate boosting\n")
  user.option <- match.call(expand.dots = TRUE)
  partial <- is.hidden.partial(user.option)
  if (!partial) {
    if (missing(x)) {
      X <- object$x
      n <- nrow(X)
      D <- object$D
      tm <- object$time
      tm.unq <- sort(unique(unlist(object$time)))
      testFlag <- FALSE
    }
    else {
      if (!missing(x) && (missing(id) || missing(tm))) {
        X <- x
        n <- nrow(X)
        tm.unq <- sort(unique(unlist(object$time)))
        tm <- lapply(1:n, function(i){tm.unq})
        id <- id.unq <- 1:n
        if (!missing(y)) {
          Y <- lapply(1:n, function(i) {y[i]})
          testFlag <- TRUE
        }
        else {
          testFlag <- FALSE
        }
      }
      else{
        if (missing(id)) {
          stop("test set id values are missing\n")
        }
        id.unq <- sort(unique(id))
        n <- length(id.unq)
        if (missing(x)) {
          stop("test set x values are missing\n")
        }
        X <- do.call(rbind, lapply(1:n, function(i) {
          x[id == id.unq[i],, drop = FALSE][1,, drop = FALSE]}))
        if (missing(tm)) {
          stop("test set time values are missing\n")
        }
        tm.unq <- sort(unique(tm))
        if (!missing(y)) {
          tm <- lapply(1:n, function(i) {tm[id == id.unq[i]]})
          Y <- lapply(1:n, function(i) {y[id == id.unq[i]]})
          testFlag <- TRUE
        }
        else {
          testFlag <- FALSE
        }
      }
    }
  }
  else {
      X <- x
      n <- nrow(X)
      tm.unq <- tm
      tm <- lapply(1:n, function(i){tm.unq})
      testFlag <- FALSE
  }
  if (object$d > 0) {
    if (length(tm.unq) > 1) {
      D <- cbind(1, bs(tm.unq, knots = attr(object$D, "knots"),
                       Boundary.knots = attr(object$D, "Boundary.knots"), degree = object$d))
    }
    else {
      stop("only one unique time point")
    }
  }
  else {
    D <- cbind(rep(1, length(tm.unq)))
  }
  if (missing(M)) {
    M <- object$M
    Mflag <- FALSE
  }
  else {
    M <- max(1, min(M, object$M))
    Mflag <- TRUE
  }
  K <- object$K
  nu <- object$nu
  ntree <- object$ntree
  p <- ncol(X)
  df.D <- ncol(D)
  xvar.names <- colnames(X)
  nu.vec <- c(nu[1], rep(nu[2], df.D - 1))
  Ymean <- object$ymean
  Ysd <- object$ysd
  gamma <- object$gamma
  baselearner <- object$baselearner
  beta <- matrix(0, n, df.D)
  if (ntree == 1) {
    beta.vimp <- beta.cov.vimp <- beta.time.vimp <- vector("list", p)
  }
  mu.list <- vector("list", M)
  forest.tol <- object$forest.tol
  vimpFlag <- testFlag && importance && ntree == 1
  vimp <- NULL
  if (ntree == 1) {
    rf.cores.old <- getOption("rf.cores")
    mc.cores.old <- getOption("mc.cores")
    membership <- mclapply(1:M, function(m) {
      options(rf.cores = 0, mc.cores = 1)
      c(predict.rfsrc(baselearner[[m]],
                      newdata = X,
                      membership = TRUE,
                      ptn.count = K,
                      importance = "none")$ptn.membership)
    })
    nullObj <- lapply(1:M, function(m) {
      orgMembership <- gamma[[m]][, 1]
      beta.m <- t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec
      Dbeta.m <- D %*% beta.m
      if (m == 1) {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      else {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          unlist(mu.list[[m-1]][i]) + Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      NULL
    })
    rm(nullObj)
    mu.list <- lapply(mu.list, function(mlist){  
      lapply(1:n,function(i) {mlist[[i]] * Ysd + Ymean})
    })
    if (testFlag) {
      err.rate <- matrix(unlist(lapply(mu.list, function(mlist) {
        c(l1Dist(Y, mlist), l2Dist(Y, mlist)) 
      })), ncol = 2, byrow = TRUE)
      colnames(err.rate) <- c("l1", "l2")
    }
    else {
      err.rate <- NULL
    }
    if (!Mflag && testFlag) {
      diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
      diff.err[is.na(diff.err)] <- 1
      if (sum(diff.err < Ysd * eps) > 0) {
        Mopt <- min(which(diff.err < eps))
      }
      else {
        Mopt <- M
      }
    }
    else {
      Mopt <- M
    }
    if (vimpFlag) {
      membershipNoise <- mclapply(1:Mopt, function(m) {
        Xnoise <- do.call(rbind, lapply(1:p, function(k) {
          X.k <- X
          X.k[, k] <- sample(X.k[, k])
          X.k
        }))
        options(rf.cores = 0, mc.cores = 1)
        c(predict.rfsrc(baselearner[[m]],
                        newdata = Xnoise,
                        membership = TRUE,
                        ptn.count = K,
                        importance = "none")$ptn.membership)
      })
      if (!is.null(rf.cores.old)) options(rf.cores = rf.cores.old)
      if (!is.null(mc.cores.old)) options(mc.cores = mc.cores.old)
      nullObj <- lapply(1:Mopt, function(m) {
        orgMembership <- gamma[[m]][, 1]
        beta.m <- t(t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd
        if (m == 1) {
          beta.m[, 1] <- beta.m[, 1] + Ymean
          beta <<- beta.m
        }
        else {
          beta <<- beta + beta.m
        }
        beta.vimp <<- lapply(1:p, function(k) {
          membership.k <- membershipNoise[[m]][((k-1) * n + 1):(k * n)]
          beta.m.k <- t(t(gamma[[m]][match(membership.k, orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd
          if (m == 1) {
            beta.m.k[, 1] <- beta.m.k[, 1] + Ymean
            beta.m.k
          }
          else {
            beta.vimp[[k]] + beta.m.k
          }
        })
        if (df.D > 1) {
          beta.cov.vimp <<- lapply(1:p, function(k) {
            b.c.v <- beta.vimp[[k]]
            b.c.v[, -1] <- beta[, -1]
            b.c.v
          })
          beta.time.vimp <<- lapply(1:p, function(k) {
            b.t.v <- beta.vimp[[k]]
            b.t.v[, 1] <- beta[, 1]
            b.t.v
          })
        }
        NULL
      })
      rm(nullObj)
    }
    else {
      nullObj <- lapply(1:Mopt, function(m) {
        orgMembership <- gamma[[m]][, 1]
        beta.m <- t(t(gamma[[m]][match(membership[[m]], orgMembership), -1, drop = FALSE]) * nu.vec) * Ysd
        if (m == 1) {
          beta.m[, 1] <- beta.m[, 1] + Ymean
          beta <<- beta.m
        }
        else {
          beta <<- beta + beta.m
        }
        NULL
      })
    }
  }
  else{
    nullObj <- lapply(1:M, function(m) {
      gm <- baselearner[[m]]$gm
      Xnew <- baselearner[[m]]$Xnew
      pen <- baselearner[[m]]$pen
      forest.wt <- predict.rfsrc(baselearner[[m]]$forest,
                                 newdata = X,
                                 importance = "none",
                                 forest.wt = TRUE)$forest.wt
      beta.m.org <- do.call("cbind", mclapply(1:n, function(i) {
        fwt.i <- forest.wt[i, ]
        fwt.i[fwt.i <= forest.tol] <- 0
        pt.i <- (fwt.i != 0)
        if (sum(pt.i) > 0) {
          fwt.i <- fwt.i / sum(fwt.i)
          YnewSum <- colSums(fwt.i[pt.i] * gm[pt.i,, drop = FALSE])
          XnewSum <- Reduce("+", lapply(which(pt.i), function(j) {fwt.i[j] * Xnew[[j]]}))
          XnewSum <- XnewSum + pen
          qr.obj <- tryCatch({qr.solve(XnewSum, YnewSum)}, error = function(ex){NULL})
          if (!is.null(qr.obj)) {
            qr.obj
          }
          else {
            rep(0, df.D)
          }
        }
        else {
          rep(0, df.D)
        }
      }))
      beta.m <- t(beta.m.org * nu.vec * Ysd)
      if (m == 1) {
        beta.m[, 1] <- beta.m[, 1] + Ymean
        beta <<- beta.m
      }
      else {
        beta <<- beta + beta.m 
      }
      Dbeta.m <- D %*% (beta.m.org * nu.vec)
      if (m == 1) {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      else {
        mu.list[[m]] <<- lapply(1:n, function(i) {
          unlist(mu.list[[m-1]][i]) + Dbeta.m[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
      }
      NULL
    })
    mu.list <- lapply(mu.list, function(mlist){  
      lapply(1:n,function(i) {mlist[[i]] * Ysd + Ymean})
    })
    if (testFlag) {
      err.rate <- matrix(unlist(lapply(mu.list, function(mlist) {
        c(l1Dist(Y, mlist), l2Dist(Y, mlist))
      })), ncol = 2, byrow = TRUE) 
      colnames(err.rate) <- c("l1", "l2")
    }
    else {
      err.rate <- NULL
    }
    if (!Mflag && testFlag) {
      diff.err <- abs(err.rate[, "l2"] - min(err.rate[, "l2"], na.rm = TRUE))
      diff.err[is.na(diff.err)] <- 1
      if (sum(diff.err < Ysd * eps) > 0) {
        Mopt <- min(which(diff.err < eps))
      }
      else {
        Mopt <- M
      }
    }
    else {
      Mopt <- M
    }
  }
  DbetaT <- D %*% t(beta)
  muhat <- lapply(1:n, function(i) {DbetaT[, i]})
  if (vimpFlag) {
    if (df.D <= 1) {
      vimp <- unlist(mclapply(1:p, function(k) {
        DbetaT.k <- D %*% t(beta.vimp[[k]])
        mu.k <- lapply(1:n, function(i) {
          DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
        l2Dist(Y, mu.k) - err.rate[Mopt, "l2"]
      }))
      names(vimp) <- xvar.names
    }
    if (df.D > 1) {
      vimp.cov <- unlist(mclapply(1:p, function(k) {
        DbetaT.k <- D %*% t(beta.cov.vimp[[k]])
        mu.k <- lapply(1:n, function(i) {
          DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
        l2Dist(Y, mu.k) - err.rate[Mopt, "l2"]
      }))
      vimp.cov.time <- unlist(mclapply(1:p, function(k) {
        DbetaT.k <- D %*% t(beta.time.vimp[[k]])
        mu.k <- lapply(1:n, function(i) {
          DbetaT.k[, i][match(tm[[i]], tm.unq, tm[[i]])]
        })
        l2Dist(Y, mu.k) - err.rate[Mopt, "l2"]
      }))
      mu.time <- sapply(1:n, function(i) {
        DbetaT[, i][sample(match(tm[[i]], tm.unq, tm[[i]]))]
      })
      vimp.time <- l2Dist(Y, mu.time) - err.rate[Mopt, "l2"]
      vimp <- c(vimp.cov, vimp.cov.time, vimp.time)
      names(vimp) <- c(xvar.names, paste(xvar.names, "time", sep=":"), "time")
    }
  }
  object$baselearner <- object$membership <- object$gamma <- NULL
  pobj <- list(
               boost.obj = object,
               x = X,
               time = tm,
               time.unq = tm.unq,
               y = if (testFlag) Y else NULL,
               mu = mu.list[[Mopt]],
               muhat = muhat,
               phi = object$phi[Mopt],
               rho = object$rho[Mopt],
               err.rate = if (!is.null(err.rate)) err.rate / Ysd else NULL,
               rmse = if (!is.null(err.rate)) err.rate[Mopt, "l2"] / Ysd else NULL,
               vimp = if (!is.null(vimp)) vimp / err.rate[Mopt, "l2"] else NULL,
               Mopt = if (testFlag) Mopt else NULL)
  class(pobj) <- c("boostmtree", "predict", class(object)[3])
  invisible(pobj)
}
