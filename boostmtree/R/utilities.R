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


blup.solve <- function(transf.data, membership, sigma, Kmax) {
  lapply(1:Kmax, function(k) {
    pt.k <- (membership == k)
    XX <- Reduce("+", lapply(which(pt.k), function(j) {
      Xnew <- transf.data[[j]]$Xnew
      t(Xnew) %*% Xnew
    }))
    XY <- Reduce("+", lapply(which(pt.k), function(j) {
      Xnew <- transf.data[[j]]$Xnew
      Ynew <- transf.data[[j]]$Ynew
      t(Xnew) %*% Ynew
    }))
    XZ <- Reduce("+", lapply(which(pt.k), function(j) {
      Xnew <- transf.data[[j]]$Xnew
      Znew <- transf.data[[j]]$Znew
      t(Xnew) %*% Znew
    }))
    ZZ <- Reduce("+", lapply(which(pt.k), function(j) {
      Znew <- transf.data[[j]]$Znew
      t(Znew) %*% Znew
    }))
    ZY <- Reduce("+", lapply(which(pt.k), function(j) {
      Znew <- transf.data[[j]]$Znew
      Ynew <- transf.data[[j]]$Ynew
      t(Znew) %*% Ynew
    }))
    Q = ZZ + diag(sigma, nrow(ZZ))
    V = XZ %*% solve(Q, diag(1, nrow(ZZ)))
    A = XX - V %*% t(XZ)
    b = XY - V %*% ZY
    fix.eff <- tryCatch({qr.solve(A, b)}, error = function(ex){NULL})
    if (is.null(fix.eff)) {
      fix.eff <- rep(0, ncol(A))
    }
    rnd.eff <- tryCatch({qr.solve(Q, ZY - t(XZ) %*% fix.eff)}, error = function(ex){NULL})
    if (is.null(rnd.eff)) {
      rnd.eff <- rep(0, ncol(Q))
    }
    return(list(fix.eff = fix.eff, rnd.eff = rnd.eff))
  })
}
gls.mse  <- function(f, dta, trn, type = c("corCompSym", "corAR1", "corSymm", "iid")) {
  f <- as.formula(f)
  type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
    gls.grow <- tryCatch({
      if (type == "corCompSym") {
        gls(f, data = dta[trn, ], correlation = corCompSymm(form = ~ 1 | id))
      }
      else if (type == "corAR1") {
        gls(f, data = dta[trn, ], correlation = corAR1(form = ~ 1 | id))
      }
      else if (type == "corSymm") {
        gls(f, data = dta[trn, ], correlation = corSymm(form = ~ 1 | id))
      }
      else {
        gls(f, data = dta[trn, ])
      }
    }, error = function(ex){NULL})
  if (!is.null(gls.grow)) {
    gls.pred  <- tapply(model.matrix(f, dta[-trn,]) %*% gls.grow$coef,
                        dta[-trn, "id"], function(x) {x})
    y.test <- tapply(dta[-trn, "y"], dta[-trn, "id"], function(x) {x})
    l2Dist(gls.pred, y.test)
  }
  else {
    NA
  }
}
is.hidden.bootstrap <-  function (user.option) {
  if (is.null(user.option$bootstrap)) {
    "by.root"
  }
  else {
    as.character(user.option$bootstrap)
  }
}
is.hidden.CVlambda <-  function (user.option) {
  if (is.null(user.option$CVlambda)) {
    FALSE
  }
  else {
    user.option$CVlambda
  }
}
is.hidden.CVrho <-  function (user.option) {
  if (is.null(user.option$CVrho)) {
    TRUE
  }
  else {
    user.option$CVrho
  }
}
is.hidden.ntree <-  function (user.option) {
  if (is.null(user.option$ntree)) {
    1
  }
  else {
    max(1, user.option$ntree)
  }
}
is.hidden.partial <-  function (user.option) {
  if (is.null(user.option$partial)) {
    FALSE
  }
  else {
    TRUE
  }
}
is.hidden.rho <-  function (user.option) {
  if (is.null(user.option$rho)) {
    TRUE
  }
  else {
    user.option$rho
  }
}
l1Dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("y1 and y2 must have same the length\n")
  }
  mean(unlist(lapply(1:length(y1), function(i) {
    mean(abs(unlist(y1[[i]]) - unlist(y2[[i]])), na.rm = TRUE)
  })), na.rm = TRUE)
}
l2Dist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("y1 and y2 must have same the length\n")
  }
  sqrt(mean(unlist(lapply(1:length(y1), function(i) {
    mean((unlist(y1[[i]]) - unlist(y2[[i]]))^2, na.rm = TRUE)
  })), na.rm = TRUE))
}
line.plot <- function(x, y, ...) {
#  n <- length(x)
#  o <- lapply(1:n, function(i) {
#    lines(x[[i]], y[[i]], col = "gray", lty = 2)
#  })
  mapply(lines, x, y = y, col = "gray", lty = 2)
}
lowess.mod <- function(x, y, ...) {
  na.pt <- is.na(x) | is.na(y)
  if (all(na.pt) || sd(y, na.rm = TRUE) == 0) {
    return(list(x = x, y = y))
  }
  else {
    lowess(x[!na.pt], y[!na.pt], ...)
  }
}
parse.depth <- function(obj) {
  obj <- stat.split(obj)[[1]]
  depth <- unlist(lapply(1:length(obj), function(k) {
    if (!is.null(obj[[k]])) {
      min(obj[[k]][, "dpthID"], na.rm = TRUE)
    }
    else {
      NA
    }
  }))
  if (!all(is.na(depth))) {
    treeDepth <- max(unlist(lapply(1:length(obj), function(k) {
      if (!is.null(obj[[k]])) {
        max(obj[[k]][, "dpthID"], na.rm = TRUE)
      }
    })), na.rm = TRUE) 
  depth[is.na(depth)] <- treeDepth + 1
  }
  depth
}
penBS <- function(d, pen.ord = 2) {
  if (d >= (pen.ord + 1)) {
    D <- diag(d)
    for (k in 1:pen.ord) D <- diff(D)
    t(D) %*% D
  }
  else {
    diag(0, d)
  }
}
penBSderiv <- function(d, pen.ord = 2) {
  if (d > 0) {
    if (d >= (pen.ord + 1)) {
      pen.matx <- penBS(d, pen.ord)
      cbind(0, rbind(0, pen.matx))
    }
    else {
      warning("not enough degrees of freedom for differencing penalty matrix: setting penalty to zero\n")
      pen.matx <- diag(1, d + 1)
      pen.matx[1, 1] <- 0
      pen.matx
    }
  }
  else {
    0
  }
}
plot.profile.prx <- function(obj, col = NULL, rnd.case = NULL, cut = .95, restrictX = TRUE) {
  if (is.null(obj$proximity)) {
    stop("this function requires proximity = TRUE in the predict call")
  }
  prx <- obj$proximity
  time <- obj$boost.obj$time
  time.unq <- sort(unique(unlist(time)))
  DbetaT <- cbind(1, obj$boost.obj$D) %*% t(obj$boost.obj$beta)
  muGrid <- lapply(1:ncol(DbetaT), function(i) {DbetaT[, i]})
  mu <- obj$boost.obj$mu
  time.hat <- obj$time
  muhat <- obj$muhat
  if (is.null(rnd.case)) {
    rnd.case <- sample(1:nrow(prx), size = 1)
  }
  rnd.prx <- prx[rnd.case, ]
  prx.cut <- quantile(rnd.prx, cut)
  rnd.match <- which(rnd.prx >= prx.cut)
  rnd.prx <- rnd.prx[rnd.match]
  rnd.time <- time.hat[[rnd.case]]
  rnd.mean <- muhat[[rnd.case]]
  prx.mu <- c(do.call(cbind, lapply(rnd.match, function(i){muGrid[[i]]}))
                  %*% rnd.prx / sum(rnd.prx))
  prx.which.time <- is.element(time.unq, unlist(lapply(rnd.match, function(i){time[[i]]})))
  prx.time <- time.unq[prx.which.time]
  prx.mu <- prx.mu[prx.which.time] 
  if (is.null(col)) col <- rep(1, length(rnd.prx))
  plot(supsmu(rnd.time, rnd.mean),
       xlim = if (restrictX) range(rnd.time) else range(c(rnd.time, prx.time)),
       ylim = range(c(rnd.mean, prx.mu, unlist(lapply(rnd.match, function(i){mu[[i]]})))),
       type = "n", xlab = "time", ylab = "mean profile")
  for (i in rnd.match) {
    lines(lowess(time[[i]], mu[[i]]), lty = 2, col = col[i])
    #points(time[[i]], mu[[i]], pch = 16, lty = 2, col = col[i])
  }
  #lines(supsmu(rnd.time, rnd.mean) ,lty = 1, lwd = 2, col = 4)
  points(rnd.time, rnd.mean, pch = 16, cex = 0.25, col = 4)
  lines(lowess(rnd.time, rnd.mean) ,lty = 1, lwd = 2, col = 4)
  points(prx.time, prx.mu, pch = 16, cex = 0.25, col = 1)
  if (restrictX) {
    lines(supsmu(prx.time, prx.mu), lty = 1, lwd = 2, col = 1)
  }
  else {
    lines(lowess(prx.time, prx.mu), lty = 1, lwd = 2, col = 1)
  }
  legend("bottomleft", bty = "n", legend = c(paste("avg prx.", format(mean(rnd.prx, na.rm = TRUE), digits=3))))
  invisible(obj$boost.obj$x[rnd.match,, drop = FALSE])
}
point.plot <- function(x, y, ...) {
#  n <- length(x)
#  o <- lapply(1:n, function(i) {
#    points(x[[i]], y[[i]], pch = 16)
#  })
  mapply(points, x, y = y, pch = 16)
}
rho.inv <- function(ni, rho, tol = 1e-2) {
  m <- ni - 1
  if (m == 0) {
    0
  }
  else if (rho < 0 && abs(rho + 1 / m) <= tol) {
    (-1 / m + tol) / (m * tol)
  } 
  else {
    rho / (1 + m * rho)
  }
}
rho.inv.sqrt <- function(ni, rho, tol = 1e-2) {
  m <- ni - 1
  if (m == 0) {
    0
  }
  else {
    if (rho < 0 && abs(rho + 1 / m) <= tol) {
      rho <- -1 / m + tol
    }
    ri <- rho / (1 + m * rho)
    as.numeric(Re(polyroot(c(ri, -2, ni))))[1]
  }
}
sigma.robust <- function(lambda, rho) {
  lambda
}
