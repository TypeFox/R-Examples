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


plot.boostmtree <- function (x, ...)
{
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {
    def.par <- par(no.readonly = TRUE) 
    n <- length(x$mu)
    M <- x$M
    univariate <- length(x$id) == length(unique(x$id))
    if (!univariate) {
      if (is.null(x$err.rate)) {
        layout(rbind(c(1, 4), c(2, 5), c(3, 6)), widths = c(1, 1))
      }
      else {
        layout(rbind(c(1, 3), c(2, 4), c(2, 5)), widths = c(1, 1))
      }
    }
    else {
      if (!is.null(x$err.rate)) {
        layout(rbind(c(1, 2)), widths = c(1, 1))
      }
    }
    if (!univariate) {
      plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
      line.plot(x$time, x$mu)
    }
    if (!univariate) {
      if (is.null(x$err.rate)) {
        plot(unlist(x$time), unlist(x$y) - unlist(x$mu), xlab = "time", ylab = "residual", type = "n")
        line.plot(x$time, lapply(1:n, function(i) {x$y[[i]] - x$mu[[i]]}))
        plot(unlist(x$y), unlist(x$mu), xlab = "y", ylab = "predicted", type = "n")
        line.plot(x$y, x$mu)
      }
      else {#error rate
        plot(1:M, x$err.rate[, "l2"],
           xlab = "iteration", 
           ylab = "In-sample estimated RMSE",
           type = "l", lty = 1)
        abline(v = x$Mopt, lty = 2, col = 2, lwd = 2)
      }
    }
    else {
      plot(unlist(x$y), unlist(x$mu), xlab = "y", ylab = "predicted", type = "n")
      point.plot(x$y, x$mu)
      abline(0, 1, col = "gray", lty = 2)
      if (!is.null(x$err.rate)) {
        plot(1:M, x$err.rate[, "l2"],
             xlab = "iteration", 
             ylab = "In-sample estimated RMSE",
             type = "l", lty = 1)
        abline(v = x$Mopt, lty = 2, col = 2, lwd = 2)
      }
    }
    if (!univariate) {
      plot(1:M, x$rho, ylim = range(lowess.mod(1:M, x$rho)$y),
           xlab = "iterations", ylab = expression(rho), type = "n")
      lines(lowess.mod(1:M, x$rho, f = 5/10))
      plot(1:M, x$phi, ylim = range(lowess.mod(1:M, x$phi)$y),
           xlab = "iterations", ylab = expression(phi), type = "n")
      lines(lowess.mod(1:M, x$phi, f = 5/10))
      plot(1:M, x$lambda, ylim = range(lowess.mod(1:M, x$lambda)$y),
           xlab = "iterations", ylab = expression(lambda), type = "n")
      lines(lowess.mod(1:M, x$lambda, f = 5/10))
    }
    par(def.par)
  }
  else {
    univariate <- length(x$boost.obj$id) == length(unique(x$boost.obj$id))
    if (!(univariate && is.null(x$err.rate))) {
      def.par <- par(no.readonly = TRUE) 
    }
    if (!univariate && is.null(x$err.rate)) {
      plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
      line.plot(x$time, x$mu)
    }
    else if (!is.null(x$err.rate)) {
      M <- x$boost.obj$M
      Mopt <- x$Mopt
      if (!univariate) {
        if (!is.null(x$vimp)) {
          layout(rbind(c(1, 3), c(1, 4), c(2, 5), c(2, 6)), widths = c(1, 1))
        }
        else {
          layout(rbind(c(1, 2), c(1, 3), c(1, 4), c(1, 5)), widths = c(1, 1))
        }
      }
      else {
        if (!is.null(x$vimp)) {
          layout(rbind(c(1, 2)), widths = c(1, 1))
        }
      }
      plot(1:M, x$err.rate[, "l2"],
           xlab = "iteration", 
           ylab = "RMSE prediction error",
           type = "l", lty = 1)
      abline(v = Mopt, lty = 2, col = 2, lwd = 2)
      if (!is.null(x$vimp)) {
        vimp <- 100 * (x$vimp / x$err.rate[Mopt, "l2"])
        barplot(vimp, las = 2, ylab = "vimp (%)", cex.names = 1.0)
      }
      if (!univariate) {
        plot(unlist(x$time), unlist(x$mu), xlab = "time", ylab = "predicted", type = "n")
        line.plot(x$time, x$mu)
        plot(1:M, x$boost.obj$rho, ylim = range(lowess.mod(1:M, x$boost.obj$rho)$y), 
             xlab = "iterations", ylab = expression(rho), type = "n")
        lines(lowess.mod(1:M, x$boost.obj$rho, f = 5/10))
        abline(v=Mopt, lty = 2, col = 2, lwd = 2)
        plot(1:M, x$boost.obj$phi, ylim = range(lowess.mod(1:M, x$boost.obj$phi)$y), 
             xlab = "iterations", ylab = expression(phi), type = "n")
        lines(lowess.mod(1:M, x$boost.obj$phi, f = 5/10))
        abline(v=Mopt, lty = 2, col = 2, lwd = 2)
        plot(1:M, x$boost.obj$lambda, ylim = range(lowess.mod(1:M, x$boost.obj$lambda)$y),
             xlab = "iterations", ylab = expression(lambda), type = "n")
        lines(lowess.mod(1:M, x$boost.obj$lambda, f = 5/10))
        abline(v=Mopt, lty = 2, col = 2, lwd = 2)
      }
    }
    if (!(univariate && is.null(x$err.rate))) {
      par(def.par)
    }
  }
}
