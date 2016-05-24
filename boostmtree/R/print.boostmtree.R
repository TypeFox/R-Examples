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


print.boostmtree <- function (x, ...)
{
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("boostmtree", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)' or '(boostmtree, predict)'")
  }
  if (sum(inherits(x, c("boostmtree", "grow"), TRUE) == c(1, 2)) == 2) {
    univariate <- length(x$id) == length(unique(x$id))
    cat("model                       :", class(x)[1], class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    if (x$ntree > 1) {
      cat("ntree                       :", x$ntree, "\n")
    }
    cat("number of K-terminal nodes  :", x$K, "\n")
    cat("regularization parameter    :", x$nu[1], "\n")
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!univariate) {
      cat("number of unique time points:", length(sort(unique(unlist(x$time)))), "\n")
      cat("avg. number of time points  :", round(mean(sapply(x$time, length), na.rm = TRUE), 2), "\n")
      cat("B-spline dimension          :", ncol(x$D), "\n")
      cat("penalization order          :", x$pen.ord, "\n")
    }
    else {
      cat("univariate family           :", TRUE, "\n")
    }
    cat("boosting iterations         :", x$M, "\n")
    if (!is.null(x$err.rate)) {
      cat("optimized number iterations :", x$Mopt, "\n")      
      if (!univariate) {
        cat("optimized rho               :", round(x$rho[x$Mopt], 4),  "\n")
        cat("optimized phi               :", round(x$phi[x$Mopt], 4),  "\n")
      }
      cat("in-sample cv RMSE           :", round(x$rmse, 4), "\n")
    }
  }
  else {
    univariate <- length(x$boost.obj$id) == length(unique(x$boost.obj$id))
    cat("model                       :", class(x)[1], class(x)[3], "\n")
    cat("fitting mode                :", class(x)[2], "\n")
    cat("sample size                 :", nrow(x$x), "\n")
    cat("number of variables         :", ncol(x$x), "\n")
    if (!univariate) {
      cat("number of unique time points:", length(sort(unique(unlist(x$time)))), "\n")
      cat("avg. number of time points  :", round(mean(sapply(x$time, length), na.rm = TRUE), 2), "\n")
      if (!is.null(x$err.rate)) {
        cat("optimized number iterations :", x$Mopt, "\n")
        cat("optimized rho               :", round(x$boost.obj$rho[x$Mopt], 4),  "\n")
        cat("optimized phi               :", round(x$boost.obj$phi[x$Mopt], 4),  "\n")
        cat("test set RMSE               :", round(x$err.rate[x$Mopt, "l2"], 4), "\n")
      }
    }
    else {
      if (!is.null(x$err.rate)) {
        cat("optimized number iterations :", x$Mopt, "\n")
        cat("test set RMSE               :", round(x$rmse, 4), "\n")
      }
      cat("univariate family           :", TRUE, "\n")
    }
  }
}
