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


partialPlot <- function (obj,
                         xvar.names,
                         tm,
                         npts = 25,
                         subset,
                         plot.it = TRUE,
                         ...)
{
  if (sum(inherits(obj, c("boostmtree", "grow"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(boostmtree, grow)'")
  }
  if (missing(xvar.names)) {
    xvar.names <- colnames(obj$x)
  }
  xvar.names <- intersect(xvar.names, colnames(obj$x))
  if (length(xvar.names) == 0) {
    stop("x-variable names provided do not match original variable names")
  }
  n.xvar <- length(xvar.names)
  tmOrg <- sort(unique(unlist(obj$time)))
  if (missing(tm)) {
    tm.q <- unique(quantile(tmOrg, (1:9)/10, na.rm = TRUE))
    tm.pt <- sapply(tm.q, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  else {
    tm.pt <- sapply(tm, function(tt) {#assign original time values
      max(which.min(abs(tmOrg - tt)))
    })
  }
  n.tm <- length(tm.pt)
  if (!missing(subset)) {
    obj$x <- obj$x[subset,, drop = FALSE]
  }
  p.obj <- lapply(xvar.names, function(nm) {
    x <- obj$x[, nm]
    n.x <- length(unique(x))
    x.unq <- sort(unique(x))[unique(as.integer(seq(1, n.x, length = min(npts, n.x))))]
    newx <- obj$x
    rObj <- t(sapply(x.unq, function(xu) {
      newx[, nm] <- rep(xu, nrow(newx))
      mu <- predict(obj, x = newx, tm = tmOrg, partial = TRUE, ...)$mu
      mn.x <- colMeans(do.call(rbind, lapply(mu, function(mm) {mm[tm.pt]})))
      c(xu, mn.x)
    }))
    colnames(rObj) <- c("x", paste("y.", 1:length(tm.pt), sep = ""))
    rObj
  })
  names(p.obj) <- xvar.names
  l.obj <- lapply(p.obj, function(pp) {
    x <- pp[, 1]
    y <- apply(pp[, -1, drop = FALSE], 2, function(yy) {
      lowess(x, yy)$y})
    rObj <- cbind(x, y)
    colnames(rObj) <- c("x", paste("y.", 1:length(tm.pt), sep = ""))
    rObj
  })
  names(l.obj) <- xvar.names
  if (plot.it) {
    def.par <- par(no.readonly = TRUE) 
    for (k in 1:n.xvar) {
      plot(range(l.obj[[k]][, 1]), range(l.obj[[k]][, -1]), type = "n",
           xlab = xvar.names[k], ylab = "predicted y (adjusted)")
      for (l in 1:n.tm) {
        lines(l.obj[[k]][, 1], l.obj[[k]][, -1, drop = FALSE][, l], type = "l", col = 1)
      }
    }
    par(def.par)
  }
  return(invisible(list(p.obj = p.obj, l.obj = l.obj, time = tmOrg[tm.pt])))
}
