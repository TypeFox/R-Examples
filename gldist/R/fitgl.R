#  Copyright (C) 2012 Yohan Chalabi
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 or 3 of
#  the License (at your option).
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

fitgl <- function(x, start, inc = FALSE, na.rm = FALSE,
                  method = c("mle", "hist", "prob", "quant", "shape"), ...) {

    method <- match.arg(method)
    if (identical(method, "shape")) inc <- FALSE

    if (na.rm) x <- na.omit(x)

    if (any(is.na(x)))
        stop("NAs are not allowed. You might want to use 'na.rm = TRUE'")

    med <- median(x)
    iqr <- IQR(x)

    # if the location and scale parameters are not included in the
    # optimization, we then scale the dataset to have med = 0 and iqr = 1.
    if (!inc)
        x <- (x - med) / iqr

    # extract additional arguments that should be passed to nlminb
    dots <- list(...)
    nm <- c('eval.max', 'iter.max', 'trace', 'abs.tol',
            'rel.tol', 'x.tol', 'step.min')
    control <- dots[names(dots) %in% nm]
    # and keep the other additional arguments for the objective function
    obj_control <- dots[!(names(dots) %in% nm)]

    obj <-
        switch(method,

               mle = {
                   function(par) {
                       if (any(is.na(par))) return(Inf)
                       if (!inc) par <- c(0, 1, par)
                       ans <-  -sum(log(dgl(x, par, maxit=1e4L)))
                       if (is.na(ans)) Inf else ans
                   }
               },

               hist = {
                   hist_args <- c(list(x = x, plot = FALSE), obj_control)
                   hh <- do.call(hist, hist_args)
                   function(par) {
                       if (any(is.na(par))) return(Inf)
                       if (!inc) par <- c(0, 1, par)
                       den <-  dgl(hh$mids, par, maxit=1e4L)
                       if (any(is.na(den))) Inf else mean((hh$density - den)^2)
                   }
               },

               prob = {
                   len <- length(x)
                   x <- sort(x)
                   p <- (1:len)/(len + 1)
                   function(par) {
                       if (any(is.na(par))) return(Inf)
                       if (!inc) par <- c(0, 1, par)
                       ans <- mean((p - pgl(x, par, maxit = 1e4L))^2)
                       if (is.nan(ans)) Inf else ans
                   }
               },

               quant = {
                   len <-
                       if (is.null(obj_control$len))
                           1000L
                       else
                           obj_control$len
                   p <- ppoints(len)
                   qs <- as.vector(quantile(x, p))
                   function(par) {
                       if (any(is.na(par))) return(Inf)
                       if (!inc) par <- c(0, 1, par)
                       q <- qgl(p, par)
                       if (any(is.na(q))) Inf else mean((q - qs)^2)
                   }
               },

               shape = {
                   S <- function(p, par) qgl(p, par)
                   glskew <- function(par) {
                       p <- (1:7)/8
                       (S(p[5], par) - S(p[3], par)) /
                           (S(p[7], par) - S(p[1], par))
                   }
                   glkurt <- function(par) {
                       p <- (1:7)/8
                       (S(p[7], par) - S(p[5], par) +
                        S(p[3], par) - S(p[1], par)) /
                            (S(p[6], par) - S(p[2], par))
                   }
                   q <- as.vector(quantile(x, (1:7)/8))
                   # med <- q[4]
                   # iqr <- q[6] - q[2]
                   skew <- (q[5] - q[3]) / (q[7] - q[1])
                   kurt <- (q[7] - q[5] + q[3] - q[1]) / (q[6] - q[2])
                   function(par) {
                       par <- c(0, 1, par)
                       #-> to ensure that fitted parameters are
                       #-> feasible for all points x
                       if (any(is.na(pgl(x, par, maxit = 1e4)))) return(Inf)
                       sample <- c(skew, kurt)
                       theoretical <- c(glskew(par), glkurt(par))
                       ans <- sum((sample - theoretical)^2)
                       if (is.na(ans)) Inf else ans
                   }
               })

    small <- 1e-4
    if (missing(start))
        start <- c(med, iqr , 0, .6)
    lower <- c(-Inf, small, -1 + small, small)
    upper <- c(Inf, Inf, 1 - small, 1 - small)

    ans <- nlminb(start[c(rep(inc, 2), TRUE, TRUE)],
                  obj,
                  lower = lower[c(rep(inc, 2), TRUE, TRUE)],
                  upper = upper[c(rep(inc, 2), TRUE, TRUE)],
                  control = control)

    # When location and scale where not included in the optimization,
    # add their sample estimates to the fitted shape parameters
    if (!inc)
        ans$par <- c(med, iqr, ans$par)

    names(ans$par) <- c("med", "iqr", "chi", "xi")
    ans

}
