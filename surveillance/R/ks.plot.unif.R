#################
# Plot the empirical distribution function of a sample from U(0,1)
# together with a confidence band of the corresponding K-S-test.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
# Parts of the code are taken from stats::ks.test, which has
# copyright 1995-2012 by The R Core Team under GPL-2 (or later).
# Furthermore, the C function calls are taken from
# http://svn.r-project.org/R/trunk/src/library/stats/src/ks.c (as at 2012-08-16),
# which similarly is Copyright (C) 1999-2009 by the R Core Team 
# and available under GPL-2. Somewhat disguised in their code is a reference
# that parts of their code uses code published in
# George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
# "Evaluating Kolmogorov's distribution".
# Journal of Statistical Software, Volume 8, 2003, Issue 18.
# URL: http://www.jstatsoft.org/v08/i18/.
#
#
# Parameters:
#    U - numeric vector containing the sample (NA's are silently removed)
#    conf.level - confindence level for the K-S-test,
#                 can also be a vector of multiple levels
#    exact - see ks.test
#    col.conf - colour of the confidence band
#    col.ref - colour of the reference line
#################

ks.plot.unif <- function (U, conf.level = 0.95, exact = NULL,
    col.conf = "gray", col.ref = "gray",
    xlab = expression(u[(i)]), ylab = "Cumulative distribution")
{
    stopifnot(is.vector(U, mode="numeric"))
    U <- U[!is.na(U)]
    n <- length(U)
    TIES <- FALSE
    if (anyDuplicated(U)) {
        warning("ties should not be present for the Kolmogorov-Smirnov test")
        TIES <- TRUE
    }
    if (is.null(exact)) exact <- (n < 100) && !TIES
    
    ## Helper function to invert the K-S test. The function
    ## pkolmogorov2x is the CDF of the Kolmogorov test statistic
    ## and is taken from the R project sources, which
    ## is (C) 1995-2009 by The R Core Team under GPL-2 
    f <- if (exact) {
        function (x, p) {               # x is the test statistic
            PVAL <- 1 - .C("pkolmogorov2x", p = as.double(x),
                           as.integer(n), PACKAGE = "surveillance")$p
            PVAL - p
        }
    } else {
        pkstwo <- function(x, tol = 1e-06) { # x is the test statistic
            ## stopifnot(length(x) == 1L)
            #Same copyright as above applies to the C code.
            if (is.na(x)) NA_real_ else if (x == 0) 0 else {
                .C("pkstwo", 1L, p = as.double(x), 
                  as.double(tol), PACKAGE = "surveillance")$p
            }
        }
        function (x, p) {
            PVAL <- 1 - pkstwo(sqrt(n) * x)
            PVAL - p
        }
    }
    
    ## Test inversion
    Dconf <- sapply(conf.level, function (level) {
        uniroot(f, lower=0, upper=1, p=1-level)$root
    })
    
    ## Small helper function to draw a line
    myabline <- function (a, b, x.grid = seq(0,1,length.out=101), ...) {
        lines(x.grid, a + b * x.grid, ...)
    }
    
    ## Figure 10 in Ogata (1988)
    plot(c(0,1), c(0,1), type="n", xlab=xlab, ylab=ylab)
    myabline(a=0, b=1, col=col.ref, lwd=2)
    rug(U)
    lines(ecdf(U), verticals=TRUE, do.points=FALSE)
    sapply(Dconf, function (D) {
    	myabline(a=D, b=1, col=col.conf, lty=2)
        myabline(a=-D, b=1, col=col.conf, lty=2)
    })
    #legend(x="topleft", col=col.conf, lty=2,
    #       legend=paste(100*conf.level,"% KS error bounds", sep=""))

    invisible()
}



######################################################################
# Check the residual process of fitted twinstim or twinSIR
# using ks.plot.unif on 1-exp(-diff(tau))
# and a scatterplot of u_i vs. u_{i+1} to inspect serial correlation
#
# Parameters:
#  object - a fitted twinSIR or twinstim model
#
# Draws the ECDF of the transformed residuals together with backtransformed
# 95% Kolmogorov-Smirnov error bounds.
######################################################################

checkResidualProcess <- function (object, plot = 1:2, mfrow = c(1,length(plot)),
                                  ...)
{
    stopifnot(inherits(object, c("twinSIR", "twinstim", "simEpidataCS")))
    
    ## check plot argument
    if (is.logical(plot)) plot <- which(rep(plot, length.out = 2)) else {
        stopifnot(is.vector(plot, mode="numeric"), plot %in% 1:2)
    }
    
    ## extract residual process
    tau <- do.call("residuals", args = list(substitute(object)),
                   envir = parent.frame())
    
    ## Transform to uniform variable
    Y <- diff(c(0,tau))
    U <- 1 - exp(-Y)
    
    ## Calculate KS test
    ks <- ks.test(U, "punif", alternative = "two.sided",
                  exact = match.call()[["exact"]])
    
    ## return value
    ret <- list(tau=tau, U=U, ks=ks)
    
    ## 2 types of plots
    plotcalls <- alist(
                 ## Investigate uniform distribution of U
                 ks.plot.unif(U, ...),
                 ## Investigate serial correlation between U_t and U_{t+1} which
                 ## corresponds to Figure 11 in Ogata (1988)
                 plot(tail(U,n=-1), head(U,n=-1),
                      xlab=expression(u[i]), ylab=expression(u[i+1]))
                 )
    
    ## eval selected plot calls
    if (length(plot) > 0L) {
        opar <- par(mfrow = mfrow); on.exit(par(opar))
        for (i in plot) eval(plotcalls[[i]])
        invisible(ret)
    } else {
        ret
    }
}
