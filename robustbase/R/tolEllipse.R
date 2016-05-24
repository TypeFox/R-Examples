#### This is from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/

##   I would like to thank Peter Filzmoser for providing the initial code of
##   this function.

tolEllipsePlot <-
    function(x, m.cov = covMcd(x), cutoff = NULL, id.n = NULL,
	     classic = FALSE, tol = 1e-07,
	     xlab = "", ylab = "", main = "Tolerance ellipse (97.5%)",
	     txt.leg = c("robust", "classical"),
	     col.leg = c("red", "blue"),
	     lty.leg = c("solid","dashed"))
{

##@bdescr
## Tolerance Ellipse Plot:
##    Plots the 97.5% tolerance ellipse of the bivariate data set (x).
##    The ellipse is defined by those data points whose distance (dist)
##    is equal to the squareroot of the 97.5% chisquare quantile with
##    2 degrees of freedom.

##@edescr
##
##@in  x	: [matrix] A data.frame or matrix, n > 2*p
##@in  m.cov	: [mcd object] An object of type mcd - its attributes
##				center and cov will be used
##@in  cutoff	: [number] Distance needed to flag data points outside the ellipse
##@in  outflag	: [logical] Whether to print the labels of the outliers
##@in  tol	: [number] tolerance to be used for computing the inverse see 'solve'.
##		    defaults to 1e-7

## MM: This is nothing else but a version  cluster::ellipsoidPoints() !! -- FIXME
    ellips <- function(loc, cov) {
	## calculates a 97,5% ellipsoid
	## input: data set, location and covariance estimate, cutoff

	dist <- sqrt(qchisq(0.975, 2))
	A <- solve(cov)
	eA <- eigen(A)
	ev <- eA$values
	lambda1 <- max(ev)
	lambda2 <- min(ev)
	eigvect <- eA$vectors[, order(ev)[2]]
	z <- seq(0, 2 * pi, by = 0.01)
	z1 <- dist/sqrt(lambda1) * cos(z)
	z2 <- dist/sqrt(lambda2) * sin(z)
	alfa <- atan(eigvect[2]/eigvect[1])
	r <- matrix(c(cos(alfa),  - sin(alfa), sin(alfa), cos(alfa)), ncol = 2)
	t(loc + t(cbind(z1, z2) %*% r))	#   xmin <- min(x, z[, 1])
    }

    ##	parameters and preconditions

    if(is.data.frame(x))
        x <- data.matrix(x)
    if(!is.matrix(x) || !is.numeric(x))
        stop("x is not a numeric dataframe or matrix.")

    n <- dim(x)[1]
    p <- dim(x)[2]

    if(p != 2)
	stop("Dimension {= ncol(x)} must be 2!")

    if(!is.numeric(m.cov$center) ||  !is.numeric(m.cov$cov))
	stop("argument 'm.cov' must have numeric components 'center' and 'cov'")

    x.loc <- m.cov$center
    x.cov <- n/(n - 1) * m.cov$cov
    xM <- colMeans(x)
    z1 <- ellips(loc = xM, cov = n/(n - 1) * cov.wt(x)$cov)
    z2 <- ellips(loc = x.loc, cov = x.cov)
    x1 <- c(min(x[, 1], z1[, 1], z2[, 1]), max(x[,1],z1[,1], z2[,1]))
    y1 <- c(min(x[, 2], z1[, 2], z2[, 2]), max(x[,2],z1[,2], z2[,2]))

    md <- sqrt(mahalanobis(x, xM, cov(x), tol=tol))
    rd <- sqrt(mahalanobis(x,m.cov$center, m.cov$cov, tol=tol))

    ## Note: the *calling* function may pass a 'missing' value
    if(missing(cutoff) || is.null(cutoff))
	cutoff <- sqrt(qchisq(0.975, df = 2))
    if(missing(id.n) || is.null(id.n))
	id.n <- sum(rd > cutoff)

    ### (2,1) is wrong for 'classic' -- we *overplot*:
    ## if(classic)
    ##  opr <- if(prod(par("mfrow"))== 1) par(mfrow=c(1,2), pty="m") else list()
    ## MM: this is *NOT* good :
    ## else par(mfrow = c(1, 1))

##  1. Robust tolerance
##  define the plot, plot a box, plot the "good" points,
##  plot the outliers either as points or as numbers depending on outflag,
##  plot the ellipse, write a title of the plot
    plot(x, xlim = x1, ylim = y1, xlab = xlab, ylab = ylab, main = main)
    box()
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    if(id.n >= 1) {
	ind <- sort(rd, index.return=TRUE)$ix[(n-id.n+1):n]
	text(x[ind, 1] + xrange/50, x[ind, 2], ind)
    }

    points(z2, type = "l", lty=lty.leg[1], col=col.leg[1])

##  2. Classical tolerance
    if(classic){
	points(z1, type = "l", lty=lty.leg[2], col=col.leg[2])
	legend("bottomright", txt.leg, lty = lty.leg, col = col.leg)

        ## par(opr)
    }

    invisible()
}
