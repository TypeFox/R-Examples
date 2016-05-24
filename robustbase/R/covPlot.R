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
## http://www.r-project.org/Licenses/

##   I would like to thank Peter Filzmoser for providing the initial code of
##   this function.

plot.mcd <-
    function(x,
             which=c("all", "dd","distance","qqchi2","tolEllipsePlot","screeplot"),
             classic= FALSE,
             ask = (which=="all" && dev.interactive()),
             cutoff = NULL, id.n, labels.id = rownames(x$X), cex.id = 0.75,
             label.pos = c(4,2), tol = 1e-7, ...)
{
    if (!inherits(x, "mcd"))
        stop("Use only with 'mcd' objects")
    covPlot(x$X, which= which, classic= classic, ask= ask, m.cov = x,
	    cutoff= cutoff, id.n = id.n, labels.id, cex.id = cex.id,
	    label.pos = label.pos, tol = tol, ...)
}

covPlot <-
    function(x, which = c("all", "dd", "distance", "qqchi2",
			  "tolEllipsePlot", "screeplot"),
	     classic = FALSE, ask = (which == "all" && dev.interactive()),
	     m.cov = covMcd(x), cutoff = NULL,
             id.n, labels.id = rownames(x), cex.id = 0.75,
             label.pos = c(4,2), tol = 1e-7, ...)
{

    ##@bdescr
    ##  Make plots based on the covariance structure of a data set:
    ##  dd       -  distance-distance plot: Robust distances versus
    ##              Mahalanobis distances
    ##  distance -  a plot of the robust distances
    ##  qqchi2   -  a qq-plot of the robust distances versus the
    ##              quantiles of the chi-squared distribution
    ##  tolEllipsePlot- a tolerance ellipse
    ##  screeplot-  a screeplot of the eigenvalues ov the covariance matrix
    ##
    ## Distance Plot:
    ## Draw a Distance-Distance Plot: Plots the robust distances
    ## versus the classical Mahalanobis distances as introduced by
    ## Rousseeuw, P. J., and van Zomeren, B. C. (1990). Unmasking
    ## Multivariate Outliers and Leverage Points. Journal of the American
    ## Statistical Association, 85, 633-639.
    ##
    ## The dashed line is the set of points where the robust distance is
    ## equal to the classical distance.
    ## The horizontal and vertical dotted lines are drawn at values equal cutoff
    ## which defaults to square root of the 97.5% quantile of a chi-squared
    ## distribution with p degrees of freedom. Points beyond these lines can
    ## be considered outliers.
    ##
    ##@edescr
    ##
    ##@in  x       : [matrix] A data.frame or matrix, n > 2*p
    ##@in  which   : [character] A plot option, one of:
    ##                classic: index plot of the classical mahalanobis distances
    ##                robust : index plot of the robust mahalanobis distances
    ##                dd :     distance-distance plot
    ##                index :  parallel index plot of classical and robust distances
    ##                all :    all three plots --- this is the default
    ##
    ##@in  classic : [logical] If true the classical plot will be displayed too
    ##                                   default is classic = FALSE
    ##@in  m.cov    : [list] An object like class "mcd" - only its attributes
    ##                                      center and cov will be used
    ##@in  cutoff  : [number] The cutoff value for the distances
    ##@in  id.n    : [number] number of observations to be identified with a label.
    ## 			Defaults to the number of observations with
    ##                  distance larger than cutoff -- missing is propagated
    ##@in  tol     : [number] tolerance to be used for computing the inverse
    ##                         - see 'solve'. defaults to 1e-7

    ## NOTE: The default tolerance 1e-7, will not work for some example
    ##       data sets, like milk or aircraft

    myscreeplot <- function(x, m.cov = covMcd(x))
    {
        erob   <- eigen(m.cov$cov,symmetric = TRUE, only.values = TRUE)$values
        eclass <- eigen(var(x), symmetric = TRUE, only.values = TRUE)$values

        leg.txt <- c("Robust", "Classical")
        leg.col <- c("green", "red")
        leg.pch  <- c(1,24)
        leg.lty  <- c("solid", "dotted")

        eall <- c(erob,eclass)
        ylim <- c( min(eall), max(eall))

        plot(erob, ylim=ylim, ylab="Eigenvalues", xlab="Index", type="n")
        legend("topright", leg.txt, pch = leg.pch, lty = leg.lty, col = leg.col)

        lines(erob,   type="o", pch= leg.pch[1], lty= leg.lty[1], col=leg.col[1])
        lines(eclass, type="o", pch= leg.pch[2], lty= leg.lty[2], col=leg.col[2])

        title(main = "Scree plot")
    }

    mydistplot <- function(x, cutoff, classic = FALSE, id.n) {
        ##  Index Plot:
        ##  Plot the vector x (robust or mahalanobis distances) against
        ##  the observation indexes. Identify by a label the id.n
        ##  observations with largest value of x. If id.n is not supplied,
        ##  calculate it as the number of observations larger than cutoff.
        ##  Use cutoff to draw a horisontal line.
        ##  Use classic = FALSE/TRUE to choose the label of the vertical axes

        n <- length(x)
	if(missing(id.n)) # maybe propagated
	    id.n <- length(which(x > cutoff))
        ylab <- paste("Square Root of",
                      if(classic) "Mahalanobis" else "Robust",
                      "distance")
	plot(x, type = "p", ylab = ylab, xlab = "Index",
	     main = "Distance Plot")
        label(1:n, x, id.n)
        abline(h = cutoff)
    }

    myddplot <- function(md, rd, cutoff, id.n) {
        ##  Distance-Distance Plot:
        ##  Plot the vector y = rd (robust distances) against
        ##  x = md (mahalanobis distances). Identify by a label the id.n
        ##  observations with largest rd. If id.n is not supplied, calculate
        ##  it as the number of observations larger than cutoff. Use cutoff
        ##  to draw a horisontal and a vertical line. Draw also a dotted line
        ##  with a slope 1.
        n <- length(md)
	if(missing(id.n)) # maybe propagated
	    id.n <- length(which(rd > cutoff))
        xlab <- "Mahalanobis distance"
        ylab <- "Robust distance"
        plot(md, rd, type = "p", xlab = xlab, ylab = ylab,
	     main = "Distance-Distance Plot")
        label(md, rd, id.n)
        abline(0, 1, lty = 2)
        abline(v = cutoff, h = cutoff)
    }

    qqplot <- function(x, p, cutoff = sqrt(qchisq(0.975, p)),
                       classic = FALSE, id.n)
    {
        ##  Chisquare QQ-Plot:
        ##  Plot the vector x (robust or mahalanobis distances) against
        ##  the square root of the quantiles of the chi-squared distribution
        ##  with p degrees of freedom.
        ##  Identify by a label the id.n observations with largest value of x.
        ##  If id.n is not supplied, calculate it as the number of observations
        ##  larger than cutoff.
        ##  Use classic = FALSE/TRUE to choose the label of the vertical axes

        ##  parameters and preconditions

        n <- length(x)
	if(missing(id.n)) # maybe propagated
	    id.n <- length(which(x > cutoff))
        qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

        x <- sort(x, index.return = TRUE)
        ix <- x$ix
        x <- x$x

	ylab <- paste(if(classic) "Mahalanobis" else "Robust", "distance")
	xlab <- "Square root of the quantiles of the chi-squared distribution"
	plot(qq, x, xlab = xlab, ylab = ylab, main = "Chisquare QQ-Plot")
	label(qq, x, id.n, ind = (n-id.n+1):n, labs = ix)
	abline(0, 1, lty = 2)
    } ## end{qqplot}

    label <- function(x, y, id.n,
		      ind = sort.list(y, decreasing = TRUE)[1:id.n],
		      labs = labels.id, adj.x = TRUE)
    {
	if(id.n > 0) { ## label the largest 'id.n' y-values
	    labpos <-
		if(adj.x) label.pos[1+ as.numeric(x > mean(range(x)))] else 3
	    text(x[ind], y[ind], labs[ind],
		 cex = cex.id, xpd = TRUE, pos = labpos, offset = 0.25)
	}
    }


    ## Begin{covPlot} -- arguments checking of preconditions

    if(is.data.frame(x))
        x <- data.matrix(x)
    if(!is.matrix(x) || !is.numeric(x))
        stop("x is not a numeric dataframe or matrix.")

    n <- dim(x)[1]
    p <- dim(x)[2]

    if(!is.numeric(m.cov$center) ||  !is.numeric(m.cov$cov))
	stop("argument 'm.cov' must have numeric components 'center' and 'cov'")
    if(length(m.cov$center) != p)
        stop("Data set and provided center have different dimensions!")

    ## ?covPlot says it only needs 'cov' and 'center'
    ## Maybe should be smarter and *test* for non-singularity
    if(is.numeric(m.cov$crit) && m.cov$crit == 0)
	stop( "The covariance matrix is singular!")

    if(is.null(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    ## now "more in line" with plot.lm()'s labeling:
    if(is.null(labels.id))
        labels.id <- as.character(1:n)
    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }
    which <- match.arg(which)

    md <- sqrt(mahalanobis(x, colMeans(x), var(x), tol = tol))
    rd <- sqrt(mahalanobis(x, m.cov$center, m.cov$cov, tol = tol))

    ## *Never* here : par(mfrow = c(1,1), pty = "m")
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    if(which == "all" || which == "distance") {
	if(classic) {
	    opr <- if(prod(par("mfrow")) == 1)
		par(mfrow = c(1,2), pty = "m") else list()
	}
        ## index plot of mahalanobis distances:
        mydistplot(rd, cutoff, id.n = id.n)
        if(classic) {
            ## index plot of robust distances:
            mydistplot(md, cutoff, classic = TRUE, id.n = id.n)
            par(opr)
        }
    }

    if(which == "all" || which == "dd") {
        myddplot(md, rd, cutoff = cutoff, id.n = id.n) # distance-distance plot
    }

    if(which == "all" || which == "qqchi2") {
        if(classic) {
            opr <- if(prod(par("mfrow")) == 1)
                par(mfrow = c(1,2), pty = "m") else list()
        }
        ## qq-plot of the robust distances versus the
        ## quantiles of the chi-squared distribution
        qqplot(rd, p, cutoff = cutoff, id.n = id.n)
        if(classic) { ## qq-plot of the mahalanobis distances

            qqplot(md, p, cutoff = cutoff, classic = TRUE, id.n = id.n)
            par(opr)
        }
    }

    if(which == "all" || which == "tolEllipsePlot") {
        if(p == 2)
            tolEllipsePlot(x, m.cov = m.cov, cutoff = cutoff,
                           id.n = id.n, classic = classic, tol = tol)
	else if(which != "all")
            warning("For tolerance ellipses the dimension 'p' must be 2!")
    }

    if(which == "all" || which == "screeplot") {
        myscreeplot(x, m.cov = m.cov)
    }

} ## end { covPlot }

## ddplot <- function(x,...) {
##     covPlot(x, which="dd", ...)
## }

## distplot <- function(x,...) {
##     covPlot(x, which="distance", ...)
## }

## chi2qqplot <- function(x,...) {
##     covPlot(x, which="qqchi2", ...)
## }

## ellipse() exists in other packages
## ellipse <- function(x,...) {
##     covPlot(x, which="tolEllipsePlot", ...)
## }
