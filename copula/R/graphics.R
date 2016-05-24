## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

##' @title Check if function 'fun' is like pcopula or like pCopula
##' @param fun function such as pCopula, dcopula _or_ F.n()
##' @return logical: TRUE if "like pCopula" _or_ F.n()
##' @author Martin Maechler
chkFun <- function(fun) {
    if(!is.function(fun)) stop("the 'fun' argument is not even a function")
    isObj <- function(nm) any(nm == c("copula", "mvdc")) ## [pdq][Cc]opula
    nf <- names(formals(fun))
    if(isObj(nf[2]) || nf[1:2] == c("x","X")) ## || is F.n()
        TRUE
    else if(isObj(nf[1])) FALSE
    else NA # and the caller will produce an error eventually
}



perspCopula <- function(x, fun, n = 51, delta = 0,
                        xlab = "x", ylab = "y",
                        zlab = deparse(substitute(fun))[1],
                        theta = -30, phi = 30, expand = 0.618,
                        ticktype = "detail", ...)
{
  stopifnot(0 <= delta, delta < 1/2) ## delta <- .Machine$double.eps^(1/4)
  xis <- yis <- seq(0 + delta, 1 - delta, len = n)
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), n, n)
  T <- persp(xis, yis, zmat, xlab=xlab, ylab=ylab, zlab=zlab,
	     theta=theta, phi=phi, expand=expand, ticktype=ticktype, ...)
  invisible(list(x = xis, y = yis, z = zmat, persp = T))
}


contourCopula <- function(x, fun, n = 51, delta = 0, box01 = TRUE, ...) {
  stopifnot(0 <= delta, delta < 1/2) ## delta <- .Machine$double.eps^(1/4)
  xis <- yis <- seq(0 + delta, 1 - delta, len = n)
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), n, n)
  contour(xis, yis, zmat, ...)
  if(box01) rect(0,0,1,1, border="gray40", lty=3)
  invisible(list(x = xis, y = yis, z = zmat))
}

perspMvdc <- function(x, fun,
                      xlim, ylim, nx = 51, ny = 51,
                      xis = seq(xlim[1], xlim[2], length = nx),
                      yis = seq(ylim[1], ylim[2], length = ny),
                      xlab = "x", ylab = "y",
                      zlab = deparse(substitute(fun))[1],
                      theta = -30, phi = 30, expand = 0.618,
                      ticktype = "detail", ...)
{
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), nx, ny)
  T <- persp(xis, yis, zmat, xlab=xlab, ylab=ylab, zlab=zlab,
	     theta=theta, phi=phi, expand=expand, ticktype=ticktype, ...)
  invisible(list(x = xis, y = yis, z = zmat, persp = T))
}


contourMvdc <- function(x, fun, xlim, ylim, nx = 51, ny = 51,
                        xis = seq(xlim[1], xlim[2], length = nx),
                        yis = seq(ylim[1], ylim[2], length = ny),
                        box01 = FALSE, ...)
{
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), nx, ny)
  contour(xis, yis, zmat, ...)
  if(box01) rect(0,0,1,1, border="gray40", lty=3)
  invisible(list(x = xis, y = yis, z = zmat))
}

setMethod("persp",   signature("copula"), perspCopula)
setMethod("contour", signature("copula"), contourCopula)

## "F.n(), C.n()" -- once w ehave empirical
## setMethod("persp", signature("mvFn"),
##           function(x, ...) {
##               perspCopula(x, F.n, ...)
##               warning("persp(<mvFn>, ..)  implementation unfinished;
##  contact maintainer(\"copula\")") ## FIXME : use trans3d() to add data; *or*
##               ## ensure seeing vertical jumps, by using {x_i-delta, x_i+delta} or ..
##           })

setMethod("persp", signature("mvdc"), perspMvdc)
setMethod("contour", signature("mvdc"), contourMvdc)


### Graphical tools for detecting dependence
### Genest and Farve (2007, Journal of Hydrologic Engineering)

ChiPlot <- function(x, plot=TRUE, pval = 0.95, ...) {
#### originally proposed by Fisher and Switzer (1985 Biometrika, 2001 Am. Stat.)
  ## x is a n by 2 matrix
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  hfg <- sapply(1:n,
                function(i) {
                  H <- (sum(x[,1] <= x[i,1] & x[,2] <= x[i,2]) - 1) / (n - 1)
                  F <- (sum(x[,1] <= x[i,1]) - 1) / (n - 1)
                  G <- (sum(x[,2] <= x[i,2]) - 1) / (n - 1)
                  c(H, F, G)
                })
  H <- hfg[1,]
  F <- hfg[2,]
  G <- hfg[3,]
  chi <-(H - F * G) / sqrt(F * (1 - F) * G * (1 - G))
  lambda <- 4 * sign( (F - 0.5) * (G - 0.5) ) * pmax( (F - 0.5)^2, (G - 0.5)^2 )
  cp <- c(1.54, 1.78, 2.18)
  idx <- pmatch(pval, c(0.9, 0.95, 0.99))
  if (is.na(idx)) stop("pval must be one of 0.9, 0.95, 0.99.")
  cp <- cp[idx]
  ymax <- max(abs(na.omit(chi)), cp / sqrt(n))
  if (plot) {
    plot(lambda, chi, xlim=c(-1, 1), ylim=c(-ymax, ymax), ...)
    abline(0, 0, lty = 3, col="gray")
    abline(cp / sqrt(n), 0, lty = 3, col="blue")
    abline(- cp / sqrt(n), 0, lty = 3, col="blue")
    lines(c(0, 0), c(-2, 2), lty = 3, col="gray")
  }
  invisible(cbind(H, F, G, chi, lambda))
}

KPlot <- function(x, plot=TRUE, ...) {
#### Genest and Boies (2003, American Statistician)
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  H <- sapply(1:n,
              function(i) (sum(x[,1] <= x[i,1] & x[,2] <= x[i,2]) - 1) / (n - 1))
  H <- sort(H)
  K0 <- function(x) x - x * log(x)
  k0 <- function(x) - log(x)
  integrand <- function(w, i) w * k0(w) * K0(w)^(i - 1) * (1 - K0(w))^(n - i)
  W <- sapply(1:n,
              function(i) integrate(integrand, 0, 1, i = i,
                                    rel.tol=.Machine$double.eps^0.25)$value)
  W <- n * choose(n - 1, 1:n - 1) * W

  if (plot) {
    plot(W, H, xlim=c(0, 1), ylim=c(0, 1))
    curve(K0(x), add=TRUE, col="blue")
    abline(0, 1, col="gray")
  }
  invisible(cbind(H, W))
}

## x <- c(-2.224, -1.538, -0.807, 0.024, 0.052, 1.324)
## y <- c(0.431, 1.035, 0.586, 1.465, 1.115, -0.847)
## ChiPlot(cbind(x, y))
## KPlot(cbind(x, y))


### Enhanced splom #############################################################

##' @title A scatter plot matrix with nice variable names
##' @param data numeric matrix or as.matrix(.)able
##' @param varnames variable names, typically unspecified
##' @param Vname character string to become "root variable name"
##' @param col.mat matrix of colors
##' @param bg.col.mat matrix of background colors
##' @param ... further arguments to splom()
##' @return a splom() object
##' @author Martin Maechler
splom2 <- function(data, varnames=NULL, Vname="U", xlab="",
                   col.mat=NULL, bg.col.mat=NULL, ...)
{
    stopifnot(is.numeric(data <- as.matrix(data)),
	      (d <- ncol(data)) >= 1)
    if(is.null(varnames)) {
	varnames <- do.call(expression,
			    lapply(1:d, function(i)
				   substitute(italic(A[I]), list(A = as.name(Vname), I=0+i))))
    }
    n <- nrow(data)
    if(is.null(col.mat))
        col.mat <- matrix(trellis.par.get("plot.symbol")$col, n,d)
    if(is.null(bg.col.mat))
        bg.col.mat <- matrix(trellis.par.get("background")$col, n,d)
    ## From Deepayan Sarkar, working around missing feature
    ##		(which should be in next release) of lattice
    my.diag.panel <- function(x, varname, ...)
        diag.panel.splom(x, varname=parse(text=varname), ...)
    ## splom
    splom(~data[,1:d], varnames=varnames, diag.panel=my.diag.panel, xlab="",
          panel = function(x, y, i, j, ...) {
              panel.fill(bg.col.mat[i,j])
              panel.splom(x, y, col=col.mat[i,j], ...)
          }, ...)
}


### Enhanced, general Q-Q plot with rugs and confidence intervals ##############

##' Q-Q plot with rugs and pointwise asymptotic confidence intervals
##'
##' @title Q-Q plot with rugs and pointwise asymptotic confidence intervals
##' @param x data (n-vector)
##' @param qF theoretical quantile function
##' @param log character string indicating whether log-scale should be used; see
##'        ?plot.default
##' @param qqline.args argument list passed to qqline(); use NULL to omit Q-Q line
##' @param rug.args argument list passed to rug(); use NULL to omit rugs
##' @param alpha significance level
##' @param CI.args argument list passed to lines() for drawing confidence intervals;
##'        use NULL to omit confidence intervals
##' @param CI.mtext argument list for information about confidence intervals; use
##'        NULL to omit information about confidence intervals
##' @param main title (can be an expression; use "" for no title)
##' @param main.args argument list passed to mtext() for drawing title
##' @param xlab x axis label
##' @param ylab y axis label
##' @param file file name (with extension .pdf) or "" (no pdf)
##' @param width width parameter of pdf()
##' @param height height parameter of pdf()
##' @param ... additional arguments passed to plot()
##' @return Q-Q plot
##' @author Marius Hofert
##' Note: - used in Genest, Hofert, Neslehova (2013)
##'       - better than pointwise asymptotic CIs would be (non-parametric)
##'         bootstrapped ones
qqplot2 <- function(x, qF, log="", qqline.args=if(log=="x" || log=="y") list(untf=TRUE) else list(),
                    rug.args=list(tcl=-0.6*par("tcl")),
                    alpha=0.05, CI.args=list(col="gray50"),
                    CI.mtext=list(text=paste0("Pointwise asymptotic ", 100*(1-alpha),
                                  "% confidence intervals"), side=4,
                                  cex=0.6*par("cex.main"), adj=0, col="gray50"),
                    main=expression(bold(italic(F)~~"Q-Q plot")),
                    main.args=list(text=main, side=3, line=1.1,
                                   cex=par("cex.main"), font=par("font.main"),
                                   adj=par("adj"), xpd=NA),
                    xlab="Theoretical quantiles", ylab="Sample quantiles",
                    file="", width=6, height=6, ...)
{
    x. <- sort(x) # drops NA
    n <- length(x.)
    p <- ppoints(n)
    q <- qF(p)
    ## plot points
    doPDF <- nchar(file)
    if(doPDF) pdf(file=file, width=width, height=height)
    plot(q, x., xlab=xlab, ylab=ylab, log=log, ...) # empirical vs. theoretical quantiles
    if(nchar(main)) do.call(mtext, main.args)
    ## plot the line (overplots points, but that's good for the eye!)
    if(!is.null(qqline.args))
        if(nchar(log)==1 && (is.null(untf <- qqline.args$untf) || !untf))
            warning("for a Q-Q line in x-log- or y-log-scale, specify 'untf = TRUE' in qqline.args")
    ## draw the line (through the first and third quartile; see ?qqline)
    ## note: - abline(a=0, b=1) only true if data is standardized (mu=0, sig2=1)
    ##       - abline(..., untf=TRUE) displays a curve (proper line in log-space)
    ##         *unless* both axes are in log-scale
        else {
            if(log=="xy") do.call(qqline, args=c(list(y=log10(x.),
                                                 distribution=function(p) log10(qF(p))),
                                          qqline.args))
            else do.call(qqline, args=c(list(y=x., distribution=qF), qqline.args))
        }
    ## rugs
    if(!is.null(rug.args)) {
        do.call(rug, c(list(q, side=1), rug.args))
        do.call(rug, c(list(x., side=2), rug.args))
    }
    ## confidence intervals
    if(!is.null(CI.args)) {
        ## Pointwise approximate (CLT) two-sided 1-alpha confidence intervals
        ## (basic idea taken from fBasics)
        ## With up/low = p +- qnorm(1-a/2)*sqrt(p(1-p)/n) it follows that
        ##   IP(F^{-1}(p) in [F_n^{-1}(low), F_n^{-1}(up)]) (p ~> F(x))
        ## = IP(x in [F_n^{-1}(low), F_n^{-1}(up)])
        ## = IP(F_n(x) in [low, up]) ~= 1-a since (CLT)
        ##   sqrt(n)*(F_n(x)-p)/sqrt(p*(1-p)) ~= N(0,1) since
        ## X_i~F => F_n(x)=bar{Z}_n with Z_i=I_{X_i<=x}~B(1, p=F(x)) => mean=p, var=p(1-p)
        SE <- sqrt(p*(1-p)/n) # standard error
        qa2 <- qnorm(1-alpha/2)
        ## lower
        low <- p - qa2 * SE
        ilow <- 0 < low & low < 1 # in (0,1)
        low.y <- quantile(x., probs=low[ilow]) # F_n^{-1}(0<low<1)
        low.x <- q[ilow] # corresponding theoretical quantile at each ppoint
        ## upper
        up <- p + qa2 * SE
        iup <- 0 < up & up < 1 # in (0,1)
        up.y <- quantile(x., probs=up[iup]) # F_n^{-1}(0<up<1)
        up.x <- q[iup] # corresponding theoretical quantile at each ppoint
        ## draw lines
        do.call(lines, c(list(low.x, low.y), CI.args))
        do.call(lines, c(list(up.x, up.y), CI.args))
        ## info
        if(!is.null(CI.mtext)) do.call(mtext, CI.mtext)
    }
    if(doPDF) dev.off()
    invisible()
}
