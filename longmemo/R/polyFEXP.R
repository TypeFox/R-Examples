####-*- mode: R; kept-old-versions: 12;	 kept-new-versions: 20; -*-

### ------ ==========================================
### 12.1.4 Approximate MLE for polynomial FEXP-models	-- book, p. 233-237
### ------ ==========================================

## Estimation of the parameters of a polynomial FEXP-model as defined
## in Beran (1993).
##
## The order of the polynomial is increased until no significant
## improvement is achieved.
## The data series may be divided into subseries for which the
## parameters are fitted separatly.
##____________________________________________________________________


## The function `per()'	 is defined in ./WhittleEst.R


### functions for plotting spectrum
### ===============================

## Below, we have
## Fourier frequencies	  fglim <- (1:nhalfm)*2*pi/n

llplot <- function(yper,spec) {
  ## Purpose: log-log plot of spectrum
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;	 modified: Martin Maechler, Date: Sep 95.
    n2 <- length(yper) ## == (n-1) %/% 2  ( = {Jan} trunc((n-1)/2) }
    if(length(spec) != n2) stop("'spec' and 'yper' must have same length")
    n <- 2*n2 + 1
    fglim <- (1:n2) * 2*pi / n
    plot (fglim,yper,log = "xy")
    lines(fglim,spec,col = 4)
}


lxplot <- function(yper,spec) {
    ## Purpose: log-x plot of spectrum
    ## -------------------------------------------------------------------------
    ## Arguments:
    ## -------------------------------------------------------------------------
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.
    n2 <- length(yper)
    if(length(spec) != n2) stop("'spec' and 'yper' must have same length")
    n <- 2*n2 + 1
    fglim <- (1:n2) * 2*pi / n
    plot (fglim,yper,log = "y")
    lines(fglim,spec,col = 4)
}


###MMMMMMMMMMMMMMMMMMMMM

### MAIN PROGRAM  (FEXP)

###MMMMMMMMMMMMMMMMMMMMM
FEXPest <- function(x, order.poly, pvalmax, verbose = FALSE)#, quiet = !verbose)
{
    ## Purpose: FEXP (Fractional EXP) longrange time-series model estimator
    ## ----------------------------------------------------------------------
    ## Arguments: x:
    ##		order.poly	 =JB= p	       maximal order of polynomial
    ##		pvalmax	 =JB= pvalmax  P-value for entering new polynomial terms
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 29 May 2004, 10:42
    ## ------  working from Jan Beran's original "main program"
    ##	Note that I've dropped the 'n.segments' loop over sub-parts of x[]

    cl <- match.call()
    n <- length(x)
    nhalfm <- (n-1) %/% 2
    ffr <- 2*pi/n* (1:nhalfm)		# Fourier frequencies
    xlong <- log(sqrt((1-cos(ffr))^2 + sin(ffr)^2)) # long-memory component
    ## = log |1 - e^{i x}|
    stopifnot((order.poly <- as.integer(order.poly)) >= 0)

    yper <- per(x)[2:(nhalfm+1)] # periodogram

    glimstep <- function(j) {

	## generalized regression

	frml <- yper ~ xlong
	if(j >= 1)
	    frml <-
		update(frml,
                       as.formula(paste("~ . +",
                                        if(j == 1) "ffr" else "poly(ffr, j)")))
	glmM <- glm(frml, family = Gamma(link = "log"))

	sRes <- summary(glmM, corr = FALSE)
	if(verbose) {
	    cat("glm() [",j,"]:\n",sep='')
            print(sRes)
	    cat("   P-val[eta]: ",
		paste(format.pval(sRes$coefficients[-1, "Pr(>|t|)"]),
		      collapse= ", "),
		"\n", sep='')
        }
        sRes$Fitted <- unname(glmM$fitted.values)
        sRes
    }

    ## loop for choosing polynomial
    for(j in 0:order.poly) {
	r <- glimstep(j)
        ##   --------
        early.stop <- j >= 1 && (maxPv <- max(r$coefficients[-1, "Pr(>|t|)"])) > pvalmax
        if(early.stop) ## condition for stopping
	    break ## poly found -- and forget "current"  r$theta etc
	## else
	r.last <- r
    }
    if(verbose && !early.stop) message("used full 'order.poly'")

    Cf <- r.last$coefficients
    dimnames(Cf)[[1]][2] <- "1 - 2*H"
    theta0 <- as.vector(Cf[, "Estimate"]) # drop names
    ## Std.Err(H) = sqrt(1/4 Var(th[2]) = 1/2 S.E.(th[2])  --> print() method
    structure(list(call = cl, n = n,
                   H = (1-theta0[2])/2,
		   coefficients = Cf, vcov = r.last$cov.scaled,## << FIXME: from R 2.14.0, use vcov(r.last)
		   order.poly = if(early.stop) j-1:1 else order.poly,
		   max.order.poly = order.poly, early.stop = early.stop,
		   maxPv = if(early.stop) maxPv, yper = yper,
		   ## spec = exp(c(cbind(1,xglim0) %*% theta0))),
		   ## spec = exp(c(r.last$X %*% theta0))),
		   spec = r.last$ Fitted),
	      class = "FEXP")
} ## FEXPest()


.ffreq <- function(n, full=FALSE) {
    stopifnot(is.numeric(n), n == (n <- as.integer(n)))
    n2 <- (if(full) n else n-1) %/% 2
    2 * pi/n * seq_len(n2)
}



## print method :
print.FEXP <-
    function (x, digits = getOption("digits"),
              ## signif.stars = getOption("show.signif.stars"),
              ...) {
    cat("'FEXP' estimator, call: ",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n",
        sprintf("polynomial order %d - %s; H = %s\n coefficients 'theta' =\n",
                as.integer(x$order.poly),
		if(x$early.stop) sprintf("selected by stopping early (P = %s > pvalmax)",
                                         format.pval(x$maxPv))
                else "user specified",
                formatC(x$H, digits= digits)))
    printCoefmat(x$coefficients, digits = digits,
		 signif.stars = FALSE, ## signif.stars = signif.stars,
		 na.print = "NA", ...)
    s.2H <- x$coefficients["1 - 2*H", "Std. Error"]
    Hround <- function(x) round(x, max(2, digits - 4))
    seForm <- function(s) sprintf(" (%g)", Hround(s))
    cat("  ==>            H =", Hround(x$H),     seForm(s.2H/2), "\n")
    cat(" <==> d := H - 1/2 =", Hround(x$H-0.5), seForm(s.2H/2), "\n\n")
    str(x[length(x)-(3:0)], no.list = TRUE)
    invisible(x)
}

plot.FEXP <-
    function (x, log = "xy", type = "l",
              col.spec = 4, lwd.spec = 2,
              xlab = NULL, ylab = expression(hat(f)(nu)),
              main = paste(deparse(x$call)[1]), sub = NULL, ...)
{
    n2 <- length(x$yper)
    n <- x$n
    fglim <- 2 * pi/n * seq_len(n2)
    if(is.null(xlab))
        xlab <- bquote(list(nu == 2*pi*j / n,
                            {} ~~ j %in% paste(1,ldots,.(n2)), n == .(n)))
    if(is.null(sub))
        sub <- paste("Data periodogram and fitted (FEXP) spectrum",
                     if(log == "xy") " [log - log]")

    plot(fglim, x$yper, log = log, type = type,
         xlab = xlab, ylab = ylab, main = main, sub = sub, ...)
    lines(fglim, x$spec, col = col.spec, lwd = lwd.spec)
    mtext(paste("used  order.poly =", x$order.poly),
          line = -1.2, col = col.spec)
}

### TODO: summary method
