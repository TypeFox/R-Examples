####-*- mode: R; kept-old-versions: 12;  kept-new-versions: 20; -*-

### ------ =========================================
### 12.1.3 Whittle estimator for fractional Gaussian
### ------ noise and fractional ARIMA(p,d,q)		-- book, p. 223-233
###	   =========================================

## Splus-functions and program for the calculation of Whittle's
## estimator and the goodness of fit statistic defined in Beran (1992).
##
## The models are fractional Gaussian noise or fractional ARIMA.
##
## The data series may be divided
## into subseries for which the parameters are fitted separatly.
##____________________________________________________________________


CetaFGN <- function(eta, m = 10000, delta = 1e-9)
{
    ## Purpose: Covariance of hat{eta} for fGn = fractional Gaussian noise
    ## -------------------------------------------------------------------------
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

    if(1 > (M <- length(eta)))
        stop("eta[] must have length at least 1") #FIXME
    if(2 > (m <- as.integer(m)))
        stop(sQuote("m")," must be at least 2")
    mhalfm <- (m-1) %/% 2:2

    ## partial derivatives of log f (at each Fourier frequency)

    lf <- matrix(ncol = M,nrow = mhalfm)
    f0 <- specFGN(eta,m)$spec

    for(j in (1:M)) { ## FIXME{MM}: better not using same delta for all [j]
        etaj <- eta
        etaj[j] <- etaj[j]+delta
        fj <- specFGN(etaj,m)$spec
        lf[,j] <- log(fj/f0)/delta
    }

    solve(crossprod(lf)) *m ## << improve: chol2inv(qr(lf)) !

}## CetaFGN()


CetaARIMA <- function(eta, p,q, m = 10000, delta = 1e-9)
{
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

  if(1 > (M <- length(eta))) stop("eta[] must have length at least 1")#FIXME
  if(2 > (m <- as.integer(m)))stop(sQuote("m")," must be at least 2")#FIXME
  mhalfm <- (m-1) %/% 2:2

  ## partial derivatives of log f (at each Fourier frequency)
  lf <- matrix(ncol = M, nrow = mhalfm)

  f0 <- specARIMA(eta, p=p, q=q, m=m)$spec

  for(j in 1:M) {
    etaj <- eta
    etaj[j] <- etaj[j] + delta
    fj <- specARIMA(etaj, p=p, q=q, m=m)$spec
    lf[,j] <- log(fj/f0)/delta
  }

  ## Calculate D

  m* chol2inv(qr.R(qr(lf)))# == m* solve(crossprod(lf))

}## CetaARIMA()

specFGN <- function(eta, m, nsum = 200)
{
  ## Purpose: Calculation of the spectral density f of
  ## normalized fractional Gaussian noise with self-similarity parameter H=H
  ## at the Fourier frequencies 2*pi*j/m (j=1,...,(m-1)/2).
  ##
  ## Remarks:
  ## -------
  ## 1. cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
  ## 2. f=theta1*spec and integral[log(spec)]=0.
  ## -------------------------------------------------------------------------
  ## INPUT: m   = "sample size", 2 * #{Fourier frequencies}
  ##        eta = (H, *);  H = self-similarity parameter
  ##
  ## OUTPUT: list(spec=spec,theta1=theta1)
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

  ##---------parameters for the calculation of f--------

  if(2 > (m <- as.integer(m)))
      stop(sQuote("m")," must be at least 2")

  if(length(eta) > 1) warning("eta[2..] are not used")
  H <- eta[1]
  hh <- -2*H-1

  mhalfm <- (m-1) %/% 2:2
  ##-- Fourier frequencies lambda_j = 2*pi*(j-1)/m (j=1,2,...,(m-1)/2)
  lambd <- 2*pi/m * (1:mhalfm)

  ##-----   calculation of f at Fourier frequencies   -------
  j <- 2*pi*(1:nsum)
  ## Using outer() instead of the loop is  *slower* (!)
  ##spec <- colSums(abs(outer(j, lambd, "+"))^hh + abs(outer(j, lambd, "-"))^hh)
  spec <- double(mhalfm)
  for(i in 1:mhalfm)
      spec[i] <- sum(abs(j + lambd[i])^hh + abs(j - lambd[i])^hh)

  spec <- sin(pi*H)/pi * gamma(-hh) * (1-cos(lambd)) * (lambd^hh + spec)

  ##----- adjust spectrum such that int{log(spec)} = 0 -------
  theta1 <- exp(2/m * sum(log(spec)))
  ## and specS = spec*() = spec / theta1 is adjusted :
  r <- list(freq = lambd, spec = spec/theta1, theta1 = theta1,
            H = H, method = paste("fGN( H=",format(H),")"))
  class(r) <- "spec"
  r
}## specFGN()

specARIMA <- function(eta, p,q, m)
{
  ## Purpose: Calculation of the spectral density of fractional ARMA
  ##	with standard normal innovations and self-similarity parameter H=H
  ##    at the Fourier frequencies 2*pi*j/m (j=1,..., (m-1)/2).
  ## cov(X(t),X(t+k)) = (sigma/(2*pi))*integral(exp(iuk)g(u)du).
  ##
  ## Remarks:
  ## -------
  ## 1. cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
  ## 2. f=theta1*spec and integral[log(spec)]=0.
  ## -------------------------------------------------------------------------
  ## INPUT: m = sample size
  ##        H = theta[1] = self-similarity parameter
  ##        phi = theta[    2:(p+1)]   = AR(p)-parameters
  ##        psi = theta[(p+2):(p+q+1)] = MA(q)-parameters
  ##
  ## OUTPUT: list(spec=spec,theta1=theta1)
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

  ##---------parameters for the calculation of f--------

  if(0 > (p <- as.integer(p)))
	stop("`p'  must be non-negative integer")
  if(0 > (q <- as.integer(q)))
        stop("`q'  must be non-negative integer")
  if(1+p+q != (M <- length(eta)))
        stop("eta[] must have length 1+p+q")
  if(2 > (m <- as.integer(m)))
        stop(sQuote("m")," must be at least 2")
  mhalfm <- (m-1) %/% 2:2

  H <- eta[1]

  ##------   Fourier frequencies   ----------------------
  ##------   x = 2*pi*(j-1)/m (j=1,2,...,(n-1)/2)   -----
  x <- 2*pi/m * (1:mhalfm)

  ##-----   calculation of f at Fourier frequencies   -------


  if(p > 0) {
    phi <- eta[2:(p+1)]
    px <- outer(x, 1:p)
    Rar <- cos(px) %*% phi
    Iar <- sin(px) %*% phi

    far <- (1-Rar)^2 + Iar^2
  } else {
      phi <- numeric(0)
      far <- 1
  }

  if(q > 0) {
    psi <- eta[(p+2):(p+q+1)]
    px <- outer(x, 1:q)
    Rma <- cos(px) %*% psi
    Ima <- sin(px) %*% psi

    fma <- (1+Rma)^2 + Ima^2
  } else {
      psi <- numeric(0)
      fma <- 1
  }

  spec <- fma/far * sqrt(2 - 2*cos(x))^(1-2*H)
  ##   			 2 - 2*cos(x) == (1-cos(x))^2 + sin(x)^2

  r <- list(freq = x, spec = spec, theta1 = 1/(2*pi),
            pq = c(p,q), eta = c(H = H, phi = phi, psi = psi),
            method = paste("ARIMA(",p,", ", H - 1/2,", ", q,")",sep = ""))
  class(r) <- "spec"
  r
}## specARIMA()


per <- function(z) {
  ## Purpose:  Computation of the periodogram via FFT
  ## -------------------------------------------------------------------------
  ## Arguments: z : time-series
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

### MM: per[-1] is the same as
### --    1/(2*pi) * spec.pgram(z, fast=FALSE, taper=0, detrend=FALSE)$spec
###    "[-1]" : spec.pgram() doesn't try to use the *wrong* value lambda=0

### Note the "fast = FALSE" needed to emulate this !
  n <- length(z)
  (Mod(fft(z))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
}



Qeta <- function(eta, model = c("fGn","fARIMA"), n, yper,
                 pq.ARIMA, give.B.only = FALSE)
{
    ## Purpose: Calculation of A, B and Tn = A/B^2
    ## where A = 2pi/n sum 2*[I(lambda_j)/f(lambda_j)],
    ##       B = 2pi/n sum 2*[I(lambda_j)/f(lambda_j)]^2  and
    ## the sum is taken over all Fourier frequencies
    ## lambda_j = 2pi*j/n (j=1,...,(n-1)/2.
    ## f is the spectral density of fractional Gaussian
    ## noise or fractional ARIMA(p,d,q) with self-similarity parameter H=H.
    ## cov(X(t),X(t+k))=integral(exp(iuk)f(u)du)
    ##
    ## NOTE: yper[1] must be the periodogram I(lambda_1) at
    ## ----  the frequency 2pi/n (i.e. not the frequency zero !).
    ## -------------------------------------------------------------------------
    ## INPUT: H
    ##        (n,nhalfm = trunc[(n-1)/2] and the
    ##         nhalfm-dimensional  GLOBAL vector `yper' must be defined.)
    ##   (n,yper): now arguments {orig. w/ defaults 'n=n, yper=yper'}

    ## OUTPUT: list(n=n,H=H,A=A,B=B,Tn=Tn,z=z,pval=pval,
    ##              theta1=theta1,spec=spec)
    ##
    ##         Tn is the goodness of fit test statistic
    ##         Tn=A/B^2 defined in Beran (1992),
    ##         z is the standardized test statistic,
    ##         pval the corresponding p-value P(w>z).
    ##         theta1 is the scale parameter such that
    ##         f=theta1*spec and integral(log[spec])=0.
    ## -------------------------------------------------------------------------
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

    model <- match.arg(model)

    stopifnot(is.numeric(yper),
              length(yper) == (n - 1) %/% 2)

    ##         spectrum at Fourier frequencies

    if (model == "fARIMA") {
        p <- pq.ARIMA[1]
        q <- pq.ARIMA[2]
        stopifnot(is.numeric(p), is.numeric(q),
                  p == round(p), q == round(q), p >= 0, q >= 0)
    }
    sp <- switch(model,
		 "fGn"	  = specFGN  (eta, n ),
		 "fARIMA" = specARIMA(eta,p,q,n))
    ## not used (and *different* from 'theta1' below) : th1 <- sp $theta1
    spec <- sp$spec
    yf <- yper/spec
    B <- 2*(2*pi/n) * sum(yf)

    if(give.B.only)
        return(B)
    ## else

    A <- 2*(2*pi/n)*sum(yf*yf)
    Tn <- A/(B^2)
    z <- sqrt(n)*(pi*Tn-1)/sqrt(2)
    theta1 <- B/(2*pi)

    list(n = n, H = eta[1], eta = eta,
	 A = A, B = B, Tn = Tn, z = z,
	 pval = pnorm(z, lower.tail=FALSE), # = 1 - Phi(.)
	 theta1 = theta1, spec = spec)
}## Qeta()



###MMMMMMMMMMMMMMMMMMMMMMM
###
### MAIN PROGRAM  (fARIMA)
###
###MMMMMMMMMMMMMMMMMMMMMMM

## Martin's new main function  *INSTEAD* of any main program :
## NOTA BENE: no "loop" over different data segments
WhittleEst <- function(x,
		       ## periodogram of data:
		       periodogr.x = per(if(scale) x/sd(x) else x)[2:((n+1) %/% 2)],
		       n = length(x), scale = FALSE,
		       model = c("fGn","fARIMA"), p, q,
		       start = list(H= 0.5, AR= numeric(), MA=numeric()),
		       verbose = getOption("verbose"))
{
    ## Author: Martin Maechler, 2004-2005; 2011
    cl <- match.call()
    stopifnot(is.list(start), sapply(start, is.numeric),
	      length(start$H) == 1,
	      n >= 4)

    ## If 'fGn', don't need p nor q  -- maybe use default p = 0, q = 0
    model <- match.arg(model)
    if(model == "fGn") {
	stopifnot(length(start) == 1 || # namely only "$ H"
		  sapply(start[-1], length) == 0)## other components have 0 length
	if(!missing(p) && p != 0)
	    stop("'p' != 0 does not make sense in \"fGn\" model")
	if(!missing(q) && q != 0)
	    stop("'q' != 0 does not make sense in \"fGn\" model")
	p <- q <- 0
    }
    else { ## fARIMA
	if(missing(p)) p <- length(start$AR)
	else {
	    stopifnot(length(start$AR) == p)
	    if(0 > (p <- as.integer(p))) stop("must have integer p >= 0")
	}
	if(missing(q)) q <- length(start$MA)
	else {
	    stopifnot(length(start$MA) == q)
	    if(0 > (q <- as.integer(q))) stop("must have integer q >= 0")
	}
    }
    pq. <- c(p,q)

    eta <- unlist(start) # c(H, ar[1:p], ma[1:q]) # where p=0, q=0 is possible
    H0 <- eta[1]
    M <- length(eta)

    ## hmm: should give a warning when doing this:
    H <- max(0.2,min(H0,0.9))		# avoid extreme initial values
    eta[1] <- H

    ##---- find estimate -----------------------------------------------------------

    Qmin <- function(eta) {
	## Purpose: function to be minimized for MLE
	## in R, we nicely have access to all the variables in "main" function
	## ---------------------------------------------------------------------
	## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

	r <- Qeta(eta, model = model, n = n, yper = periodogr.x, pq.ARIMA = pq.,
		  give.B.only = TRUE)
	if(verbose) cat("Qmin(eta =",eta,") -> B =", r,"\n")
	drop(r)
    }

    etatry <- eta
    ##S+ s <- 2*(1 - H)
    ##S+ result <- nlmin(Qmin, etatry, xc.tol = 0.0000001, init.step = s)
    ##S+	  ----- FIXME: use optim() instead
    ##S+ eta <- result$x		   # = arg min _{eta} Qmin(eta, ..)
    if(M == 1) { # one-dimensional -- only 'H'
	res <- optimize(Qmin, lower = 0.1, upper = 0.99)
	eta <- c(H = res$minimum)

    } else { ## M > 1
	res <- optim(par = etatry, fn = Qmin)
	eta <- res$par
	## names(eta) <- c("H",
	##		paste("ar", seq_len(p), sep=""),
	##		paste("ma", seq_len(q), sep=""))
	if(res$convergence) {
	    warning("optim(<negative Log.Lik.>) did not converge: 'convergence'=",
		    res$convergence,"\n message:", res$message)
	}
    }


    ## calculate goodness of fit statistic
    Qresult <- Qeta(eta, model = model, n = n, yper = periodogr.x, pq.ARIMA = pq.)
    theta1 <- Qresult$theta1

    ## output
    if(verbose) {
	cat("	H    =", eta[1], fill = TRUE)
	cat("theta1  =", theta1, fill = TRUE)
	if(length(eta) > 1)
	    cat("eta[2:.]=", eta[-1], fill = TRUE)
    }
    Vcov <- switch(model,
		 "fGn"	 = CetaFGN  (eta),
		 "fARIMA" = CetaARIMA(eta,p,q)
		 ) / n
    ssd <- sqrt(diag(Vcov))
    Cf <- cbind(Estimate = eta, "Std. Error" = ssd,
		"z value" = eta/ssd,
		"Pr(>|z|)" = 2 * pnorm(-abs(eta/ssd)))
    dimnames(Vcov) <- list(names(eta), names(eta))
    ##	 Hlow <- eta[1]-1.96*sqrt(Vcov[1,1])
    ##	 Hup  <- eta[1]+1.96*sqrt(Vcov[1,1])
    ##	 if(verbose) cat("95%-C.I. for H: [",Hlow,",",Hup,"]\n")

    ## etalow <- eta - 1.96*ssd
    ## etaup <-	 eta + 1.96*ssd

    spec <- theta1 * Qresult$spec

    structure(list(call = cl, model = model, n = n, p=p, q=q,
		   coefficients = Cf, theta1 = theta1, vcov = Vcov,
		   periodogr.x = periodogr.x, spec = spec),
	      class = "WhittleEst")

}# WhittleEst()

if(getRversion() < "2.13") nobs <- function (object, ...) UseMethod("nobs")
nobs.WhittleEst <- nobs.FEXP <- function (object, ...) object$n

vcov.WhittleEst <- vcov.FEXP <- function (object, ...) object$vcov

## Define coef() methods --> confint.default() e.g., work correctly
coef1 <- function(object, ...) { ## care to keep {row}names :
    cf <- object$coefficients[, "Estimate", drop=FALSE]
    structure(as.vector(cf), names=dimnames(cf)[[1]])
}
coef.WhittleEst <- coef.FEXP <- coef1

linesSpec <- function (x, type = "l", col = 4, lwd = 2, ...) {
    ffr <- .ffreq(x$n)
    lines(ffr, x$spec, type=type, col=col, lwd=lwd, ...)
}

lines.WhittleEst <- lines.FEXP <- linesSpec


.WhittleModel <- function(mod)
    c("fGn" = "fractional Gaussian noise",
      "fARIMA" = "fractional ARIMA")[mod]

print.WhittleEst <-
    function (x, digits = getOption("digits"),
	      ## signif.stars = getOption("show.signif.stars"),
	      ...)
{
    stopifnot(is.character(mod <- x$model))
    longMod <- .WhittleModel(mod)
    cf <- x$coefficients
    stopifnot(is.numeric(H.hat <- cf["H","Estimate"]))
    cat(sprintf("'WhittleEst' Whittle estimator for  %s ('%s');	 call:\n",
		longMod, mod),
	paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n",
	if(mod == "fARIMA") sprintf("ARMA order (p=%d, q=%d); ", x$p, x$q),
	"\t  time series of length  n = ", x$n,
	".\n\n",
	sprintf("H = %s\ncoefficients 'eta' =\n", ## FIXME "theta = (theta1,eta) ???
		formatC(H.hat, digits= digits)), sep="")
    printCoefmat(cf, digits = digits,
		 signif.stars = FALSE, ## signif.stars = signif.stars,
		 na.print = "NA", ...)
    s.H <- cf["H", "Std. Error"]
    Hround <- function(x) round(x, max(2, digits - 4))
    seForm <- function(s) sprintf(" (%g)", Hround(s))
    cat(" <==> d := H - 1/2 = ", Hround(H.hat-0.5), seForm(s.H), "\n\n", sep="")
    ## FIXME: "theta" here, too
    str(x[length(x)-(2:0)], no.list = TRUE)
    invisible(x)
}


plot.WhittleEst <-
    function (x, log = "xy", type = "l",
              col.spec = 4, lwd.spec = 2,
              xlab = NULL, ylab = expression(hat(f)(nu)),
              main = paste(deparse(x$call)[1]), sub = NULL, ...)
{
    ### FIXME: improve 'main' -- rather indicate *estimated* H ! --
    n2 <- length(x$periodogr.x)
    n <- x$n
    ffr <- 2 * pi/n * seq_len(n2)
    if(is.null(xlab))
        xlab <- bquote(list(nu == 2*pi*j / n,
                            {} ~~ j %in% paste(1,ldots,.(n2)), n == .(n)))
    if(is.null(sub))
	sub <-
	    sprintf("Data periodogram and fitted (Whittle Estimator '%s') spectrum%s",
		    paste(x$model, if(x$model == "fARIMA")
			  sprintf("(p=%d, q=%d)", x$p, x$q)),
		    if(log == "xy") " [log - log]")

    plot(ffr, x$periodogr.x, log = log, type = type,
         xlab = xlab, ylab = ylab, main = main, sub = sub, ...)
    lines(ffr, x$spec, col = col.spec, lwd = lwd.spec)
}

linesSpec <- function (x, type = "l", col = 4, lwd = 2, ...) {
    ffr <- .ffreq(x$n)
    lines(ffr, x$spec, type=type, col=col, lwd=lwd, ...)
}

lines.WhittleEst <- lines.FEXP <- linesSpec


### TODO: summary method

