#### Functions used for testing, demos  of   fitCopula()
####                    -------              -----------

##' Error catching version of fitCopula()
##' @param ... further arguments passed to fitCopula()
fitC <- function(copula, u, p, ...) {
    stopifnot(length(p) > 0, p == round(p))
    r <- tryCatch(fitCopula(copula, u, ...), error = function(e)e)
    if(is(r, "error")) {
	cat("Caught fitCopula() error: '", r$message,"'\n*-*-*\n", sep="")
	new("fitCopula", estimate = rep.int(NA_real_,p),
	    var.est = matrix(NA_real_, p,p), copula = copula)
    }
    else r
}

##' @title Fit one sample of size n
##' @param cop copula with *one* parameter to estimate {and specified if generating 'x'}
##' @param n sample size; only used if \code{x} is not specified
##' @param x numeric matrix (n x d) "as from the copula"
##' @param ... further arguments passed to fitC()
##' @return a matrix with 2 columns, "est" | "se"  and one row for each method
fit1 <- function(cop, n, x = rCopula(n, cop), ...) {
  u <- pobs(x)
  f.itau <- fitC(cop, u, p=1, method= "itau", ...)
  f.irho <- fitC(cop, u, p=1, method= "irho", ...)
  f.mpl	 <- fitC(cop, u, p=1, method= "mpl",  ...)
  f.ml	 <- fitC(cop, x, p=1, method= "ml",   ...)
  if(length(f.mpl@estimate) != 1)
    stop("'mpl' gave estimate of length != 1, but ", length(f.mpl@estimate))
  cbind(est = c(itau = f.itau@estimate, irho = f.irho@estimate,
		mpl  = f.mpl@estimate,	ml   = f.ml@estimate),
	se = c(itau = sqrt(f.itau@var.est), irho = sqrt(f.irho@var.est),
		mpl = sqrt(f.mpl@var.est),    ml = sqrt(f.ml@var.est)))
}


## N replicates of fit1()  --> mean of estimates and se, and sd of estimates
##' @param ... further arguments passed to fit1()
replFitCop <- function(cop, n, N, verbose=TRUE, na.rm=FALSE, ...) {
    if(verbose) { k <- 0
		  pb <- txtProgressBar(max = N, style = if(isatty(stdout())) 3 else 1)
		  on.exit(close(pb)) } # setup progress bar and close at end
    x.try <- rCopula(2, cop)
    if(!is.matrix(x.try) || any(is.na(x.try)))
        stop("copula 'cop' cannot be used to generate data; maybe missing parameter?")
    ## below replicate(N, { .. } ) now fails because of the "..." passing:
    sim <- sapply(seq_len(N), function(i, ...) {
        if(verbose) setTxtProgressBar(pb, (k <<- k + 1)) # update progress bar
        fit1(cop,n, , ...)
    }, ..., simplify="array")
    m <- rowMeans(sim, na.rm=na.rm, dims=2)
    cbind(bias = abs(m[,"est"] - cop@parameters), ## abs of bias
          dSE  = abs(m[,"se"] - apply(sim[,1,],1,sd))) ## abs of mean of se minus sd of estimates
}

##' Test of the methods in fitCopula for a bivariate one-parameter copula family
##'
##' @title Test of the methods in fitCopula for a bivariate one-parameter copula family
##' @param cop is the copula family
##' @param tau.set is the set of tau values for which the parameters of cop are set
##' @param n.set is the set of n values used in the simulation
##' @param N is the number of repetitions for computing the bias and dSE
##' @param verbose logical indicating if progress bar is desired
##' @param ... further arguments passed to replFitCop()
##' @return an array bias and dSE in different tau and n scenarios
##' @author Martin
tstFit1cop <- function(cop, tau.set, n.set, N, verbose=TRUE, ...) {
  theta <- structure(iTau(cop, tau.set),
		     names = paste("tau", format(tau.set), sep="="))
  names(n.set) <- paste("n",format(n.set),sep="=")
  ## setTheta() also works for tCopula(), tevCopula() with non-fixed df:
  cop.set <- lapply(theta, setTheta, x=cop)
  sapply(cop.set,
         function(cop) {
	      if(verbose) cat("tau = ", format(tau(cop)), "; copula par: ",
			      format(cop@parameters), "\n", sep="")
	      validObject(cop)
              sapply(n.set, function(n)
                {
                     if(verbose) cat("  n = ",n,"\n", sep="")
                     replFitCop(cop, n, N=N, verbose=verbose, ...)
                     ##--------
                }, simplify="array")
         }, simplify = "array")
}


##' Reshapes the results for processing by xyplot
##' @title Reshapes the results for processing by xyplot
##' @param res object returned by tstFit1cop
##' @param taunum logical indicating whether the taus are numeric or character
##' @return a matrix containing the reshaped results
##' @author Martin
reshape.tstFit <- function(res, taunum=FALSE) {
  names(dimnames(res)) <- c("method","stat","n","tau")
  d <- cbind(expand.grid(dimnames(res),KEEP.OUT.ATTRS=FALSE), x=as.vector(res))
  d[,"n"  ] <- as.numeric(vapply(strsplit(as.character(d[, "n" ]),split="="),`[`,"",2))
  if (taunum)
    d[,"tau"] <- as.numeric(vapply(strsplit(as.character(d[,"tau"]),split="="),`[`,"",2))
  d
}

## plot the results (after reshaping)
plots.tstFit <- function(d, log=TRUE, auto.key=TRUE, ...) {
  require(lattice)
  xyplot(x ~ n | stat*tau, groups=method, data=d, type="b",
         scales = list(x=list(log=log), y=list(log=log)),
         auto.key=TRUE, ...)
}
