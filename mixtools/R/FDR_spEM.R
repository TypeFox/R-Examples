#################################################################
## FUNCTIONS FOR FALSE DISCOVERY RATE (FDR) ESTIMATION         ##
## USING SEMI-PARAMETRIC EM-LIKE ALGORITHMS                    ##
## D. CHAUVEAU                                                 ##
## mixtools 1.0 addition                                       ##
#################################################################


#################################################################
## EM-like algorithm for a nonparametric univariate mixture model
## - Component 1 known = N(0,1) 
##   = the pdf of a probit transform of pvalue under H0
## - component 2 symmetric shifted from a location parameter mu
## NB: stochastic=TRUE not implemented, parameter removed here

spEMsymlocN01 <- function(x, mu0=2, bw = bw.nrd0(x), h=bw, eps = 1e-8,
						maxiter=100, verbose = FALSE, plotf=FALSE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  n <- length(x)
#  if (length(mu0)>1) m <- length(mu0) else m <- mu0 
  m <- 2 # fixed number of components in this model
  z.hat <- matrix(0, nrow=n, ncol=m)
  fkernel <- matrix(0, nrow=n, ncol=m)
  tt0 <- proc.time()
  o <- order(x) # for plotting fhat's if requested
  kmeans <- kmeans(x, mu0)  # is this a good init for probit data?
  for(j in 1:m) {
    z.hat[kmeans$cluster==j, j] <- 1
  }
  iter <- 0
  finished <- FALSE
  lambda <- matrix(0,maxiter,m)
  mu <- rep(0,maxiter) # only component 2 mean used in this case
  while (!finished) {
  #while (max(abs(change)) > eps & iter < maxiter) {
    iter <- iter + 1
    t0 <- proc.time()

    ## M-Step
    lambda[iter,] <- colMeans(z.hat) 
    # mu[iter,] <- apply(sweep(z.hat, 1, x, "*"), 2, mean)/lambda[iter,]
    mu[iter] <- sum(x*z.hat[,2])/(n*lambda[iter,2]) #

    ## second component density estimation step evaluated at x_i-mu's  
      ans <- .C("KDEsymloc1comp", n=as.integer(n), m=as.integer(m),
                mu=as.double(mu[iter]), lbd2=as.double(lambda[iter,2]),
                x=as.double(x), bw=as.double(bw),
                z=as.double(z.hat), f = double(n)) 
                #add PACKAGE="mixtools" option when integrated

    # successive plots of fhat's (for debugging mostly)
    if (plotf) {
   		if (iter==1) plotfunc <- plot else plotfunc <- lines 
    	plotfunc(x[o],ans$f[o],type="l", col=iter) 
    	}
    
    # version lambda_j f_j(xi) specific for component one = N(0,1)
    # NB: this is the only place where the known component pdf is used
    lambda.f <- cbind(lambda[iter,1]*dnorm(x), lambda[iter,2]*ans$f)

    ## E-step (for next iteration)
    z.hat <- lambda.f/rowSums(lambda.f)    
    finished <- iter >= maxiter
    if (iter>1) { # This convergence criterion is too simplistic:
      change <- c(lambda[iter,] - lambda[iter-1,], mu[iter]-mu[iter-1])
      finished <- finished | (max(abs(change)) < eps)
    }
    if (verbose) {
      t1 <- proc.time()
      cat("iteration ", iter, "  lambda ", round(lambda[iter,], 4), 
          "  mu ", round(mu[iter], 4))
      cat(" time", (t1 - t0)[3], "\n")
    }
  }
  if (verbose) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter,], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }
    return(structure(list(data=x, posteriors=z.hat, lambda=lambda[1:iter,],
                          bandwidth=bw, lambdahat=lambda[iter,], 
                          mu = mu[1:iter], muhat = mu[iter], symmetric=TRUE),
                          class="spEMN01"))
}





	
#######################################################
# plot mixture pdf for the semiparametric mixture model
# with component 1 pdf passed in the knownpdf parameter
# uses a weighted kernel density estimate of the nonparametric component 2 
# a = object of class "spEMN01" as returned by spEMsymlocNorm
plot.spEMN01 <- function(x, bw=x$bandwidth, knownpdf=dnorm, add.plot=FALSE, ...) {	
	t <- seq(min(x$data), max(x$data), len=200)
	f1 <- x$lambdahat[1]*knownpdf(t)
	f2 <- x$lambdahat[2]*wkde(x$data-x$muhat, u=t-x$muhat, w=x$post[,2], bw=bw, sym=TRUE)
	f <- f1+f2
	if (!add.plot) plot(t,f1+f2, type="l", ...) else lines(t,f1+f2, ...)
	lines(t,f1, col=2); lines(t,f2, col=3)
	}



####################################################
# plot and compare FDR for 1 or 2 EM-like strategies
# for mixtools 1.0
# post NEEDS to be sorted by p	
plotFDR <- function(post1, post2=NULL, lg1="FDR 1", lg2=NULL, title=NULL, 
					compH0=1, alpha=0.1, complete.data =NULL, pctfdr=0.3) {
	n <- dim(post1)[1]
	cs1 <- cumsum(post1[,compH0]) # local FDR(p_i)'s
	fdr1 <- cs1/(1:n) # FDR(p_i)'s
	if (is.null(title)) title <- paste("FDR estimate(s), n=",n)
	if (!is.null(post2)) {
			cs2 <- cumsum(post2[,compH0]) # local FDR(p_i)'s
			fdr2 <- cs2/(1:n)
			if (is.null(lg2)) lg2 <- "FDR 2"
			}
	i1 <- sum(fdr1<pctfdr) # plot up to pctfdr % FDR
	if (i1 == 0) i1 <- n   # for very bad fit, fdr[1] > pctfdr
	# cat("index",i1)
	plot(fdr1[1:i1], type="l", main=title, col=1, ylim=c(0,fdr1[i1]),
			xlab="index", ylab="probability")
	if (!is.null(post2)) lines(fdr2[1:i1], col=2)
	abline(alpha, 0, lty=3)
	if (!is.null(complete.data)) { # true complete data available
		V <- cumsum(complete.data[,1]==1) # cumulative nb of items under H0 
		trueFDR <- V/(1:n)
		lines(trueFDR[1:i1], lty=2, col=3)	
		if (!is.null(post2)) 
			legend("topleft", c(lg1,lg2,"True FDR"),col=1:3, lty=c(1,1,2))
		if (is.null(post2)) 
			legend("topleft", c(lg1,"True FDR"),col=c(1,3), lty=c(1,2))
				} else {
			if (!is.null(post2)) 
				legend("topleft", c(lg1,lg2), col=1:2, lty=c(1,1))
			if (is.null(post2)) 
				legend("topleft", lg1, col=1, lty=1)
			}
	}



