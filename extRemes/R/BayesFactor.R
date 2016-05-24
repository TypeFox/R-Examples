BayesFactor <- function(m1, m2, burn.in = 499, FUN = "postmode", method = c("laplace", "harmonic"), verbose = FALSE) {

    if(class(m1) != "fevd") stop("BayesFactor: invalid m1 argument.")
    if(class(m2) != "fevd") stop("BayesFactor: invalid m2 argument.")

    method <- tolower(method)
    method <- match.arg(method)

    if(m1$method != "Bayesian") stop("BayesFactor: m1 did not use Bayesian estimation.")
    if(m2$method != "Bayesian") stop("BayesFactor: m2 did not use Bayesian estimation.")

    if(method == "laplace") {

	theta1 <- postmode(m1, burn.in = burn.in, verbose = verbose)
	theta2 <- postmode(m2, burn.in = burn.in, verbose = verbose)

	if(verbose) cat("Grabbing the data sets.\n")
    	y1 <- datagrabber(m1, cov.data = FALSE)
    	y2 <- datagrabber(m2, cov.data = FALSE)
	data1 <- datagrabber(m1, response = FALSE)
	data2 <- datagrabber(m2, response = FALSE)

	des1 <- setup.design(m1)
	des2 <- setup.design(m2)

	np1 <- dim(m1$results)[2] - 1
	np2 <- dim(m2$results)[2] - 1

	H1 <- try(optimHess(theta1, oevd, gr = grlevd, o = m1, des = des1, x = y1, data = data1, u = m1$threshold,
		npy = m1$npy, phi = m1$par.models$log.scale, blocks = m1$blocks))

	if(class(H1) == "try-error") {

	    warning("BayesFactor: unable to estimate the Hessian at the posterior mode for m1.  Using the covariance of the MCMC sample instead.")
	    if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

	        Sigma1 <- cov(m1$results[-(1:burn.in), -(np1 + 1)])

	    } else {

		Sigma1 <- cov(m1$results[, -(np1 + 1)])

	    } # end of if else 'burn.in' stmts.

	} else {

	    Sigma1 <- try(solve(H1))

            if(class(Sigma1) == "try-error") {

                warning("BayesFactor: unable to find the inverse Hessian at the posterior mode for m1.  Using the covariance of the MCMC sample instead.")
		if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

                    Sigma1 <- cov(m1$results[-(1:burn.in), -(np1 + 1)])

                } else {

                    Sigma1 <- cov(m1$results[, -(np1 + 1)])

                } # end of if else 'burn.in' stmts.
                    
            } # end of if finding the inverse failed for m1 stmts.

	} # end of if else 'try error' for m1 stmts

	H2 <- try(optimHess(theta2, oevd, gr = grlevd, o = m2, des = des2, x = y2, data = data2, u = m2$threshold,
                npy = m2$npy, phi = m2$par.models$log.scale, blocks = m2$blocks))

	if(class(H2) == "try-error") {

            warning("BayesFactor: unable to estimate the Hessian at the posterior mode for m2.  Using the covariance of the MCMC sample instead.")
            if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

                Sigma2 <- cov(m2$results[-(1:burn.in), -(np2 + 1)])

            } else {

                Sigma2 <- cov(m2$results[, -(np2 + 1)])

            } # end of if else 'burn.in' stmts.

        } else {

            Sigma2 <- try(solve(H2))

            if(class(Sigma2) == "try-error") {

                warning("BayesFactor: unable to find the inverse Hessian at the posterior mode for m2.  Using the covariance of the MCMC sample instead.")
                if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

                    Sigma2 <- cov(m2$results[-(1:burn.in), -(np2 + 1)])

                } else {

                    Sigma2 <- cov(m2$results[, -(np2 + 1)])

                } # end of if else 'burn.in' stmts.  

            } # end of if finding the inverse failed for m1 stmts.

        } # end of if else 'try error' for m1 stmts

	K1 <- try(chol(Sigma1))

	if(class(K1) == "try-error") det1 <- det(Sigma1) 
	else det1 <- exp(2 * sum(log(diag(K1)), na.rm = TRUE))

	K2 <- try(chol(Sigma2))

        if(class(K2) == "try-error") det2 <- det(Sigma2)
        else det2 <- exp(2 * sum(log(diag(K2)), na.rm = TRUE))

	ll1 <- oevd(p = theta1, o = m1, des = des1, x = y1, data = data1, u = m1$threshold,
		span = m1$span, npy = m1$npy, phi = m1$par.models$log.scale, blocks = m1$blocks)

	ll2 <- oevd(p = theta2, o = m2, des = des2, x = y2, data = data2, u = m2$threshold,
                span = m2$span, npy = m2$npy, phi = m2$par.models$log.scale, blocks = m2$blocks)

	p1 <- do.call(m1$priorFun, c(list(theta = theta1), m1$priorParams))
	p2 <- do.call(m2$priorFun, c(list(theta = theta2), m2$priorParams))

	pr1 <- (2 * pi)^(np1 / 2) * sqrt(det1) * ll1 * p1
	pr2 <- (2 * pi)^(np2 / 2) * sqrt(det2) * ll2 * p2

    } else {

	# if(verbose) cat("Grabbing the MCMC samples.\n")
	# par1 <- findAllMCMCpars(x = m1, burn.in = burn.in)
	# par2 <- findAllMCMCpars(x = m2, burn.in = burn.in)

	# if(verbose) cat("Calculating the likelihoods.\n")
	# l1 <- apply(par1, 1, lfun, x = m1, data = y1)
	# l2 <- apply(par2, 1, lfun, x = m2, data = y2)

	if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

	    l1 <- m1$chain.info[-(1:burn.in), "loglik"]
	    l2 <- m2$chain.info[-(1:burn.in), "loglik"]

	} else {

	    l1 <- m1$chain.info[, "loglik"]
            l2 <- m2$chain.info[, "loglik"]

	}

	pr1 <- 1 / mean(1 / l1, na.rm = TRUE)
	pr2 <- 1 / mean(1 / l2, na.rm = TRUE)

    } # end of if else 'method' stmts.

    dname <- c(deparse(substitute(m1)), " ", deparse(substitute(m2)))
    names(dname) <- c("model 1", "model 2")

    STATISTIC <- pr1 / pr2
    names(STATISTIC) <- "Bayes Factor (B12)"

    structure(list(statistic = STATISTIC, method = method, data.name = dname), class = "htest")

} # end of 'BayesFactor' function.

postmode <- function(x, burn.in = 499, verbose = FALSE, ...) {

    UseMethod("postmode")

} # end of 'postmode' function.

postmode.fevd <- function(x, burn.in = 499, verbose = FALSE, ...) {

    if(x$method != "Bayesian") stop("postmode.fevd: invalid estimation method.")

    if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

        if(verbose) cat("Removing first ", burn.in, " parameter vectors from the MCMC sample.\n")
        theta <- x$results[-(1:burn.in),]

    }

    pdim <- dim(theta)
    N <- pdim[2]
    m <- pdim[1]
    theta <- theta[,-N]
    thnames <- colnames(theta)

    if(is.element("log.scale", thnames)) {

	theta[ , thnames == "log.scale" ] <- exp(theta[, thnames == "log.scale" ])
	thnames[ thnames == "log.scale" ] <- "scale"
	colnames(theta) <- thnames

    }

    # if(verbose) cat("Obtaining all of the parameter values for each data point.\n")
    # par <- findAllMCMCpars(x = x, burn.in = burn.in)

    # if(verbose) cat("Grabbing the response data.\n")
    # y <- datagrabber(x, cov.data = FALSE)

    # if(verbose) cat("Finding the log-prior density for each iteration of the chain.\n")
    # priorD <- apply(theta, 1, function(p) do.call(x$priorFun, c(list(theta = p), x$priorParams)))

#     if(is.fixedfevd(x)) {
# 
# 	x2 <- x
# 
# 	if(verbose) cat("Finding log-likelihood for each iteration of the chain.\n")
# 
# 	res <- apply(par, 1, function(p) return(levd(x = y, threshold = p[4], location = p[1],
# 		    scale = p[2], shape = p[3], type = x2$type, negative = FALSE,
# 		    span = x2$span, npy = x2$npy)))
# 
# 	res <- res + priorD
# 
#     } else {
# 
# 	# Function to calculate the prior and likelihood
#         # for each iteration of the MCMC sample.
# 
#         hfun <- function(id, p, x, data, ind, priorD, vb) {
# 
# 	    if(vb && (id %% 100 == 0)) cat(id, " ")
# 	    id2 <- ind == id
# 
# 	    llik <- levd(x = data, threshold = p[id2,"threshold"], location = p[id2, "location"],
# 		    scale = p[id2, "scale"], shape = p[id2, "shape"], type = x$type,
# 		    negative = FALSE, span = x$span, npy = x$npy)
# 
# 	    res <- priorD[id] + llik
# 
# 	    return(res)
# 
#         } # end of internal 'hfun' function.
# 
#	if(verbose) cat("Finding h = log(likelihood * prior) for each iteration of the chain.\n")
#	res <- apply(matrix(1:m, ncol = 1), 2, p = par, x = x, data = y, ind = rep(1:x$n, each = dim(theta)[1]), priorD = priorD, vb = verbose)
#	if(verbose) cat("\nFinding the paramters equal to the maximum value of h.\n")
# 
#     } # end of if else model is fixed or varies with covariates stmt.

    if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) res <- x$chain.info[-(1:burn.in),"loglik"] + x$chain.info[-(1:burn.in),"prior"]
    else res <- x$chain.info[,"loglik"] + x$chain.info[,"prior"]

    ind <- res == max(res, na.rm = TRUE)
    if(sum(ind) > 1) out <- colMeans(theta[ind, ], na.rm = TRUE)
    else if(sum(ind) == 1) out <- theta[ind, ]
    else {

	warning("postmode.fevd: sorry, something went terribly wrong.  Returning NULL.")
	out <- NULL
    }

    return(out)

} # end of 'postmode.fevd' function.

findAllMCMCpars <- function(x, burn.in = 499, qcov = NULL, ...) {

    np <- dim(x$results)[2] - 1
    if(burn.in > 0) p <- x$results[-(1:burn.in), 1:np]
    else p <- x$results[, 1:np]
    pnames <- colnames(p)

    ni <- dim(p)[1]

    if(is.fixedfevd(x)) {

	if(!is.element("location", pnames)) loc <- rep(0, ni)
	else loc <- p[, "location"]

	if(is.element("log.scale", pnames)) scale <- exp(p[,"log.scale"])
	else scale <- p[,"scale"]

	if(!is.element("shape", pnames)) shape <- rep(0, ni)
	else shape <- p[, "shape"]

	if(is.null(x$threshold)) u <- rep(0, ni)
	else u <- x$threshold

	if(length(u) == 1) u <- rep(u, ni)

    } else {

	if(!is.null(qcov)) {

    	    nq <- dim(qcov)[1]
    
    	    # Set up location matrix
            if(is.element("location", pnames)) {
    
                loc <- matrix(p[, "location"], ni, nq)
                nloc <- 1
    
            } else if(is.element("mu", substring(pnames, 1, 2))) {
    
                nloc <- sum(substring(pnames, 1, 2) == "mu")
                loc <- matrix(NA, ni, nq)
    
                for(i in 1:nq) loc[,i] <- rowSums(p[, 1:nloc, drop = FALSE] * matrix(qcov[i, 1:nloc], ni, nloc, byrow = TRUE))
    
    
            } else {
    
                loc <- matrix(0, ni, nq)
                nloc <- 0
    
            } # end of setting up location matrix stmts.
    
    	    # Set up scale matrix.
            if(is.element("log.scale", pnames)) {
    
                nsc <- 1
                scale <- matrix(exp(p[, "log.scale"]), ni, nq)
    
            } else if(is.element("scale", pnames)) {
    
                nsc <- 1
                scale <- matrix(p[, "scale"], ni, nq)
    
            } else if(is.element("sig", substring(pnames, 1, 3))) {
    
                nsc <- sum(substring(pnames, 1, 3) == "sig")
                scale <- matrix(NA, ni, nq)
    
                for(i in 1:nq) scale[,i] <- rowSums(p[, (nloc + 1):(nloc + nsc), drop = FALSE] * matrix(qcov[i, (nloc + 1):(nloc + nsc)], ni, nsc, byrow = TRUE))
    
            } else if(is.element("phi", substring(pnames, 1, 3))) {
    
                nsc <- sum(substring(pnames, 1, 3) == "phi")
                scale <- matrix(NA, ni, nq)
    
                for(i in 1:nq) scale[,i] <- rowSums(p[, (nloc + 1):(nloc + nsc), drop = FALSE] * matrix(qcov[i, (nloc + 1):(nloc + nsc)], ni, nsc, byrow = TRUE))
    
                scale <- exp(scale)
    
            }
    
    	    # Set up shape matrix
            if(is.element("shape", pnames)) {
    
                nsh <- 1
                shape <- matrix(p[, "shape"], ni, nq)
    
            } else if(is.element("xi", substring(pnames, 1, 2))) {
    
                nsh <- sum(substring(pnames, 1, 2) == "xi")
                shape <- matrix(NA, ni, nq)
    
                for(i in 1:nq) rowSums(p[, (nloc + nsc + 1):(nloc + nsc + nsh), drop = FALSE] * matrix(qcov[i, (nloc + nsc + 1):(nloc + nsc + nsh)],
    				ni, nsh, byrow = TRUE))
    
            } else {
    
                nsh <- 0
                shape <- matrix(0, ni, nq)
    
            }
    
            u <- matrix(qcov[, "threshold"], ni, nq, byrow = TRUE)
    
    	    out <- cbind(c(loc), c(scale), c(shape), c(u))

        } else {

    	    designs <- setup.design(x)
    
            X.loc <- designs$X.loc
            if(!is.null(X.loc)) nloc <- ncol(X.loc)
            else nloc <- 0
    
            X.sc <- designs$X.sc
            nsc <- ncol(X.sc)
    
            X.sh <- designs$X.sh
            if(!is.null(X.sh)) nsh <- ncol(X.sh)
            else nsh <- 0
    
            # loc <- scale <- shape <- res <- matrix(NA, x$n, ni)
    
            if(is.null(X.loc)) loc <- matrix(0, x$n, ni)
            if(is.null(X.sh)) shape <- matrix(0, x$n, ni)
    
            if(is.null(x$threshold)) u <- matrix(0, x$n, ni)
            else u <- matrix(x$threshold, x$n, ni)
    
    
            parfinder <- function(z, y) {
    
                return(rowSums(t(z * t(y)), na.rm = TRUE))
    
            } # end of internal 'parfinder' function.
    
            if(is.null(X.loc)) loc <- matrix(0, ni, x$n)
            else loc <- apply(p[, 1:nloc, drop = FALSE], 1, parfinder, y = X.loc)
    
            scale <- apply(p[, (nloc+1):(nloc+nsc), drop = FALSE], 1, parfinder, y = X.sc)
            if(is.element("log.scale", pnames) || is.element("phi", substring(pnames, 1, 3))) scale <- exp(scale)
    
            if(is.null(X.sh)) shape <- matrix(0, ni, x$n)
            else shape <- apply(p[, (nloc+nsc+1):np, drop = FALSE], 1, parfinder, y = X.sh)

        } # end of if else 'qcov' stmts.

    } # end of if else fixed model stmts.

    # out <- list(loc = loc, scale = scale, shape = shape, threshold = u)
    out <- cbind(c(loc), c(scale), c(shape), c(u))
    colnames(out) <- c("location", "scale", "shape", "threshold")
    return(out)

} # end of 'findAllMCMCpars' function.
