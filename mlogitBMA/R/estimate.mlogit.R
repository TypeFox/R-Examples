mnl.loglikelihood <- function(x, y, coefficients) {
	U <- mnl.utilities(x, coefficients)
	P <- mnl.probabilities(U)
	sum(y*log(P))
}
	
mnl.utilities <- function(x, coefficients, npars=dim(x)[2], nobs=dim(x)[1], nalts=dim(x)[3]) {
	pars3d <- aperm(array(coefficients, c(npars, nobs, nalts)), c(2,1,3))
	U <- x*pars3d
    apply(U, c(1,3), sum)
}
	
mnl.probabilities <- function(U) {
	eU <- exp(U)
    sum.eU <- apply(eU, 1, sum)
    eU/sum.eU
}

do.estimate.mlogit <- function(x, y, specification=NULL, method='BHHH', ...) {
	# x is in 3d format
	# y is either a vector of indices of the chosen alternatives, or a 0-1 matrix (nobs x nalts)
	P <- 0
	
	mnl.logLik <- function(param) {
		U <- mnl.utilities(x, param, npars, nobs, nalts)
		P <<- mnl.probabilities(U)
		sum(y*log(P))
	}

	mnl.logLik.sum <- function(param) {
		sum(mnl.logLik(param))
	}
		
	mnl.grad <- function(param) {
        d <- y - P
        apply(aperm(array(d, c(nobs,nalts,npars)), c(1,3,2))*x, c(1,2), sum)
	}
	    
	mnl.hess <- function(param) {
		g <- mnl.grad(param)
		-t(g)%*%g
	}
	
	npars <- dim(x)[2]
	nobs <- dim(x)[1]
	nalts <- dim(x)[3]

	params <- rep(0, npars)
	ll0 <- mnl.logLik(params)
	s <- Sys.time()
	maxlik <- maxLik(mnl.logLik, grad=mnl.grad, hess=mnl.hess, 
						start=params, 
						method=method, ...
						)
	time <- Sys.time() - s
	
	coefs <- maxlik$estimate
	hessian <- maxlik$hessian
	names(coefs) <- rownames(hessian) <- colnames(hessian) <- dimnames(x)[[2]]

    U <- mnl.utilities(x, coefs, npars, nobs, nalts)
    P <- mnl.probabilities(U)
    fitted.values <- P
    residuals <- y - fitted.values
    g <- maxlik$gradient
    conv <- (t(g)%*%solve(-hessian))%*%g
	result <- list(coefficients=coefs, logLik=maxlik$maximum, logLik0=ll0,
					aic=-2 * maxlik$maximum + 2 * npars,
					bic=-2 * maxlik$maximum + npars * log(nobs),
					iter=maxlik$iterations, hessian=hessian,
					gradient=g,
					fitted.values=fitted.values,
					residuals=residuals,
					specification = specification,
					convergence = conv,
					method=method,
					time=time,
					code = maxlik$code,
					message=maxlik$message,
					last.step=maxlik$last.spep
				)
	attr(result, 'class') <- 'mnl'
	invisible(result)
}

"estimate.mlogit" <- function (object, ...) UseMethod("estimate.mlogit")

"estimate.mlogit.formula" <- function(f, data, method='BHHH', choices=NULL, base.choice=1,
									varying = NULL, sep='.', ...) {
	mnlspec <- mnl.spec(f, data=data, choices=choices, base.choice=base.choice, 
						varying=varying, sep=sep)
	estimate.mlogit(mnlspec, data, method=method, ...)
}

"estimate.mlogit.mnl.spec" <- function(object, data, method='BHHH', ...) {
	mnl.data <- get.mnl.data(data, spec=object)
	invisible(do.estimate.mlogit(mnl.data$x, mnl.data$y, specification=object, method=method, ...))
}

"estimate.mlogit.bic.mlogit" <- function(object, ...) {
	estimate.mlogit(object$bma.specifications, ...)
}

"estimate.mlogit.list" <- function(object, data, verbose=TRUE, ...) {
	result <- list()
	for (ispec in 1:length(object)){
		if(verbose) cat('\n\nEstimating specification #', ispec, 
						'\n===============================\n')
		result[[ispec]] <- estimate.mlogit(object[[ispec]], data, ...)
		if(verbose) print(summary(result[[ispec]]))
	}
	invisible(result)
}

get.mnl.data <- function(x, y=NULL, spec, intercept.prefix='asc.') {
	# Returns a list with components: 
	#	x - 3d array (nobs x nvars x nchoices)
	#	y - 0-1 response matrix (nobs x nchoices) 
	data <- list()
	data$y <- y
	nchoices <- length(spec$choices)
	nx <- nrow(x)
	if(is.null(y)) data$y <- x[,spec$response]
	unique.y <- unique(data$y)
	if(!all(is.element(unique.y, 1:nchoices))) { # replace values of y by its index within choices
		ynew <- rep(NA, length(data$y))
		for(ichoice in 1:nchoices)
			ynew[data$y==spec$choices[ichoice]] <- ichoice
		keep.idx <- !is.na(ynew) 
		data$y <- ynew[keep.idx]
		x <- x[keep.idx,]
		nx <- nrow(x)
	}
	
	# create a response matrix from a vector
	ym <- matrix(0, ncol=nchoices, nrow=nx)
	index <- matrix(c(1:nx,data$y), ncol=2)
	ym[index] <- 1
	data$y <- ym
	
	other.choices <- spec$choices[-spec$base.choice]
	other.choices.idx <- (1:nchoices)[-spec$base.choice]
	
	lvars <- sum(spec$same.coefs) + sum(spec$variable.used[,!spec$same.coefs])
		
	col.names <- colnames(x)
	data$x <- array(NA, c(nx, lvars, nchoices))
	var.count <- 0
	x.col.names <- c()
	for (var.sh in colnames(spec$variable.used)) {
		if (spec$same.coefs[var.sh]) { #variables with the same coefficient
			var.count <- var.count + 1
			if (is.element(var.sh, col.names)){ # is the variable non-alternative specific name in the data?
				data$x[,var.count,] <- x[,var.sh]
			} else { # get alternative specific name from the data
				which.choice <- (1:nchoices)[is.element(spec$full.var.names[,var.sh], col.names)]
				alt.spec.vars <- spec$full.var.names[which.choice, var.sh]
				d <- x[,alt.spec.vars]
				if(is.element(spec$base.choice, which.choice)) {
					# generate differences
					d <- d - x[,spec$full.var.names[spec$base.choice, var.sh]]
				}
				data$x[,var.count,which.choice] <- as.matrix(d)
			}
			x.col.names <- c(x.col.names, var.sh)
		} else { # different coefficients for different alternatives
			var.full <- spec$full.var.names[,var.sh]
			var.base <- spec$full.var.names[spec$base.choice, var.sh]
			get.difference <- FALSE
			if (is.element(var.base, col.names)) get.difference <- TRUE
			for (ichoice in other.choices.idx) {
				if(spec$variable.used[ichoice,var.sh]) {
					var.count <- var.count + 1
					if (is.element(var.full[ichoice], col.names)) {#is the variable alt.-specific name in the data?
						# generate differences
						d <- x[,var.full[ichoice]]
						if (get.difference) d <- d - x[,var.base]
						data$x[,var.count,ichoice] <- d
					} else data$x[,var.count,ichoice] <- x[,var.sh]
					data$x[,var.count,-ichoice] <- 0
					x.col.names <- c(x.col.names, var.full[ichoice])
				}
			}
		}	
	}
	if(any(spec$intercepts)) {
		int.array <- array(0, c(nx, sum(spec$intercepts), nchoices))
		dimnames(int.array)[[3]] <- spec$choices
		j<-1
		for (ichoice in (1:nchoices)[spec$intercepts]) {
			int.array[,j,(1:nchoices)[ichoice]] <- 1
			j <- j+1
		}
		data$x <- abind(int.array, data$x, along=2)
		x.col.names <- c(paste(intercept.prefix, spec$choices[spec$intercepts], sep=''), x.col.names)
	}
	data$x[,,spec$base.choice] <- 0
	dimnames(data$x) <- list(1:nx, x.col.names, dimnames(data$x)[[3]])

	return(data)
}