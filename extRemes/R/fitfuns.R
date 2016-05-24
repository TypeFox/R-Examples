fevd <- function(x, data, threshold=NULL, threshold.fun=~1, location.fun=~1, scale.fun=~1, shape.fun=~1, use.phi=FALSE,
    type=c("GEV", "GP", "PP", "Gumbel", "Exponential"),
    method=c("MLE", "GMLE", "Bayesian", "Lmoments"), initial=NULL, span, units=NULL, time.units="days", period.basis="year",
    na.action=na.fail, optim.args=NULL, priorFun=NULL, priorParams=NULL, proposalFun=NULL, proposalParams=NULL, 
    iter=9999, weights=1, blocks=NULL, verbose=FALSE) {
  
    if(verbose) begin.tiid <- Sys.time()
    out <- list()

    inout <- list() # object to return initial parameter values with likelihood.

    out$call <- match.call()
    if(!missing(data)) {

        out$data.name <- c(deparse(substitute(x)), deparse(substitute(data)))
	# cov.pointer <- as.character(substitute(data))

    } else {

	out$data.name <- c(deparse(substitute(x)), "")
	# cov.pointer <- NULL

    }

    # data.pointer <- as.character(substitute(x))
    # if( length(data.pointer) > 1) out$x <- x	
    # else out$data.pointer <- data.pointer
    # if(length(cov.pointer) > 1) out$cov.data <- data
    # else out$cov.pointer <- cov.pointer

    type <- match.arg(type)
    method <- match.arg(method)

    out$weights <- weights

    if(!missing(data)) {

        if(is.element(out$data.name[1], colnames(data))) {

            out$in.data <- TRUE
            wc <- out$data.name[1] == colnames(data)
            x <- c(data[,wc])
            x.fun <- ifelse(out$data.name[1] == "substitute(x)", deparse(x), out$data.name[1])
            x.fun <- formula(paste(x.fun, "~ 1"))
            out$x.fun <- x.fun

        } else if(is.formula(x)) {

            out$in.data <- TRUE
            x.fun <- x
            out$x.fun <- x.fun
            x <- model.response(model.frame(x.fun, data=data))

        } else out$in.data <- FALSE

        if(length(x) != nrow(data)) stop("fevd: data must have same number of rows as the length of x.")
        if(!identical(weights, 1) && length(x) != length(weights)) # CJP
          stop("fevd: weights should be the same length as x.") 

        tmp <- cbind(x, data)
        tmp <- na.action(tmp)
        x <- tmp[,1]
        data <- tmp[,-1, drop=FALSE]

    } else {

        if(is.formula(x)) stop("fevd: Must provide data argument if you supply a formula to the x argument.")
        x <- na.action(x)
        out$in.data <- FALSE

    }

    if(!out$in.data) {
        data.pointer <- as.character(substitute(x))
        if(length(data.pointer) > 1) out$x <- x
        else out$data.pointer <- data.pointer
    }

    if(!is.null(blocks)){  # CJP

      # test that blocks is of correct type
      if(type == "PP"){

        if(!is.element('nBlocks', names(blocks))){

          if(is.element('data', names(blocks))) {

            blocks$nBlocks <- nrow(blocks$data)

          } else stop("fevd: When supplying blocks, must provide 'blocks$nBlocks' if 'blocks$data' is  not provided.")

        }

        if(!is.element('weights', names(blocks))) blocks$weights <- 1

        if(!is.element('proportionMissing', names(blocks))) blocks$proportionMissing <- 0

        if(!is.element('threshold', names(blocks)) && !is.null(threshold)) {

          if(length(threshold) == 1) {

            blocks$threshold <- threshold

          } else {

            stop("fevd: No blocks$threshold specified and threshold is not constant. User must supply the threshold for each block via blocks$threshold.")

          }

        }

      } else {

        warning("fevd: Blocks are used only for type 'PP'; ignoring blocks argument.")
        blocks <- NULL # code will now operate as if user did not supply 'blocks'

      }

    }

    out$x <- x
    if(!missing(data)) out$cov.data <- data
    
    if(method == "MLE" && !is.null(priorFun)) method <- "GMLE"
    else if(method == "GMLE" && is.null(priorFun)) {

	if(shape.fun != ~1) stop("fevd: must supply a prior function for GMLE method when shape parameter varies.")

	if(is.element(type, c("GEV","GP","PP"))) {

	    priorFun <- "shapePriorBeta"
	    if(is.null(priorParams)) priorParams <- list(q=6, p=9)

	} else {

	    warning("fevd: method GMLE selected, but no priorFun given (and no default for the desired df type).  Switching to ML estimation.")
	    method <- "MLE"

	}

    } # end of if else method is GMLE but no prior function input stmts.

    if(method=="GMLE") {

	out$priorFun <- priorFun
	out$priorParams <- priorParams

    }

    out$method <- method
    method <- tolower(method)

    out$type <- type
    type <- tolower(type)

    out$period.basis <- period.basis
    out$optim.args <- optim.args

    out$units <- units

    if(method=="bayesian" && missing(use.phi)) {
	use.phi <- TRUE
	if(verbose) cat("\n", "Setting use.phi argument to TRUE for greater stability in estimation (default with Bayesian method).  Use use.phi=FALSE if you prefer that.\n")
    }

    # if(is.element(type, c("gp","pp","exponential","beta","pareto"))) out$npy <- npy 
    out$par.models <- list(threshold=threshold.fun, location=location.fun, scale=scale.fun, log.scale=use.phi, shape=shape.fun,
			    term.names=list(threshold=all.vars(threshold.fun), location=all.vars(location.fun),
			    scale=all.vars(scale.fun), shape=all.vars(shape.fun)))
    pars <- list()

    if(is.element(type, c("gp","pp","exponential","beta","pareto"))) {
	const.thresh <- check.constant(threshold.fun) & check.constant(threshold)
	out$const.thresh <- const.thresh
    }
    if(is.element(type, c("gev","pp","gumbel","weibull","frechet"))) {
	const.loc <- check.constant(location.fun)
	out$const.loc <- const.loc
    }
    const.scale <- check.constant(scale.fun)
    out$const.scale <- const.scale
    const.shape <- check.constant(shape.fun)
    out$const.shape <- const.shape
   
    if(is.element(type, c("pp", "gp", "exponential", "beta", "pareto"))) {
        if(missing(span)) {
          if(is.null(blocks)) {
            tiden <- attributes(x)$times
            n <- length(x)
            if(is.null(tiden)) {
              tiden <- 1:n
              start <- 1
              end <- n
              span <- end - start
            } else {
              start <- tiden[1]
              end <- tiden[n]
              span <- as.numeric(difftime(as.POSIXlt(tiden)[n], as.POSIXlt(tiden)[1], units=time.units))
            }
          } else {
            span <- blocks$nBlocks
          }
        } # end of if 'missing span' stmts.
        if(time.units=="days") npy <- 365.25
        else if(time.units=="months") npy <- 12
        else if(time.units=="years") npy <- 1
	else if(time.units == "hours") npy <- 24 * 365.25
	else if(time.units == "minutes") npy <- 60 * 24 * 365.25
	else if(time.units == "seconds") npy <- 60 * 60 * 24 * 365.25
        else {

            tmp.units <- unlist(strsplit(time.units, split="/")) 
            if(length(tmp.units) != 2) stop("fevd: invalid time.units argument.")
            numper <- as.numeric(tmp.units[1])
            if(is.na(numper)) stop("fevd: invalid time.units argument.")
            pertiid <- tmp.units[2]
            if(!is.element(pertiid, c("day","month","year","hour","minute","second"))) stop("fevd: invalid time.units argument.")
            if(pertiid=="year") npy <- numper
            else if(pertiid=="month") npy <- numper*12
            else if(pertiid=="day") npy <- numper*365.25
	    else if(pertiid == "hour") npy <- numper * 24 * 365.25
	    else if(pertiid == "minute") npy <- numper * 60 * 24 * 365.25
	    else if(pertiid == "second") npy <- numper * 60 * 60 * 24 * 365.25

        } # end of if complicated time units stmts.

        if(!is.null(blocks)) span <- span*npy # CJP - so is on the same scale as span w/o blocks
        out$time.units <- time.units
        out$span <- span/npy
        out$npy <- npy
        if(verbose)
          cat("\n", "Data span ", span/npy, "years", "\n") 
      } else npy <- NULL # end of if 'type is of the POT variety' stmts.
#     if(!missing(data)) {
# 	if(is.element(out$data.name[1], colnames(data))) {
# 	    out$in.data <- TRUE
# 	    wc <- out$data.name[1] == colnames(data)
# 	    x <- c(data[,wc])
# 	    x.fun <- ifelse(out$data.name[1] == "substitute(x)", deparse(x), out$data.name[1])
# 	    x.fun <- formula(paste(x.fun, "~ 1"))
# 	    out$x.fun <- x.fun
# 	} else if(is.formula(x)) {
# 	    out$in.data <- TRUE
# 	    x.fun <- x
# 	    out$x.fun <- x.fun
# 	    x <- model.response(model.frame(x.fun, data=data))
#  	} else out$in.data <- FALSE
# 	if(length(x) != nrow(data)) stop("fevd: data must have same number of rows as the length of x.")
# 	tmp <- cbind(x, data)
# 	tmp <- na.action(tmp)
# 	x <- tmp[,1]
# 	data <- tmp[,-1]
#     } else {
# 	if(is.formula(x)) stop("fevd: Must provide data argument if you supply a formula to the x argument.")
# 	x <- na.action(x)
# 	out$in.data <- FALSE
#     }
# 
#     if(!out$in.data) {
# 	data.pointer <- as.character(substitute(x))
# 	if(length(data.pointer) > 1) out$x <- x
# 	else out$data.pointer <- data.pointer
#     }

    n <- length(x)
    out$n <- n
    out$na.action <- deparse(substitute(na.action))

    if(!is.null(initial)) {
	if(!is.list(initial)) stop("fevd: initial must be NULL or a list object.")
	find.init <- FALSE
	if(is.null(initial$location) && is.element(type, c("gev","pp","gumbel","weibull","frechet"))) find.init <- TRUE
	if(use.phi && is.null(initial$log.scale)) find.init <- TRUE
	if(!use.phi && is.null(initial$scale)) find.init <- TRUE
	if(!is.element(type, c("gumbel","exponential")) && is.null(initial$shape)) find.init <- TRUE
    } else {
	initial <- list()
	find.init <- TRUE
    }

    if(method != "lmoments") {
	if(verbose) cat("Setting up parameter model design matrices.\n")
	designs <- list()
	if(!missing(data)) {
	    if(is.element(type, c("gp", "pp", "exponential", "beta", "pareto"))) X.u <- setup.design(x=threshold.fun, data=data, n=n, dname="threshold.fun") 
	    if(is.element(type, c("gev", "pp", "gumbel", "weibull", "frechet"))) {
	 	X.loc <- setup.design(x=location.fun, data=data, n=n, const=const.loc, dname="location.fun")
		designs$X.loc <- X.loc
	    }
	    X.sc <- setup.design(x=scale.fun, data=data, n=n, const=const.scale, dname="scale.fun") 
	    designs$X.sc <- X.sc
	    if(!is.element(type,c("gumbel","exponential"))) {
		X.sh <- setup.design(x=shape.fun, data=data, n=n, const=const.shape, dname="shape.fun")
		designs$X.sh <- X.sh
	    }
	} else {
	    if(is.element(type, c("gp", "pp", "exponential", "beta", "pareto"))) X.u <- setup.design(x=threshold.fun, n=n, dname="threshold.fun")
            if(is.element(type, c("gev", "pp", "gumbel", "weibull", "frechet"))) {
	 	X.loc <- setup.design(x=location.fun, n=n, const=const.loc, dname="location.fun")
		designs$X.loc <- X.loc
	    }
            X.sc <- setup.design(x=scale.fun, n=n, const=const.scale, dname="scale.fun")
	    designs$X.sc <- X.sc
            if(!is.element(type,c("gumbel","exponential"))) {
		X.sh <- setup.design(x=shape.fun, n=n, const=const.shape, dname="shape.fun")
		designs$X.sh <- X.sh
	    }
	} # end of if else missing data stmts.

        if(!is.null(blocks)){  # CJP
          blocks$designs <- list()
          if(is.element('data', names(blocks))){
            blocks$X.u <- setup.design(x = threshold.fun, data=blocks$data, n=blocks$nBlocks, dname="threshold.fun")
            blocks$designs$X.loc <- setup.design(x=location.fun, data=blocks$data, n=blocks$nBlocks, const=const.loc, dname="location.fun")
       	    blocks$designs$X.sc <- setup.design(x=scale.fun, data=blocks$data, n=blocks$nBlocks, const=const.scale, dname="scale.fun") 
            blocks$designs$X.sh <- setup.design(x=shape.fun, data=blocks$data, n=blocks$nBlocks, const=const.shape, dname="shape.fun")
          } else {
	    blocks$X.u <- setup.design(x=threshold.fun, n=blocks$nBlocks, dname="threshold.fun")
            blocks$designs$X.loc <- setup.design(x=location.fun, n=blocks$nBlocks, const=const.loc, dname="location.fun")
            blocks$designs$X.sc <- setup.design(x=scale.fun, n=blocks$nBlocks, const=const.scale, dname="scale.fun")
            blocks$designs$X.sh <- setup.design(x=shape.fun, n=blocks$nBlocks, const=const.shape, dname="shape.fun")
          }
	}
        
	if(verbose) cat("Parameter model design matrices set up.\n")
    } # end of if method is not Lmoments to set up design matrices stmts.

    if(is.element(type, c("gp", "pp", "exponential", "beta", "pareto"))) {

	if(method != "lmoments") threshold <- rowSums(matrix(threshold, n, ncol(X.u), byrow=TRUE) * X.u)
        if(!is.null(blocks)) # CJP
          blocks$threshold <- rowSums(matrix(blocks$threshold, blocks$nBlocks, ncol(blocks$X.u), byrow=TRUE) * blocks$X.u)
        
	excess.id <- x > threshold
	if(all(threshold==threshold[1])) out$threshold <- threshold[1]
	else out$threshold <- threshold

	out$rate <- mean(excess.id)
    } # end of if '!missing(threshold)' stmts.

    out$blocks <- blocks # CJP

    if(method=="lmoments" || find.init) {
	if(method=="lmoments") {
	    if(verbose) cat("Beginning estimation procedure.\n")
	    is.constant <- unlist(lapply(list(u=threshold, loc=location.fun, scale=scale.fun, sh=shape.fun), check.constant))
	    if(!all(is.constant)) warning("fevd: For method Lmoments, this function does not handle covariates in the parameters.  Fitting w/o covariates.")
	    if(!is.element(type, c("gev","gp"))) stop("fevd: currently, Lmoments are only handled for estimation of GEV and GP distribution parameters.")
	}
	xtemp <- x
	class(xtemp) <- "lmoments"
	ipars1 <- try(initializer(xtemp, model=type, threshold=threshold, npy=npy, blocks=blocks))
	if(class(ipars1) != "try-error") { # Eric -- 8/6/13
	    if(ipars1["scale"] <= 0) ipars1["scale"] <- 1e-8
	    if(method=="lmoments") {
	        out$results <- ipars1
	        class(out) <- "fevd"
	        return(out)
	    }
	} else {

	    ipars1 <- NULL
	    if(method == "lmoments") stop("fevd: Sorry, could not find L-moments estimates.")

	}
    } # end of if 'Lmoments method or find them anyway as potential initial estimates' stmts.

    if((method != "lmoments") && find.init) {

	xtemp <- x
	class(xtemp) <- "moms"
	ipars2 <- try(initializer(xtemp, model=type, threshold=threshold, npy=npy, blocks=blocks))

	if(class(ipars2) != "try-error") {
	    if(ipars2["scale"] <= 0) ipars2["scale"] <- 1e-8
	} else ipars2 <- NULL

	if(!is.element(type, c("pp", "gp","exponential","beta","pareto","gumbel"))) {

	    if(!is.null(ipars1)) testLmoments <- levd(x, location=ipars1["location"], scale=ipars1["scale"], shape=ipars1["shape"], type=out$type, npy=npy)
	    else testLmoments <- Inf

	    if(!is.null(ipars2)) testMoments <- levd(x, location=ipars2["location"], scale=ipars2["scale"], shape=ipars2["shape"], type=out$type, npy=npy)
	    else testMoments <- Inf

	} else if(type=="pp") { # CJP

	if(!is.null(ipars1)) { # Eric -- 8/6/13

          if(!is.null(blocks)) {
            blocks$location=ipars1["location"]; blocks$scale=ipars1["scale"]; blocks$shape=ipars1["shape"]
          }
          testLmoments <- levd(x, threshold=threshold, location=ipars1["location"], scale=ipars1["scale"], shape=ipars1["shape"], type=out$type, npy=npy, blocks=blocks)

	} else testLmoments <- Inf
	# end of if 'ipars1' not NULL stmts.

	if(!is.null(ipars2)) { # Eric -- 8/6/13

          if(!is.null(blocks))  {
            blocks$location=ipars2["location"]; blocks$scale=ipars2["scale"]; blocks$shape=ipars2["shape"]
          }

          testMoments <- levd(x, threshold=threshold, location=ipars2["location"], scale=ipars2["scale"], shape=ipars2["shape"], type=out$type, npy=npy, blocks=blocks)

	} else testMoments <- Inf
        # end of if 'ipars2' not NULL stmts.

          if(!is.null(blocks)) # CJP
            blocks$location <- blocks$scale <- blocks$shape <- NULL

	} else if(!is.element(type, c("gumbel","exponential"))) {

	    if(!is.null(ipars1)) testLmoments <- levd(x, threshold=threshold, scale=ipars1["scale"], shape=ipars1["shape"], type=out$type, npy=npy)
	    else testLmoments <- Inf

            if(!is.null(ipars2)) testMoments <- levd(x, threshold=threshold, scale=ipars2["scale"], shape=ipars2["shape"], type=out$type, npy=npy)
	    else testMoments <- Inf

	} else if(type=="gumbel") {

	    if(!is.null(ipars1)) testLmoments <- levd(x, location=ipars1["location"], scale=ipars1["scale"], type=out$type, npy=npy)
	    else testLmoments <- Inf

            if(!is.null(ipars2)) testMoments <- levd(x, location=ipars2["location"], scale=ipars2["scale"], type=out$type, npy=npy)
	    else testMoments <- Inf

        } else if(type=="exponential") {

	    if(!is.null(ipars1)) testLmoments <- levd(x, threshold=threshold, scale=ipars1["scale"], shape=0, type=out$type, npy=npy)
	    else testLmoments <- Inf

	    if(!is.null(ipars2)) testMoments <- levd(x, threshold=threshold, scale=ipars2["scale"], shape=0, type=out$type, npy=npy)
	    else testMoments <- Inf

	}

	if(is.finite(testLmoments) || is.finite(testMoments)) { # Eric -- 8/6/13 
	    if(testLmoments < testMoments) {

	        if(is.null(initial$location) && !is.element(type, c("gp","exponential","beta","pareto"))) initial$location <- ipars1["location"]

	        if(is.null(initial$log.scale) && use.phi) initial$log.scale <- log(ipars1["scale"])
	        else if(is.null(initial$scale)) initial$scale <- ipars1["scale"]

	        if(!is.element(type, c("gumbel","exponential")) && is.null(initial$shape)) initial$shape <- ipars1["shape"]

	        if(verbose) cat("Using Lmoments estimates as initial estimates.  Initial value =", testLmoments, "\n")

	    } else {

	        if(is.null(initial$location) && !is.element(type, c("gp","exponential","beta","pareto"))) initial$location <- ipars2["location"]

                if(is.null(initial$log.scale) && use.phi) initial$log.scale <- log(ipars2["scale"])
	        else if(is.null(initial$scale)) initial$scale <- ipars2["scale"]

                if(!is.element(type, c("gumbel","exponential")) && is.null(initial$shape)) initial$shape <- ipars2["shape"]

	        if(verbose) cat("Initial estimates found where necessary (not from Lmoments).  Initial value =", testMoments, "\n")

	    }
	} else { # Eric -- 8/6/13; if both L-moments and Moments fails, set initial values to mu = 0, sig = 1, xi = 0.01

	    if(is.null(initial$location) && !is.element(type, c("gp","exponential","beta","pareto"))) initial$location <- 0 

	    if(is.null(initial$log.scale) && use.phi) initial$log.scale <- 0
            else if(is.null(initial$scale)) initial$scale <- 1

	    if(!is.element(type, c("gumbel","exponential")) && is.null(initial$shape)) initial$shape <- 0.01

	    warning("fevd: L-moments and Moment initial parameter estimates could not be calculated.  Using arbitrary starting values.")
	}

	inout <- list(Lmoments=list(pars=ipars1, likelihood=testLmoments), MOM=list(pars=ipars2, likelihood=testMoments))

    } # end of if 'method is not lmoments and need to find initial estimates' stmts.

    # Set up initial estimates.
    if(!is.null(initial$location)) {
	if(ncol(X.loc) != length(initial$location)) {
	    if((length(initial$location)==1) && ncol(X.loc) > 1) initial$location <- c(initial$location, rep(0, ncol(X.loc)-1))
	    else stop("fevd: initial parameter estimates must have length 1 or number of parameters present.  Incorrect number for location parameter.")
	}
	if(length(initial$location)==1) names(initial$location) <- "location"
	else names(initial$location) <- paste("mu", 0:(ncol(X.loc)-1), sep="")
    }

    if(use.phi && (ncol(X.sc) != length(initial$log.scale))) {
	if((length(initial$log.scale)==1) && ncol(X.sc) > 1) initial$log.scale <- c(initial$log.scale, rep(0, ncol(X.sc)-1))
	else stop("fevd: initial parameter estimates must have length 1 or number of parameters present.  Incorrect number for log(scale) parameter.")
    } else if(!use.phi && (ncol(X.sc) != length(initial$scale))) {
	if((length(initial$scale)==1) && ncol(X.sc) > 1) initial$scale <- c(initial$scale, rep(0, ncol(X.sc)-1))
        else stop("fevd: initial parameter estimates must have length 1 or number of parameters present.  Incorrect number for scale parameter.")
    }

    if(use.phi) {
        if(length(initial$log.scale)==1) names(initial$log.scale) <- "log.scale"
        else names(initial$log.scale) <- paste("phi", 0:(ncol(X.sc)-1), sep="")
    } else {
	if(length(initial$scale)==1) names(initial$scale) <- "scale"
	else names(initial$scale) <- paste("sigma", 0:(ncol(X.sc)-1), sep="")
    }
    

    if(!is.element(type, c("gumbel","exponential"))) {
        if(ncol(X.sh) != length(initial$shape)) {
	    if((length(initial$shape)==1) && ncol(X.sh) > 1) initial$shape <- c(initial$shape, rep(0, ncol(X.sh)-1))
            else stop("fevd: initial parameter estimates must have length 1 or number of parameters present.  Incorrect number for shape parameter.")
        }
        if(length(initial$shape)==1) names(initial$shape) <- "shape"
        else names(initial$shape) <- paste("xi", 0:(ncol(X.sh)-1), sep="")
    } # end of if type is not light tailed stmts.

    if(is.element(method, c("mle","gmle"))) {
	if(use.phi) init.pars <- c(initial$location, initial$log.scale, initial$shape)
	else init.pars <- c(initial$location, initial$scale, initial$shape)
	if(type=="exponential" && const.scale) {

	    res <- list()
	    excess.id <- x > threshold
	    mle <- mean(x[excess.id] - threshold[excess.id])
	    names(mle) <- "scale"
	    res$par <- mle
	    k <- sum(excess.id)
	    res$n <- k
	    res$value <- k*(log(mle) + 1)

	} else {

	    if(!is.null(a <- optim.args)) {
	        anam <- names(a)
	        if(!is.element("gr",anam)) {
		    if(method=="mle") opt.gr <- grlevd
		    else opt.gr <- NULL
	        } else opt.gr <- a$gr
	        if(is.null(a$method) && use.phi) opt.method <- ifelse(is.element(type,c("gev","gp","pp", "gumbel")), "BFGS", "L-BFGS-B")
	        else opt.method <- a$method
	        if(is.element(type, c("weibull","beta","frechet","pareto"))) opt.method <- "L-BFGS-B"

		if(is.element(opt.method, c("L-BFGS-B","Brent"))) {
	            if(is.null(a$lower)) {
	                if(!is.element(type,c("frechet","pareto"))) opt.lower <- -Inf
	                else opt.lower <- c(rep(-Inf, length(init.pars)-1), 0)
	            } else {
	                    if(is.element(type,c("frechet","pareto"))) opt.lower <- c(a$lower[1:(length(init.pars)-1)], 0)
	                    else opt.lower <- a$lower
	            }
	            if(is.null(a$upper)) {
	               if(!is.element(type,c("weibull","beta"))) opt.upper <- Inf
	               else opt.upper <- c(rep(Inf, length(init.pars)-1), 0)
	            } else {
	               if(is.element(type,c("weibull","beta"))) opt.upper <- c(a$upper[1:(length(init.pars)-1)], 0)
	               else opt.upper <- a$upper
	            }
	        } else {
		    opt.lower <- -Inf
		    opt.upper <- Inf
		} # end of if 'optimization method allows bounds' stmts.

	        if(is.null(a$control)) opt.control <- list()
	        else opt.control <- a$control
	        anam <- names(a$control)
	        if(!is.element("trace", anam) && verbose) opt.control$trace <- 6

	       if(is.null(a$hessian)) opt.hessian <- TRUE
	       else opt.hessian <- a$hessian
	    } else {
	        if(method=="mle") opt.gr <- grlevd
		else opt.gr <- NULL
	        if(is.element(type, c("gev","gp","pp","gumbel"))) opt.method <- "BFGS"
	        else opt.method <- "L-BFGS-B"
	        if(!is.element(type,c("frechet","pareto"))) opt.lower <- -Inf
	        else opt.lower <- c(rep(-Inf, length(init.pars)-1), 0)
	        if(!is.element(type,c("weibull","beta"))) opt.upper <- Inf
	        else opt.upper <- c(rep(Inf, length(init.pars)-1), 0)
	        if(verbose) opt.control <- list(trace=6)
	        else opt.control <- list()
	        opt.hessian <- TRUE
	    } # end of if else 'optim.args' is passed stmts.

	    parnames <- names(init.pars)
	    out$parnames <- parnames
	    if(verbose && (method != "lmoments")) {
                cat("Initial estimates are:\n")
                print(init.pars)
                cat("Beginning estimation procedure.\n")
            }

	    if(type=="pp" && find.init) {

	 	# Yet a third attempt to find good starting values--this time using the MLE for the GP parameters.
	 	if(verbose) cat("\n", "First fitting a GP-Poisson model in order to try to get a good initial estimate as PP likelihoods can be very unstable.\n")
	 	look <- out
	 	look$type <- "GP"
		des2 <- designs
		des2$X.loc <- NULL
 	 	if(!missing(data)) resGP <- optim(init.pars[-(1:ncol(X.loc))], oevd, gr=opt.gr, o=look, des=des2, x=x, data=data, u=threshold, npy=npy, phi=use.phi, method=opt.method, lower=opt.lower, upper=opt.upper, control=opt.control, hessian=opt.hessian)
                 else resGP <- optim(init.pars[-(1:ncol(X.loc))], oevd, gr=opt.gr, o=look, des=des2, x=x, u=threshold, npy=npy, phi=use.phi, method=opt.method, lower=opt.lower, upper=opt.upper, control=opt.control, hessian=opt.hessian)
	 	tmpi <- resGP$par
                if(is.null(blocks)) { # CJP
                  lrate <- npy * mean(x > threshold)
                } else {
                  lrate <- sum(x > threshold)/(blocks$nBlocks * mean(blocks$weights))
                }
	 	xi3 <- tmpi[(ncol(X.sc)+1):length(tmpi)]
	 	if(!use.phi) sigma3 <- exp(tmpi[1:ncol(X.sc)] + xi3*log(lrate))
		else sigma3 <- tmpi[1:ncol(X.sc)] + xi3*log(lrate)
	 	lp <- lrate^(-xi3) - 1
	 	if(all(is.finite(lp))) mu3 <- mean(threshold) - (sigma3/xi3)*lp # CJP fixed bug - "mean(threshold)" not "threshold"
	 	else mu3 <- mean(x) # CJP: caution here as mean(x) will differ when only exceedances supplied... 
		nloc <- ncol(X.loc)
		if(length(mu3) < nloc) mu3 <- c(mu3, rep(0, nloc-length(mu3)))
		else mu3 <- mu3[1]

                if(!is.null(blocks)) { # CJP
                  blocks$location <- rowSums(matrix(mu3, blocks$nBlocks, nloc)*blocks$designs$X.loc)
                  blocks$scale=rowSums(matrix(sigma3, blocks$nBlocks, ncol(blocks$designs$X.sc))*blocks$designs$X.sc)
                  blocks$shape=rowSums(matrix(xi3, blocks$nBlocks, ncol(blocks$designs$X.sh))*blocks$designs$X.sh)
                }

		if(all(is.finite(c(mu3, sigma3, xi3)))) {

		    testGPmle <- try(levd(x=x, threshold=threshold, location=rowSums(matrix(mu3, n, nloc)*X.loc),
				    scale=rowSums(matrix(sigma3, n, ncol(X.sc))*X.sc),
				    shape=rowSums(matrix(xi3, n, ncol(X.sh))*X.sh), type="PP", npy=npy, blocks=blocks))

		    if(class(testGPmle) == "try-error") testGPmle <- Inf

		} else testGPmle <- Inf

                if(!is.null(blocks)) # CJP
                  blocks$location <- blocks$scale <- blocks$shape <- NULL

		if(is.finite(testLmoments) || is.finite(testMoments) || is.finite(testGPmle)) {
		    if((testGPmle < testLmoments) && (testGPmle < testMoments)) {
		        if(verbose) cat("\n", "Changing initial estimates to those based on GP MLEs.  They are: \n") 
		        if(use.phi) init.pars <- c(mu3, log(sigma3), xi3)
		        else init.pars <- c(mu3, sigma3, xi3)
		        names(init.pars) <- parnames
		        if(verbose) print(init.pars)
		    } else if(verbose) cat("\n", "Sticking with originally chosen initial estimates.\n")
		} # end of if any one of the three initial estimate attempts was successfull stmts.

		inout$PoissonGP <- list(pars=c(mu3, sigma3, xi3), likelihood=testGPmle)

	    } # end of if 'type is PP, then try to improve initial estimates' stmts.

	    if(method=="mle") {

	        if(!missing(data)) {

		    res <- optim(init.pars, oevd, gr=opt.gr, o=out, des=designs, x=x, data=data, u=threshold,
			span=span/npy, npy=npy, phi=use.phi, blocks=blocks, method=opt.method, lower=opt.lower,
			upper=opt.upper, control=opt.control, hessian=opt.hessian)  # CJP

	        } else {

		    res <- optim(init.pars, oevd, gr=opt.gr, o=out, des=designs, x=x, u=threshold,
			span=span/npy, npy=npy, phi=use.phi, blocks=blocks, method=opt.method, lower=opt.lower,
			upper=opt.upper, control=opt.control, hessian=opt.hessian)

		}

	    } else if(method == "gmle") {

		if(!missing(data)) {

		    res <- optim(init.pars, oevdgen, gr = opt.gr, o = out, des = designs, x = x, data = data,
			u = threshold, span = span / npy, npy = npy, phi = use.phi, blocks = blocks,
			priorFun = priorFun, priorParams = priorParams, method = opt.method, lower = opt.lower,
			upper = opt.upper, control = opt.control, hessian = opt.hessian)

                } else {

		    res <- optim(init.pars, oevdgen, gr = opt.gr, o = out, des = designs, x = x, u = threshold,
			span = span / npy, npy = npy, phi = use.phi, blocks = blocks, priorFun = priorFun,
			priorParams = priorParams, method = opt.method, lower = opt.lower, upper = opt.upper,
			control = opt.control, hessian = opt.hessian)

		}

	    } # end of if method is GMLE stmts.

	} # end of if else type is exponential.

    if(is.element("shape", names(res$par))) {
        if(is.element(type, c("frechet","pareto"))) {
	    res$par["shape"] <- abs(res$par["shape"])
	    if(res$par["shape"] == 0) {
		warning("fevd: shape parameter estimated to be zero.  Re-setting to be 1e16.")
		res$par["shape"] <- 1e16
	    }
	} else {
	    if(is.element(type, c("weibull","beta"))) res$par["shape"] <- -abs(res$par["shape"])
	    if(res$par["shape"] == 0) {
                warning("fevd: shape parameter estimated to be zero.  Re-setting to be -1e16.")
                res$par["shape"] <- -1e16
            }
	}
    }
    res$num.pars <- list(location=ncol(designs$X.loc), scale=ncol(designs$X.sc), shape=ncol(designs$X.sh))
    out$results <- res
    } else if(method=="bayesian") {
	
	if(is.element(type, c("gev","gumbel","weibull","frechet","pp"))) {
	    nloc <- ncol(X.loc)
	    loc.names <- names(initial$location)
	} else {
	    nloc <- 0
	    loc.names <- NULL
	}

	nsc <- ncol(X.sc)
	if(use.phi && is.null(initial$log.scale)) {
	    initial$log.scale <- log(initial$scale)
	    if(nsc == 1) names(initial$log.scale) <- "log.scale"
	    else names(initial$log.scale) <- paste("phi", 0:(nsc-1), sep="")
	}
	sc.names <- names(initial$log.scale)

	if(!is.element(type, c("gumbel","exponential"))) {
	    nsh <- ncol(X.sh)
	    sh.names <- names(initial$shape)
	} else {
	    nsh <- 0
	    sh.names <- NULL
	}

	np <- nloc + nsc + nsh

	find.priorParams <- FALSE
	if(is.null(priorFun) && is.null(priorParams)) find.priorParams <- TRUE
	else if(is.null(priorFun) && (is.null(priorParams$m) || is.null(priorParams$v))) find.priorParams <- TRUE
	else if(!is.null(priorFun)) {
	    if(priorFun == "fevdPriorDefault") {
		if(is.null(priorParams)) find.priorParams <- TRUE
		else if(is.null(priorParams$m) || is.null(priorParams$v)) find.priorParams <- TRUE
	    }
	}

	if(is.null(priorFun) || find.priorParams) {

	    if(is.null(priorFun)) priorFun <- "fevdPriorDefault"

	    if(find.priorParams) {

		xtemp <- x
		class(xtemp) <- "mle"
		if(verbose) cat("\n", "Finding MLE to obtain prior means and variances.\n")
		if(missing(data)) {

		    if(missing(span)) hold <- initializer(xtemp, u=threshold, use.phi=use.phi, type=out$type, time.units=time.units, period.basis=period.basis, blocks=blocks)
		    else hold <- initializer(xtemp, u=threshold, use.phi=use.phi, type=out$type, span=span, time.units=time.units, period.basis=period.basis, blocks=blocks)

		} else {

		    if(missing(span)) hold <- initializer(xtemp, data=data, u=threshold, u.fun=threshold.fun, loc.fun=location.fun, sc.fun=scale.fun, sh.fun=shape.fun,
					    use.phi=use.phi, type=out$type, time.units=time.units, period.basis=period.basis, blocks=blocks)
		    else hold <- initializer(xtemp, data=data, u=threshold, u.fun=threshold.fun, loc.fun=location.fun, sc.fun=scale.fun, sh.fun=shape.fun,
                                            use.phi=use.phi, type=out$type, span=span, time.units=time.units, period.basis=period.basis, blocks=blocks)

		}

		if(is.null(priorParams)) priorParams <- list(m=hold[1:np], v=rep(10,np))
		else if(is.null(priorParams$m)) priorParams$m <- hold[1:np]
		else if(is.null(priorParams$v)) priorParams$v <- rep(10, np)

	    } # end of if 'priorParams' is not provided by user.

	} # end of if 'priorFun' is not provided by user.

	out$priorFun <- priorFun
	out$priorParams <- priorParams

	if(is.null(proposalFun)) {

	    proposalFun <- "fevdProposalDefault"
	    if(is.null(proposalParams)) proposalParams <- list(sd=rep(0.1,np))

	} # end of if 'proposalFun' not supplied by user stmts.

	out$proposalFun <- proposalFun
	out$proposalParams <- proposalParams

	# Eric 05/30/14
	chain.info <- matrix(NA, iter, np + 2)
	colnames(chain.info) <- c(loc.names, sc.names, sh.names, "loglik", "prior")
	chain.info[2:iter, 1:np] <- 0

	res <- matrix(NA, iter, np+1)
	res[,np+1] <- 0
	colnames(res) <- c(loc.names, sc.names, sh.names, "new")

	if(nloc > 0) res[1,1:nloc] <- initial$location

	if(use.phi) res[1,(nloc+1):(nloc+nsc)] <- initial$log.scale
	else res[1,(nloc+1):(nloc+nsc)] <- initial$scale

	res[1,(nloc+nsc+1):np] <- initial$shape

	theta.i <- res[1,1:np]

	if(verbose) {
	    cat("\n", "Finding log-Likelihood of initial parameter values:\n")
	    print(theta.i)
	} 
	if(!missing(data)) ll.i <- -oevd(p=res[1,], o=out, des=designs, x=x, data=data, u=threshold, span=span, npy=npy, phi=use.phi, blocks=blocks) # CJP
	else ll.i <- -oevd(p=res[1,], o=out, des=designs, x=x, u=threshold, span=span, npy=npy, phi=use.phi, blocks=blocks)

	if(verbose) cat("\n", "Finding prior df value of initial parameter values.\n")
	p.i <- do.call(priorFun, c(list(theta=theta.i), priorParams))

	chain.info[1, np + 1] <- ll.i
	chain.info[1, np + 2] <- p.i

	if(verbose) cat("\n", "Beginning the MCMC iterations (", iter, " total iterations)\n")
	for(i in 2:iter) {

	    if(verbose && i <= 10) cat(i, " ")
	    if(verbose && i %% 100 == 0) cat(i, " ")
	    ord <- sample(1:np, np)
	    theta.star <- theta.i
	    acc <- 0

	    for(j in ord) {

		# par.star <- do.call(proposalFun, c(list(p=theta.i[j], ind=j), proposalParams))

		# theta.star[j] <- par.star

		par.star <- do.call(proposalFun, c(list(p=theta.i, ind=j), proposalParams))
		theta.star[j] <- par.star[j]

		if(!missing(data)) ll.star <- -oevd(p=theta.star, o=out, des=designs, x=x, data=data, u=threshold, span=span, npy=npy, phi=use.phi, blocks=blocks) # CJP
        	else ll.star <- -oevd(p=theta.star, o=out, des=designs, x=x, u=threshold, span=span, npy=npy, phi=use.phi, blocks=blocks)

		prior.star <- do.call(priorFun, c(list(theta=theta.star), priorParams))

		look <- will.accept(ll.i=ll.i, prior.i=p.i, ll.star=ll.star, prior.star=prior.star, log=TRUE)

		if(look$accept) {

		    p.i <- prior.star
		    ll.i <- ll.star
		    theta.i <- theta.star
		    acc <- acc+1
		    chain.info[i, j] <- 1

		}

	    } # end of inner 'j' loop.

	    res[i,] <- c(theta.i, acc)
	    chain.info[i, (np+1):(np+2)] <- c(ll.i, p.i)

	} # end of for 'i' loop.

	if(verbose) cat("\n", "Finished MCMC iterations.\n")
	out$results <- res
	out$chain.info <- chain.info

    } else stop("fevd: invalid method argument.")

    out$initial.results <- inout

    if(verbose) print(Sys.time() - begin.tiid)
    class(out) <- "fevd"

    return(out)

} # end of 'fevd' function.

levd <- function(x, threshold, location, scale, shape,
                type=c("GEV", "GP", "PP", "Gumbel", "Weibull", "Frechet", "Exponential", "Beta", "Pareto"),
                log=TRUE, negative=TRUE, span, npy=365.25, infval=Inf, weights=1, blocks=NULL) {

# CJP added blocks argument
    type <- match.arg(type)
    type <- tolower(type)
    n <- length(x)

    w <- weights   
    if(length(w) == 1) w <- rep(w, n)
 
    # if(type=="pp") {
# 	if(!missing(span)) {
# 	    if((length(threshold)==1) && (length(location)==1) && (length(scale)==1)) shortcut <- TRUE
# 	    else shortcut <- FALSE
 #        } else shortcut <- FALSE
  #   }
    # CJP: if we do implement 'shortcut', can't we check the const.thresh, const.location, etc. variables? I guess the issue is that we need to pass them in to levd()
    shortcut <- FALSE

    if(!missing(threshold)) if(length(threshold)==1) threshold <- rep(threshold, n)
    if(missing(threshold) && missing(location)) stop("levd: must supply one of threshold or location.")

    if(missing(location) && !missing(threshold) && is.element(type, c("gev", "weibull", "gumbel", "frechet", "pp"))) location <- threshold
    if(!missing(location)) if(length(location)==1) location <- rep(location, n)

    if(missing(scale)) stop("levd: must supply a scale argument.")
    else if(length(scale)==1) scale <- rep(scale, n)

    if(!missing(shape)) {
	if(length(shape)==1) shape <- rep(shape, n)
	if(any(is.na(shape))) {
	    if(!negative) return(-infval)
            else return(infval)
	}
	zero.shape <- shape==0
    	if(all(zero.shape) && !is.element(type, c("gumbel", "exponential"))) {
	    if(type=="gev") type <- "gumbel"
	    else if(type=="gp") type <- "exponential"
	}
	if(any(zero.shape) && any(!zero.shape)) {
            warning("levd: some shape parameters are zero, and some not.  Re-setting those that are to -1e-10.")
            shape[zero.shape] <- -1e-10
        }
    } else shape <- 0
    # end of if '!missing(shape)' stmts.

    if(!is.null(blocks)){ # CJP
      if(any(is.na(blocks$shape))) {
        if(!negative) return(-infval)
        else return(infval)
      }
      blocks$zero.shape <- blocks$shape==0
      if(any(blocks$zero.shape) && any(!blocks$zero.shape)) {
        warning("levd: some shape parameters are zero, and some not.  Re-setting those that are to -1e-10.")
        blocks$shape[blocks$zero.shape] <- -1e-10
      }
    }

    if(any(scale <= 0)) {
	if(!negative) return(-infval)
	else return(infval)
    }

    if(!is.null(blocks) && any(blocks$scale <= 0)) {  # CJP
	if(!negative) return(-infval)
	else return(infval)
    }
    
# stop("levd: invalid scale parameter.  Must be positively valued.")
    if(is.element(type, c("gp","exponential","beta","pareto"))) {
	excess.id <- x > threshold
	threshold <- threshold[excess.id]
	scale <- scale[excess.id]
	shape <- shape[excess.id]
	y <- (x[excess.id] - threshold)/scale
	m <- length(y)
    }
    if(type=="gumbel") {
	y <- (x - location)/scale
	res <- sum(colSums(w * cbind(log(scale), y, exp(-y))))
    } else if(is.element(type, c("gev","weibull","frechet"))) {
	y <- (x - location)/scale
	y <- 1 + shape*y
	if(any(y <= 0)) {
	    if(!negative) return(-infval)
	    else return(infval)
	}
	if(type=="frechet") shape <- abs(shape)
	else if(type=="weibull") shape <- -abs(shape)
	res <- sum(w * cbind(log(scale), y^(-1/shape), log(y) * (1/shape + 1))) # CJP removed colSums() here
    } else if(type=="exponential") res <- sum(w[excess.id] * cbind(log(scale), y)) # CJP removed colSums() here
    else if(is.element(type,c("gp","beta","pareto"))) {
	y <- 1 + shape*y
	if(any(y <= 0)) {
	    if(!negative) return(-infval)
            else return(infval)
	}
	if(type=="pareto") shape <- abs(shape)
	else if(type=="beta") shape <- -abs(shape)
	res <- sum(w[excess.id] * cbind(log(scale), log(y) * (1/shape + 1))) # CJP removed colSums()
    } else if(type=="pp") {
	if(any(zero.shape)) {
	    warning("levd: some shape parameters are zero.  Re-setting them to -1e-10.")
	    shape[zero.shape] <- -1e-10
	}
        if(!is.null(blocks) && any(blocks$zero.shape)) { # CJP
          warning("levd: some shape parameters are zero.  Re-setting them to -1e-10.")
          blocks$shape[blocks$zero.shape] <- -1e-10
	}

	excess.id <- x > threshold
	y <- (x - location)/scale
	y <- 1 + shape * y

	if(!shortcut) z <- 1 + ((shape*(threshold - location))/scale) # CJP made this produce a vector not a matrix
	else z <- 1 + ((shape[1]*(threshold[1] - location[1]))/scale[1])

	if(any(y[excess.id] <= 0) || any(z < 0)) {
                if(!negative) return(-infval)
                else return(infval)
        }

        if(!is.null(blocks)){  # CJP
          blocks$z <- 1 + ((blocks$shape*(blocks$threshold - blocks$location))/blocks$scale)
          if(any(blocks$z < 0)) {
            if(!negative) return(-infval)
            else return(infval)
          }
        } 

	th.u <- threshold[excess.id]
	sc.u <- scale[excess.id]
	sh.u <- shape[excess.id]
	y.u <- y[excess.id]

        if(is.null(blocks)) {  # CJP: sum(w*vec1 + w*vec2) may be faster than cbind
          if(!shortcut) res <- sum(w[excess.id] * cbind(log(sc.u), log(y.u) * (1/sh.u + 1))) + sum(w * z^(-1/shape))/npy # CJP removed colSums
          else res <- sum(w[excess.id] * cbind(log(sc.u), log(y.u) * (1/sh.u + 1))) + span * z^(-1/shape[1]) # CJP removed colSums
        } else {
          res <- sum(w[excess.id] * cbind(log(sc.u), log(y.u) * (1/sh.u + 1))) # pp density of exceedances
          res <- res + blocks$nBlocks * mean(blocks$weights * (1 - blocks$proportionMissing) * blocks$z^(-1/blocks$shape)) # this works either if stationary and each element is a scalar or if one or more vectors, though at the moment, I think $z and $shape will always be passed in as vectors
        }
    } # end of which type stmts.
    if(!negative) res <- -res
    if(!log) res <- exp(res)
    if(is.na(res)) {
	if(!negative) return(-infval)
	else return(infval)
    }
    return(res)
} # end of 'levd' function.

parcov.fevd <- function(x) {
    cov.theta <- NULL
    theta.hat <- x$results$par
    theta.names <- names(theta.hat)
    if(is.element("log.scale",theta.names)) {
        theta.hat[theta.names=="log.scale"] <- exp(theta.hat[theta.names=="log.scale"])
        theta.names[theta.names=="log.scale"] <- "scale"
        names(theta.hat) <- theta.names
        phiU <- FALSE
    } else phiU <- x$par.models$log.scale
    y <- datagrabber(x)
    if(x$data.name[2] != "") {
	xdat <- y[,1]
	data <- y[,-1]
    } else {
	xdat <- y[,1]
	data <- NULL
    }
    designs <- setup.design(x)
    if(x$method != "GMLE") {

    hold <- try(suppressWarnings(optimHess(theta.hat, oevd, gr=grlevd, o=x, des=designs, x=xdat, data=data, u=x$threshold, npy=x$npy, phi=phiU, blocks=x$blocks)), silent=TRUE) # CJP

        if((class(hold) != "try-error") && all(!is.na(hold))) {

            cov.theta <- try(suppressWarnings(solve(hold)), silent=TRUE)
	    if(any(diag(cov.theta) <= 0)) re.do <- TRUE
	    else re.do <- FALSE 

        } else re.do <- TRUE

    } else re.do <- TRUE


    if(re.do) {

       hold <- try(optimHess(theta.hat, oevd, o=x, des=designs, x=xdat, data=data, u=x$threshold, npy=x$npy, phi=phiU, blocks=x$blocks), silent=TRUE) # CJP
       if(class(hold) != "try-error" && all(!is.na(hold))) {

           cov.theta <- try(solve(hold), silent=TRUE)
	   if(class(cov.theta) == "try-error" || any(diag(cov.theta) <= 0)) cov.theta <- NULL

	}

    } # end of if 're.do' stmts.

    return(cov.theta)

} # end of 'parcov.fevd' function.

distill.fevd <- function(x, ...) {

    if(x$method=="GMLE") newcl <- "fevd.mle"
    else newcl <- tolower(x$method)
    class(x) <- paste("fevd.", newcl, sep="")
    UseMethod("distill",x)

} # end of 'distill.fevd' function.

distill.fevd.lmoments <- function(x, ...) {

    return(x$results)

} # end of 'distill.fevd.lmoments' function.

distill.fevd.bayesian <- function(x, cov=TRUE, FUN="mean", burn.in=499, ...) {

    f <- match.fun(FUN)
    p <- x$results
    np <- dim(p)[2] - 1
    p <- p[,-(np + 1)]
    pnames <- colnames(p)

    if(is.element("log.scale", pnames)) {

	id <- pnames == "log.scale"
	p[,"log.scale"] <- exp(p[,"log.scale"])
	pnames[id] <- "scale"
	colnames(p) <- pnames

    }

    if(burn.in != 0) {

	n <- dim(p)[1] 
	if(missing(burn.in)) if(burn.in <= 2 * n - 1) burn.in <- floor(n/4)
	else if(burn.in <= n - 1) stop("distill: number of MCMC iterations too small compared to burn.in")
	p <- p[-(1:burn.in),]

    }

    if(FUN == "mean") out <- colMeans(p, na.rm = TRUE)
    else if(FUN == "postmode" || FUN == "mode") out <- postmode(x, burn.in = burn.in, ...)
    else out <- apply(p, 2, f, ...)

    if(cov) {

	cov.theta <- cov(p)
	out <- c(out, c(cov.theta))
	names(out) <- c(pnames, paste(pnames[rep(1:np,np)], pnames[rep(1:np,each=np)], sep="."))

    }

    return(out)

} # end of 'distill.fevd.bayesian' function.

distill.fevd.mle <- function(x, cov=TRUE, ...) {
    theta.hat <- x$results$par
    theta.names <- names(theta.hat)
    np <- length(theta.hat)

    if(is.element("log.scale",theta.names)) {
         theta.hat[theta.names=="log.scale"] <- exp(theta.hat[theta.names=="log.scale"])
         theta.names[theta.names=="log.scale"] <- "scale"
         names(theta.hat) <- theta.names
     }

    out <- c(theta.hat, x$results$value)
    if(cov) {
	cov.theta <- parcov.fevd(x)

        if(!is.null(cov.theta)) {
            out <- c(out, cov.theta)
            names(out) <- c(theta.names, "nllh", paste(theta.names[rep(1:np,np)], theta.names[rep(1:np,each=np)], sep="."))
        } else names(out) <- c(theta.names, "nllh")
    } else names(out) <- c(theta.names, "nllh")

    return(out)
} # end of 'distill.fevd.mle' function.

summary.fevd <- function(object, ...) {

    if(object$method == "GMLE") newcl <- "mle"
    else newcl <- tolower(object$method)
    class(object) <- paste("fevd.", newcl, sep="")
    UseMethod("summary", object)

} # end of 'summary.fevd' function.

summary.fevd.lmoments <- function(object, ...) {
    x <- object
    a <- list(...)
    if(!is.null(a$silent)) silent <- a$silent
    else silent <- FALSE

    theta.hat <- x$results

    if(!silent) {
	print(x$call)
	if(x$data.name[2]=="") print(paste(x$type, " Fitted to ", x$data.name[1], " using L-moments estimation."))
	else  print(paste(x$type, " Fitted to ", all.vars(x$x.fun), " of ", x$data.name[2], "data frame, using L-moments estimation."))
	print(theta.hat)
    }
    invisible(theta.hat)
} # end of 'summary.fevd.lmoments' function.

summary.fevd.bayesian <- function(object, FUN="mean", burn.in=499, ...) {

    x <- object
    a <- list(...)
    if(!is.null(a$silent)) silent <- a$silent
    else silent <- FALSE

    np <- dim(x$results)[2] - 1

    if(!is.na(burn.in) && !is.null(burn.in) && burn.in > 0) {

	ll <- x$chain.info[-(1:burn.in), np + 1]
	mth <- colMeans(x$results[-(1:burn.in), -(np + 1)], na.rm = TRUE)

    } else {

	ll <- x$chain.info[, np + 1]
	mth <- colMeans(x$results[, -(np + 1)], na.rm = TRUE)

    } # end of if else 'burn.in' stmts.

    const <- is.fixedfevd(x)

    if(!silent) {
        cat("\n")
        print(x$call)
        cat("\n")
        print(paste("Estimation Method used: ", x$method, sep=""))
        cat("\n")
    }

    out <- list()
    out$par <- ci(x, type="parameter", FUN=FUN, burn.in=burn.in)
    # if(const) out$rl <- ci(x, return.period=c(2,20,100), FUN=FUN, burn.in=burn.in)
    pdim <- dim(x$results)
    # out$acceptance.rate <- colMeans(x$results[2:pdim[1], -pdim[2]] != x$results[1:(pdim[1] - 1), -pdim[2] ])
    out$acceptance.rate <- colMeans(x$chain.info[, 2:(pdim[2] - 1)], na.rm = TRUE)

    if(!silent) {
	cat("\n", "Acceptance Rates:\n")
	print(out$acceptance.rate)
	print(out$par)
	# if(const) print(out$rl)
    }

    cov.theta <- distill(x, FUN=FUN, burn.in=burn.in)
    cov.theta <- matrix(cov.theta[(np + 1):length(cov.theta)], np, np)
    colnames(cov.theta) <- rownames(cov.theta) <- colnames(x$results)[1:np]
    if(!silent) {
	cat("\n", "Estimated parameter covariance matrix.\n")
	print(cov.theta)
    }
    out$cov.theta <- cov.theta

    out$se.theta <- sqrt(diag(cov.theta))

    designs <- setup.design(x)
    y <- c(datagrabber(x, cov.data=FALSE))
    if(x$data.name[2] == "") data <- NULL
    else data <- datagrabber(x, response=FALSE)

    if(is.element(x$type, c("GEV", "Gumbel", "Weibull", "Frechet"))) n <- x$n
    else n <- sum(y > x$threshold)

    nllh <- oevd(p=out$par, o=x, des=designs, x=y, data=data, u=x$threshold, span=x$span, npy=x$npy, phi=x$par.models$log.scale, blocks=x$blocks)  # CJP
    out$nllh <- nllh
    # out$AIC <- 2 * nllh + 2 * np
    if(FUN == "postmode" || FUN == "mode") out$BIC <- 2 * nllh + np * log(n)
    Dev <- -2 * ll
    Dbar <- mean(Dev, na.rm = TRUE)
    Dmth <- -2 * oevd(p=mth, o=x, des=designs, x=y, data=data, u=x$threshold, span=x$span, npy=x$npy, phi=x$par.models$log.scale, blocks=x$blocks)

    pd <- Dbar - Dmth

    out$DIC <- Dmth + 2 * pd

    if(!silent) {

# 	cat("\n", "AIC =", out$AIC, "\n")
 #        cat("\n", "BIC =", out$BIC, "\n")
	cat("\n", "DIC = ", out$DIC, "\n")

    }

    invisible(out)

} # end of 'summary.fevd.bayesian' function.

summary.fevd.mle <- function(object, ...) {
    x <- object
    a <- list(...)

    out <- list()

    cov.theta <- se.theta <- NULL
    if(!is.null(a$silent)) silent <- a$silent
    else silent <- FALSE

    if(!silent) {
        cat("\n")
        print(x$call)
        cat("\n")
        print(paste("Estimation Method used: ", x$method, sep=""))
        cat("\n")
    }

    if(!silent) cat("\n", "Negative Log-Likelihood Value: ", x$results$value, "\n\n")
    theta.hat <- x$results$par
    theta.names <- names(theta.hat)
    if(is.element("log.scale",theta.names)) {
	theta.hat[theta.names=="log.scale"] <- exp(theta.hat[theta.names=="log.scale"])
	theta.names[theta.names=="log.scale"] <- "scale"
	names(theta.hat) <- theta.names
	phiU <- FALSE
    } else phiU <- x$par.models$log.scale

    out$par <- theta.hat
    np <- length(theta.hat)

    designs <- setup.design(x)

    if(!is.null(x$data.pointer)) xdat <- get(x$data.pointer)
    else xdat <- x$x

    cov.theta <- parcov.fevd(x)
    out$cov.theta <- cov.theta
    if(!is.null(cov.theta)) {
	se.theta <- sqrt(diag(cov.theta))
	names(se.theta) <- theta.names
	out$se.theta <- se.theta
    }

    if(!silent) {
	cat("\n", "Estimated parameters:\n")
        print(theta.hat)
    }

    if(!is.null(se.theta)) {
	if(!silent) {
	    cat("\n", "Standard Error Estimates:\n")
	    print(se.theta)
	}
	theta.hat <- rbind(theta.hat, se.theta)
	if(is.matrix(theta.hat) && dim(theta.hat)[1]==2) rownames(theta.hat) <- c("Estimate", "Std. Error")
	if(is.matrix(theta.hat)) colnames(theta.hat) <- theta.names
	else names(theta.hat) <- theta.names
	if(!silent && !is.null(cov.theta)) {
	    cat("\n", "Estimated parameter covariance matrix.\n")
	    print(cov.theta)
	}
    }

    nllh <- x$results$value
    out$nllh <- nllh

    if(is.element(x$type, c("GEV","Gumbel","Weibull","Frechet"))) n <- x$n
    else {
	y <- c(datagrabber(x, cov.data=FALSE))
	n <- sum(y > x$threshold)
    }
    out$AIC <- 2 * nllh + 2 * np
    out$BIC <- 2 * nllh + np * log(n)

    if(!silent) {
	cat("\n", "AIC =", out$AIC, "\n")
	cat("\n", "BIC =", out$BIC, "\n")
    }

    invisible(out)
} # end of 'summary.fevd.mle' function.

print.fevd <- function(x, ...) {
   summary(x)
   invisible()
} # end of 'print.fevd' function.

grlevd <- function(p, o, des, x, data, u=NULL, span, npy, phi=FALSE, blocks=NULL) {
    ##
    ## Function to estimate the likelihood
    ## parameter gradients.
    ##
    npar <- length(p)
    type <- o$type
    n <- length(x)

    w <- o$weights
    if(is.null(w)) w <- 1
    if(length(w) == 1) w <- rep(w, n)

    shortcut <- FALSE
    dpar <- matrix(NA, n, npar)

    if(!is.null(blocks)) # CJP
      blocks$dpar <- matrix(NA, blocks$nBlocks, npar)

    if(is.element(type, c("PP","GEV","Gumbel","Weibull","Frechet"))) {
	X.loc <- des$X.loc
	nloc <- ncol(X.loc)
	loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc)
	dpar[,1:nloc] <- X.loc

        if(!is.null(blocks)) {  # CJP
          blocks$loc <- rowSums(matrix(p[1:nloc], blocks$nBlocks, nloc, byrow=TRUE)*blocks$designs$X.loc)
          blocks$dpar[,1:nloc] <- blocks$designs$X.loc
        }
    } else {
	nloc <- 0
	loc <- NULL
        if(!is.null(blocks)) # CJP
          blocks$loc <- NULL
    }

    X.sc <- des$X.sc
    nsc <- ncol(X.sc)
    if(phi) scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc))
    else scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)

    dpar[,(nloc+1):(nloc+nsc)] <- X.sc
    
    if(phi) phi.gradTerm <- scale else phi.gradTerm <- 1 # CJP2

    if(!is.null(blocks)) { # CJP
      if(phi) {
        blocks$scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE)*blocks$designs$X.sc))
        blocks$phi.gradTerm <- blocks$scale
      } else {
        blocks$scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE)*blocks$designs$X.sc)
        blocks$phi.gradTerm <- 1
      }
      blocks$dpar[,(nloc+1):(nloc+nsc)] <- blocks$designs$X.sc
    }

    if(!is.element(type, c("Gumbel","Exponential"))) {
	X.sh <- des$X.sh
	nsh <- ncol(X.sh)
	shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh)
	dpar[,(nloc+nsc+1):(nloc+nsc+nsh)] <- X.sh

        if(!is.null(blocks)) {  # CJP
          blocks$shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], blocks$nBlocks, nsh, byrow=TRUE)*blocks$designs$X.sh)
          blocks$dpar[,(nloc+nsc+1):(nloc+nsc+nsh)] <- blocks$designs$X.sh
        }
    } else {
	nsh <- 0
	shape <- NULL
        if(!is.null(blocks))  # CJP
          blocks$shape <- NULL
    }

    if(is.element(type, c("Frechet","Pareto"))) shape <- abs(shape)
    else if(is.element(type, c("Weibull","Beta"))) shape <- -abs(shape)

    
    res <- p+NA
    if(is.element(type, c("GEV","Weibull","Frechet"))) {

        sh.sc <- shape/scale
        sh.1 <- 1/shape + 1
        x.fmu <- x - loc
        z <- 1 + sh.sc*x.fmu

        for(i in 1:npar) {  # CJP: don't need colSums and sum(w*vec1 + w*vec2 + w*vec3) may be faster than cbind
            if(i <= nloc) res[i] <- sum(colSums(w * cbind(-sh.1 * (sh.sc) * dpar[,i]/z, (dpar[,i]/scale)*z^(-sh.1))))
            else if(i <= nloc + nsc) res[i] <- sum(colSums(w * cbind(phi.gradTerm*dpar[,i]/scale, -phi.gradTerm * sh.1 * (shape/(scale^2))*x.fmu*dpar[,i]/z, (z^(-sh.1)/(scale^2)) * dpar[,i] * x.fmu * phi.gradTerm)))  # CJP2
            else res[i] <- sum(colSums(w * cbind(sh.1 * (dpar[,i]/scale) * x.fmu/z, -dpar[,i]/(shape^2)*log(z),
                                    z^(-1/shape) * ((dpar[,i]/(shape^2))*log(z) - (dpar[,i]/scale)*x.fmu/(z*shape)))))
        } # end of for 'i' loop.

    } else if(type == "Gumbel") {
        x.fmu <- x - loc
        z <- x.fmu/scale
        for(i in 1:npar) {
          if(i <= nloc) res[i] <- sum(colSums(w * cbind(-dpar[,i]/scale, (dpar[,i]/scale)*exp(-z))))
          else res[i] <- sum(colSums(w * cbind(phi.gradTerm*dpar[,i]/scale, -((x.fmu/(scale^2))*dpar[,i])*phi.gradTerm, phi.gradTerm*(dpar[,i]/(scale^2))*x.fmu*exp(-z))))  # CJP2 - Eric, please check this
       } # end of for 'i' loop.
   } else if(is.element(type, c("GP","Beta","Pareto"))) {
       excess.id <- x > u
       x <- x[excess.id]
       scale <- scale[excess.id]
       shape <- shape[excess.id]
       u <- u[excess.id]
       dparu <- dpar[excess.id,]
       x.fu <- x - u
       z <- 1 + (shape/scale)*x.fu
       sh.1 <- 1/shape + 1
       if(phi) phi.gradTermExc <- scale else phi.gradTermExc <- 1   # CJP2
       for(i in 1:npar) {
           if(i <= nsc) res[i] <- sum(colSums(w[excess.id] * cbind(phi.gradTermExc * dparu[,i]/scale, - ((phi.gradTermExc*sh.1*dparu[,i]*shape*x.fu/(scale^2))/z)))) # CJP2 - Eric, please check this
           else res[i] <- sum(colSums(w[excess.id] * cbind((-dparu[,i]/(shape^2))*log(z), sh.1*((dparu[,i]/scale)*x.fu)/z)))
       } # end of for 'i' loop.
   } else if(type=="Exponential") {
       excess.id <- x > u
       if(!is.matrix(dpar)) dpar <- cbind(dpar)
       dparu <- dpar[excess.id,]
       if(!is.matrix(dparu)) dparu <- cbind(dparu)
       scale <- scale[excess.id]
       xu <- x[excess.id] - u[excess.id]
       if(phi) phi.gradTermExc <- scale else phi.gradTermExc <- 1   # CJP2
       for(i in 1:npar) res[i] <- sum(colSums(w[excess.id] * cbind(phi.gradTermExc*dparu[,i]/scale, - phi.gradTermExc * xu * dparu[,i]/(scale^2)))) # CJP2 - Eric, please check this
   } else if(type=="PP") {
       excess.id <- x > u
       mu <- loc[excess.id]
       sc <- scale[excess.id]
       sh <- shape[excess.id]
       shu.scu <- sh/sc
       uu <- u[excess.id]
       xu <- x[excess.id]
       dparu <- dpar[excess.id,]
       x.fmu <- xu - mu
       if(!shortcut) u.l <- u - loc
       else u.l <- u[1] - loc[1]
       zu <- 1 + shu.scu*x.fmu
       if(!shortcut) Z <- matrix(1 + (shape/scale)*u.l, ncol=1)
       else Z <- 1 + (shape[1]/scale[1])*u.l
       shu.1 <- 1/sh + 1
       if(!shortcut) sh.1 <- 1/shape + 1
       else sh.1 <- 1/shape[1] + 1
       zpsh1 <- Z^(-sh.1)
       if(phi) phi.gradTermExc <- sc else phi.gradTermExc <- 1   # CJP2

       if(!is.null(blocks)) {
         blocks$u.l <- blocks$threshold - blocks$loc         
         blocks$Z <- 1 + (blocks$shape/blocks$scale)*blocks$u.l
         blocks$sh.1 <- 1/blocks$shape + 1
         blocks$zpsh1 <- blocks$Z^(-blocks$sh.1)
       }
       for(i in 1:npar) { # CJP: I changed sum(colSums(...)) here to sum()
	    if(!shortcut) { # CJP: same comment about sum(w*vec1 + w*vec2)
              if(is.null(blocks)) { # CJP
                if(i <= nloc) res[i] <- sum(w[excess.id] * cbind(-shu.1*dparu[,i]*shu.scu/zu)) + colSums(w * zpsh1 * dpar[,i]/scale)/npy # CJP: I think this can be sum not colSums
                else if(i <= nloc + nsc) res[i] <- sum(w[excess.id] * cbind(phi.gradTermExc * dparu[,i]/sc, -(shu.1*dparu[,i]*sh*x.fmu/(zu*sc^2))*phi.gradTermExc)) + colSums(w * phi.gradTerm * dpar[,i]*u.l*zpsh1/(scale^2))/npy  # CJP: I think this can be sum not colSums; CJP2 - added phi.gradTerm
                else res[i] <- sum(w[excess.id] * cbind(-(dparu[,i]/(sh^2))*log(zu), (shu.1*(dparu[,i]/sc)*x.fmu)/zu)) +
                                            colSums(w * (exp(-log(Z)/shape)*((dpar[,i]/(shape^2))*log(Z) - (dpar[,i]/(shape*scale)*u.l)/Z)))/npy
              } else {  # w/ blocks
                if(i <= nloc) res[i] <- sum(w[excess.id] * cbind(-shu.1*dparu[,i]*shu.scu/zu)) +
                  blocks$nBlocks * mean(blocks$weights * (1 - blocks$proportionMissing) * blocks$zpsh1 * blocks$dpar[,i]/blocks$scale)
                else if(i <= nloc + nsc) res[i] <- sum(w[excess.id] * cbind(phi.gradTermExc*dparu[,i]/sc, -(shu.1*dparu[,i]*sh*x.fmu/(zu*sc^2))*phi.gradTermExc)) +
                  blocks$nBlocks * mean(blocks$weights * (1 - blocks$proportionMissing) * blocks$phi.gradTerm * blocks$dpar[,i]*blocks$u.l*blocks$zpsh1/(blocks$scale^2)) # CJP2: added blocks$phi.gradTerm
                else res[i] <- sum(w[excess.id] * cbind(-(dparu[,i]/(sh^2))*log(zu), (shu.1*(dparu[,i]/sc)*x.fmu)/zu)) +
                  blocks$nBlocks * mean(blocks$weights * (1 - blocks$proportionMissing) * (exp(-log(blocks$Z)/blocks$shape)*((blocks$dpar[,i]/(blocks$shape^2))*log(blocks$Z) - (blocks$dpar[,i]/(blocks$shape*blocks$scale)*blocks$u.l)/blocks$Z)))
              }
	    } else {
	        if(i <= nloc) res[i] <- -sum(w[excess.id] * cbind(shu.1*dparu[,i]*shu.scu/zu)) + span * (zpsh1*dpar[1,i]/scale[1])
                else if(i <= nloc + nsc) res[i] <- sum(w[excess.id] * cbind(phi.gradTermExc*dparu[,i]/sc, -(shu.1*dparu[,i]*sh*x.fmu/(zu*sc^2))*phi.gradTermExc)) + span * (phi.gradTerm*dpar[1,i]*u.l*zpsh1/(scale[1]^2)) # CJP2: added phi.gradTerm
                else res[i] <- sum(w[excess.id] * cbind(-(dparu[,i]/(sh^2))*log(zu), shu.1*(dparu[,i]/sc)*x.fmu/zu)) +
                                            span * (exp(-log(Z)/shape[1])*((dpar[1,i]/(shape[1]^2))*log(Z) - (dpar[1,i]/(shape[1]*scale[1])*u.l)/Z))
	    }
       } # end of for 'i' loop.
    } # end of which gradient likelihood to use stmts.
    return(res)
} # end of 'grlevd' function.

grlevdTracer <- function(x, p=c(0,1,0), which.vary=1, p1.range, p2.range=NULL, threshold=NULL,
			    threshold.fun=~1, location.fun=~1, scale.fun=~1, shape.fun=~1, data, phi=FALSE,
			    type=c("GEV","GP","PP","Gumbel","Weibull","Frechet","Exponential","Beta","Pareto"),
			    nint=100, npy=365.25, weights = 1, blocks=NULL, na.action=na.fail, par1.name="", par2.name="", zl.lh=NULL, 
			    plot=TRUE, ...) {
# CJP note: this function now assumes that blocks contains the blocks data and design matrices already constructed within fevd and passed out as the $blocks element of the output from fevd - I didn't want to repeat all of that code here  

  # CJP2 : changed 'u' to 'threshold' and 'threshold' to 'threshold.fun'
  type <- match.arg(type)
    id <- which.vary
    data.name <- c(deparse(substitute(x)), deparse(substitute(data)))

    if(!is.function(na.action)) na.action <- get(na.action)

    if(!missing(data)) {
        if(length(x) != nrow(data)) stop("fevd: data must have same number of rows as the length of x.")
        tmp <- cbind(x, data)
        tmp <- na.action(tmp)
        x <- tmp[,1]
        data <- tmp[,-1]
	n <- length(x)
    } else {
	x <- na.action(x)
	n <- length(x)
    }

    w <- weights
    if(length(w) == 1) w <- rep(w, n)

    if(!is.null(threshold)) {
	const.thresh <- check.constant(threshold.fun)
	if(!missing(data)) X.u <- setup.design(x=threshold.fun, data=data, n=n, const=const.thresh, dname="threshold")
	else X.u <- setup.design(x=threshold.fun, n=n, const=const.thresh, dname="threshold")
	threshold <- rowSums(threshold*X.u)

    } # end of if 'threshold (u)' stmts.

    p1 <- seq(p1.range[1],p1.range[2],,nint)
    if(!is.null(p2.range)) p2 <- seq(p2.range[1],p2.range[2],,nint)

    loc.exists <- is.element(type, c("GEV","PP","Gumbel","Weibull","Frechet"))
    sh.exists <- !is.element(type, c("Gumbel","Exponential"))

    designs <- list()
    if(loc.exists) {
	const.loc <- check.constant(location.fun)
	if(!missing(data)) X.loc <- setup.design(x=location.fun, data=data, n=n, const=const.loc, dname="location.fun")
	else X.loc <- setup.design(x=location.fun, n=n, const=const.loc, dname="location")
	designs$X.loc <- X.loc
	nloc <- ncol(X.loc)
      } else {
	nloc <- 0
	const.loc <- NULL
    }

    const.scale <- check.constant(scale.fun)
    if(!missing(data)) X.sc <- setup.design(x=scale.fun, data=data, n=n, const=const.scale, dname="scale.fun")
    else X.sc <- setup.design(x=scale.fun, n=n, const=const.scale, dname="scale.fun")
    designs$X.sc <- X.sc
    nsc <- ncol(X.sc)

    if(sh.exists) {
        const.shape <- check.constant(shape.fun)
	if(!missing(data)) X.sh <- setup.design(x=shape.fun, data=data, n=n, const=const.shape, dname="shape.fun")
        else X.sh <- setup.design(x=shape.fun, n=n, const=const.shape, dname="shape")
	designs$X.sh <- X.sh
	nsh <- ncol(X.sh)
	if(is.element(type, c("Frechet","Pareto"))) shape <- abs(shape)
	else if(is.element(type, c("Weibull","Beta"))) shape <- -abs(shape)

    } else {
	nsh <- 0
	const.shape <- NULL
	shape <- 0
    } # end of 'shape.range' stmts.

    obj <- list(type=type, const.loc=const.loc, const.scale=const.scale, const.shape=const.shape,
		par.models=list(location=location.fun, scale=scale.fun, shape=shape.fun, phi=phi, weights=w,
		term.names=list(location=all.vars(location.fun), scale=all.vars(scale.fun), shape=all.vars(shape.fun))))

    num.vary <- length(id)
    if(num.vary==1) {
	lh <- lhgr <- numeric(nint)+NA
	for(i in 1:nint) {
            p[id] <- p1[i]
	    p.tmp <- p
	    scale <- rowSums(matrix(p.tmp[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)
	    if(phi) scale <- exp(scale)

            if(!is.null(blocks)) {
              blocks$location <- rowSums(matrix(p[1:nloc], blocks$nBlocks, nloc, byrow=TRUE)*blocks$designs$X.loc)
              blocks$scale <- rowSums(matrix(p.tmp[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE)*blocks$designs$X.sc)
              if(phi) blocks$scale <- exp(blocks$scale)
              blocks$shape=rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], blocks$nBlocks, nsh, byrow=TRUE)*blocks$designs$X.sh)
            }
            
            if(loc.exists && sh.exists) lh[i] <- levd(x=x, threshold=threshold, location=rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc),
						    scale=scale, shape=rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh),
						    type=type, npy=npy, weights=w, blocks=blocks)
            else if(loc.exists) lh[i] <- levd(x=x, threshold=threshold, location=rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc),
						    scale=scale, type=type, npy=npy, weights=w, blocks=blocks)
            else if(sh.exists) lh[i] <- levd(x=x, threshold=threshold, scale=scale, shape=rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh),
						    type=type, npy=npy, weights=w)
	    else lh[i] <- levd(x=x, threshold=threshold, scale=scale, type=type, npy=npy, weights=w, blocks=blocks)
            if(!missing(data)) lhgr[i] <- grlevd(p=p.tmp, o=obj, des=designs, x=x, data=data, u=threshold, npy=npy, phi=phi, blocks=blocks)[id]
            else lhgr[i] <- grlevd(p=p.tmp, o=obj, des=designs, x=x, u=threshold, npy=npy, phi=phi, blocks=blocks)[id]
        } # end of for 'i' loop.
	res <- cbind(p1, lh, lhgr)
	colnames(res) <- c("par.values","nllh","nllh.gr")
	attributes(res)$data.name <- data.name
	attributes(res)$par.name <- par1.name
	attributes(res)$model.name <- type
    } else if(num.vary==2) {
	lh <- lhgr1 <- lhgr2 <- matrix(NA, nint, nint)
	for(i in 1:nint) {
	    p[id[1]] <- p1[i]
	    p.tmp <- p
	    if(const.scale && phi) p.tmp[nloc+1] <- exp(p.tmp[nloc+1])
	    for(j in 1:nint) {

	        p[id[2]] <- p2[j]

		p.tmp <- p
		scale <- rowSums(matrix(p.tmp[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)
                if(phi) scale <- exp(scale)

                if(!is.null(blocks)) {
                  blocks$location <- rowSums(matrix(p[1:nloc], blocks$nBlocks, nloc, byrow=TRUE)*blocks$designs$X.loc)
                  blocks$scale <- rowSums(matrix(p.tmp[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE)*blocks$designs$X.sc)  # fixed by CJP2 6/11/13 (was 'blocks$X.sc')
                  if(phi) blocks$scale <- exp(blocks$scale)
                  blocks$shape=rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], blocks$nBlocks, nsh, byrow=TRUE)*blocks$designs$X.sh)
                }

		if(loc.exists && sh.exists) lh[i,j] <- levd(x=x, threshold=threshold, location=rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc),
							    scale=scale, shape=rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh),
							    type=type, npy=npy, weights=w, blocks=blocks)
                else if(loc.exists) lh[i,j] <- levd(x=x, threshold=threshold, location=rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc),
							    scale=scale, type=type, npy=npy, weights=w, blocks=blocks)
                else if(sh.exists) lh[i,j] <- levd(x=x, threshold=threshold, scale=scale, shape=rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh),
							    type=type, npy=npy, weights=w)
		else lh[i,j] <- levd(x=x, threshold=threshold, scale=scale, type=type, npy=npy, weights=w, blocks=blocks)

                if(!missing(data)) hold <- grlevd(p=p.tmp, o=obj, des=designs, x=x, data=data, u=threshold, npy=npy, blocks=blocks)[id]
                else hold <- grlevd(p=p.tmp, o=obj, des=designs, x=x, u=threshold, npy=npy, blocks=blocks)[id]
		lhgr1[i,j] <- hold[1]
		lhgr2[i,j] <- hold[2]
	    } # end of for 'j' loop.
	} # end of for 'i' loop.
	par.values <- list(p1,p2)
	names(par.values) <- c(paste("gr1.", par1.name, sep=""), paste("gr2.", par2.name, sep=""))
	res <- list(nllh=lh, nllh.gr1=lhgr1, nllh.gr2=lhgr2, p1=p1, p2=p2, call=match.call(), data.name=data.name, par.name=c(par1.name, par2.name), model=type)
	parvar <- cbind(rep(p1, nint), rep(p2, each=nint))
    } # end of 1 or two varied parameters stmts.
    class(res) <- "grlevdTracer"
    if(plot) plot(res, type="both", which.gr=1:2, ...)
    invisible(res)
} # end of 'grlevdTracer' function.

plot.grlevdTracer <- function(x, type=c("likelihood", "gradient", "both"), which.gr=1, log=TRUE, negative=TRUE, ...) {
    type <- match.arg(type)
    a <- list(...)

    if(is.element(type, c("likelihood","both")) && is.null(a$ylab)) {
	if(log && negative) nllh.ylab <- "Negative Log-Likelihood"
	else if(log) nllh.ylab <- "Log-Likelihood"
	else if(!negative) nllh.ylab <- "Likelihood"
	else nllh.ylab <- "Negative Likelihood"
    } 
    if(is.element(type, c("gradient", "both"))) gr.ylab <- "Negative Log-Likelihood Gradient"

    if(is.matrix(x)) {

	tmp <- attributes(x)
	dname <- tmp$data.name[1]
	pname <- tmp$par.name
	mname <- tmp$model.name

	if(is.element(type, c("likelihood","both"))) {
	    X <- x[,"nllh"]
	    if(!any(is.finite(X))) {
		skip <- TRUE
	    } else {
		skip <- FALSE
	        if(!log && !negative) X <- exp(-X)
	        else if(!negative) X <- -X
	        else if(!log) X <- -exp(-X)
	    }
	}

        if(is.null(a$ylim)) {
	    if(is.element(type, c("likelihood","both"))) yl <- range(X, finite=TRUE)
	    if(type=="both") par(mfrow=c(2,1))
	    if(is.element(type, c("likelihood","both"))) {
		if(skip) plot(1:length(X), rep(0, length(X)), type="n", yaxt="n", main="No finite negative log-likelihood values for this range")
		else {
		    if(is.null(a$ylab) && is.null(a$main)) {
		        if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", ylab=nllh.ylab, main=paste(mname, dname, sep="\n"), ylim=yl, ...)
		        else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, ylab=nllh.ylab, main=paste(mname, dname, sep="\n"), ylim=yl, ...)
		        else plot(x[,"par.values"], X, type="l", ylab=nllh.ylab, main=paste(mname, dname, sep="\n"), ylim=yl, ...)
		    } else if(is.null(a$ylab)) {
		        if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", ylab=nllh.ylab, ylim=yl, ...)
		        else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, ylab=nllh.ylab, ylim=yl, ...)
                        else plot(x[,"par.values"], X, type="l", ylab=nllh.ylab, ylim=yl, ...)
	            } else if(is.null(a$main)) {
		        if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", main=paste(mname, dname, sep="\n"), ylim=yl, ...)
                        else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, main=paste(mname, dname, sep="\n"), ylim=yl, ...)
                        else plot(x[,"par.values"], X, type="l", ylim=yl, main=paste(mname, dname, sep="\n"), ...)
		    } else {
		        if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", ylim=yl, ...)
                        else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, ylim=yl, ...)
                        else plot(x[,"par.values"], X, type="l", ylim=yl, ...)
		    }
	        } # end of if else skip stmts.
	    } # end of if type is likelihood or both stmts.

            if(is.element(type, c("gradient", "both"))) {
		if(!any(is.finite(x[,"nllh.gr"]))) {
		    plot(1:length(x[,"nllh.gr"]), rep(0, length(x[,"nllh.gr"])), type="n", yaxt="n", main="No finite negative log-likelihood gradient values for this range")
		} else {

	            if(is.null(a$main) && type=="both") m1 <- ""
	            else if(is.null(a$main)) m1 <- paste(mname, dname, sep="\n")

	            if(is.null(a$ylab) && is.null(a$main)) {
	        	# gr.ylab <- paste(gr.ylab, " wrt ", pname, sep="")
	        	if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, ylab=gr.ylab, main=m1, ...)
	        	else plot(x[,"par.values"], x[,"nllh.gr"], type="l", ylab=gr.ylab, main=m1, ...)
	            } else if(is.null(a$ylab)) {
		        # gr.ylab <- paste(gr.ylab, " wrt ", pname, sep="")
                        if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, ylab=gr.ylab, ...)
                        else plot(x[,"par.values"], x[,"nllh.gr"], type="l", ylab=gr.ylab, ...)
	            } else if(is.null(a$main)) {
	        	if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, main=m1, ...)
                        else plot(x[,"par.values"], x[,"nllh.gr"], type="l", main=m1, ...)
	            } else {
	        	if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, ...)
                        else plot(x[,"par.values"], x[,"nllh.gr"], type="l", ...)
	            } # end of if ylab and main supplied or not stmts.
		    abline(h=0, lty=2)
	        } # end of if else no finite gradient values stmts.
	    } # end of if type is gradient or both stmts.

	} else {
            if(type=="both") par(mfrow=c(2,1))
            if(is.element(type, c("likelihood","both"))) {
                if(is.null(a$ylab) && is.null(a$main)) {
                    if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", ylab=nllh.ylab, main=paste(mname, dname, sep="\n"), ...)
                    else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, ylab=nllh.ylab, main=paste(mname, dname, sep="\n"), ...)
                    else plot(x[,"par.values"], X, type="l", ylab=nllh.ylab, main=paste(mname, dname, sep="\n"), ...)
                } else if(is.null(a$ylab)) {
                    if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", ylab=nllh.ylab, ...)
                    else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, ylab=nllh.ylab, ...)
                    else plot(x[,"par.values"], X, type="l", ylab=nllh.ylab, ...)
                } else {
                    if((type == "both") && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab="", ...)
                    else if(type == "likelihood" && is.null(a$xlab)) plot(x[,"par.values"], X, type="l", xlab=pname, ...)
                    else plot(x[,"par.values"], X, type="l", ...)
                }
            } # end of if type is likelihood or both stmts.

            if(is.element(type, c("gradient", "both"))) {
                if(is.null(a$main) && type=="both") m1 <- ""
                else if(is.null(a$main)) m1 <- paste(mname, dname, sep="\n")
                if(is.null(a$ylab) && is.null(a$main)) {
                    gr.ylab <- paste(gr.ylab, " wrt ", pname, sep="")
                    if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, ylab=gr.ylab, main=m1, ...)
                    else plot(x[,"par.values"], x[,"nllh.gr"], type="l", ylab=gr.ylab, main=m1, ...)
                } else if(is.null(a$ylab)) {
                    gr.ylab <- paste(gr.ylab, " wrt ", pname, sep="")
                    if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, ylab=gr.ylab, ...)
                    else plot(x[,"par.values"], x[,"nllh.gr"], type="l", ylab=gr.ylab, ...)
                } else if(is.null(a$main)) {
                    if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, main=m1, ...)
                    else plot(x[,"par.values"], x[,"nllh.gr"], type="l", main=m1, ...)
                } else {
                    if(is.null(a$xlab)) plot(x[,"par.values"], x[,"nllh.gr"], type="l", xlab=pname, ...)
                    else plot(x[,"par.values"], x[,"nllh.gr"], type="l", ...)
                } # end of if ylab and main supplied or not stmts.
		abline(h=0, lty=2)
            } # end of if type is gradient or both stmts.
	} # end of if ylim is passed or not stmts.

    } else {
	dname <- x$data.name[1]
        p1name <- x$par.name[1]
        p2name <- x$par.name[2]
        p1 <- x$p1
        p2 <- x$p2
	nint <- length(p1)
	if(is.element(type, c("likelihood","both"))) {
	    X <- x$nllh
	    if(!log && !negative) X <- exp(-X)
            else if(!negative) X <- -X
            else if(!log) X <- -exp(-X)
	    minloc <- cbind(rep(p1, nint), rep(p2, each=nint))[which.min(c(X)),]
	}

	if(!(which.gr[1]==1 || which.gr[1]==2 || all(which.gr==1:2))) stop("plot.grlevdTracer: invalid which.gr argument.  Must be 1, 2 or 1:2.")

	if(is.null(a$xlab)) xl <- p1name
	else xl <- a$xlab
	
	if(is.null(a$ylab)) yl <- p2name
	else yl <- a$ylab

	if(!is.null(a$zlim)) zl.lh <- zl.gr1 <- zl.gr2 <- a$zlim
        else {
            zl.lh <- range(c(X), finite=TRUE)
            zl.gr1 <- range(c(x$nllh.gr1), finite=TRUE)
            zl.gr2 <- range(c(x$nllh.gr2), finite=TRUE)
        }

	# if(is.null(a$col)) zcol <- tim.colors(64)
	# else zcol <- a$col

	if(is.null(a$horizontal)) hz <- TRUE
	else hz <- a$horizontal

	if(is.null(a$legend.mar)) legmar <- ifelse(hz, 3.1, 5.1)
	else legmar <- a$legend.mar

	if(type=="both") {
	    if(length(which.gr)==1) par(mfrow=c(1,2))
	    else par(mfrow=c(1,3))
	}

	if(is.null(a$main)) {
	    if(is.element(type, c("likelihood","both"))) m1 <- paste(nllh.ylab, "\n", x$model, ": ", x$data.name[1], sep="")
	    else m1 <- ""
	    if(length(which.gr)==1) m2 <- m3 <- paste(x$model, ": Negative Log-Likelihood\nGradient wrt ", x$par.name[which.gr], " (", x$data.name[1], ")", sep="")
            else {
	        m2 <- paste(x$model, ": Negative Log-Likelihood\nGradient wrt ", p1name, " (", x$data.name[1], ")", sep="")
                m3 <- paste(x$model, ": Negative Log-Likelihood\nGradient wrt ", p2name, " (", x$data.name[1], ")", sep="")
	    }
        } else m1 <- m2 <- m3 <- a$main

	if(is.element(type,c("likelihood","both"))) {
	    image(p1, p2, X, main=m1, xlab=xl, ylab=yl, zlim=zl.lh)
	    points(minloc[1], minloc[2], pch="x", col="darkorange")
	    # image.plot(p1, p2, X, zlim=zl.lh, horizontal=hz, legend.only=TRUE, legend.mar=legmar)
	} # end of if likelihood or both stmts.

	if(is.element(type, c("gradient","both"))) {
	    if(is.element(1, which.gr)) {
		image(p1, p2, x$nllh.gr1, main=m2, xlab=xl, ylab=yl, zlim=zl.gr1)
		# image.plot(p1, p2, x$nllh.gr1,  zlim=zl.gr1, horizontal=hz, legend.only=TRUE, legend.mar=legmar)
	    }

	    if(is.element(2, which.gr)) {
                image(p1, p2, x$nllh.gr2, main=m3, xlab=xl, ylab=yl, zlim=zl.gr2)
                # image.plot(p1, p2, x$nllh.gr2,  zlim=zl.gr2, horizontal=hz, legend.only=TRUE, legend.mar=legmar)
            }

	} # end of if type is gradient or both stmts.

    }# end of if one or two parameters varied stmts.
    invisible()
} # end of 'plot.grlevdTracer' function.

check.constant <- function(x) {
    # internal function to determine whether a parameter is to be modeled
    # as constant (TRUE) or not (FALSE).
    if(is.formula(x)) {
        if(length(as.character(x))==2 && as.character(x)[2] == "1") return(TRUE)
        else return(FALSE)
    } else {
        if(length(x)==1) return(TRUE)
        else if(all(x == x[1])) return(TRUE)
        else return(FALSE)
    }
} # end of 'check.constant' function.

setup.design <- function(x, ...) {
    UseMethod("setup.design", x)
} # end of 'setup.design' function.

setup.design.default <- function(x, ..., data, n, const, dname) {
    if(missing(data)) data <- NULL
    if(missing(const)) const <- check.constant(x)
    if(const || is.vector(x)) return(matrix(1, n, 1))
    else if(is.formula(x)) {
        if(is.null(data)) return(model.matrix(x))
        else return(model.matrix(x, data))
    } else stop(paste("setup.design: invalid ", dname, " parameter.", sep=""))
} # end of 'setup.design' function.

setup.design.fevd <- function(x, ...) {
    out <- list()
    n <- x$n
    mods <- x$par.models
    dname <- x$data.name[2]
    data <- datagrabber(x, response=FALSE)

    out$X.u <- setup.design(mods$threshold, data=data, n=n, dname="threshold")
    out$X.loc <- setup.design(mods$location, data=data, n=n, dname="location")
    out$X.sc <- setup.design(mods$scale, data=data, n=n, dname="scale")
    out$X.sh <- setup.design(mods$shape, data=data, n=n, dname="shape")
    return(out)
} # end of 'setup.design.fevd' function.

oevd <- function(p, o, des, x, data=NULL, u=NULL, span, npy, phi=FALSE, blocks=NULL) {
    ##
    ## Function to be optimized by 'optim'.
    ## Takes the parameter design matrices and
    ## calculates the resulting vector of parameters,
    ## then calls 'levd' to get the negative log-likelihood.
    type <- o$type
    n <- length(x)

    w <- o$weights
    if(length(w) == 1) w <- rep(w, n)

    if(!is.element(type, c("GP","Beta","Pareto","Exponential"))) {
        X.loc <- des$X.loc
        nloc <- ncol(X.loc)
        loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE)*X.loc)
    } else {
        nloc <- 0
        loc <- 0
    }

    X.sc <- des$X.sc
    nsc <- ncol(X.sc)
    if(phi) scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc))
    else scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE)*X.sc)

    if(!is.element(o$type, c("Gumbel","Exponential"))) {
        X.sh <- des$X.sh
        nsh <- ncol(X.sh)
        shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE)*X.sh)
    } else {
        nsh <- 0
        shape <- 0
    }

    if(!is.null(blocks)) {  # could check for const.loc, etc. and just do sum(p[1:nloc]*blocks$designs$X.param[1,])
      if(is.element(type, c("PP"))) {
        blocks$location <- rowSums(matrix(p[1:nloc], blocks$nBlocks, nloc, byrow=TRUE)*blocks$designs$X.loc)
      }
      if(phi) blocks$scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE)*blocks$designs$X.sc))
      else blocks$scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE)*blocks$designs$X.sc)
      blocks$shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], blocks$nBlocks, nsh, byrow=TRUE)*blocks$designs$X.sh)
    }

    if(missing(span)) res <- levd(x=x, threshold=u, location=loc, scale=scale, shape=shape, type=type, npy=npy, infval=1e10, weights=w, blocks=blocks)
    else res <- levd(x=x, threshold=u, location=loc, scale=scale, shape=shape, type=type, span=span, npy=npy, infval=1e10, weights=w, blocks=blocks)
    return(res)

} # end of 'oevd' function.

oevdgen <- function(p, o, des, x, data = NULL, u = NULL, span, npy, phi = FALSE, priorFun, priorParams = NULL, blocks = NULL) {

    ##
    ## Function to be optimized by 'optim'.
    ## Takes the parameter design matrices and
    ## calculates the resulting vector of parameters,
    ## then calls 'levd' to get the negative log-likelihood.

    type <- o$type
    n <- length(x)

    w <- o$weights
    if(length(w) == 1) w <- rep(w, n)

    if(!is.element(type, c("GP","Beta","Pareto","Exponential"))) {

        X.loc <- des$X.loc
        nloc <- ncol(X.loc)
        loc <- rowSums(matrix(p[1:nloc], n, nloc, byrow=TRUE) * X.loc)
	if(nloc == 1) lname <- "location"
	else lname <- paste("mu", 0:(nloc - 1), sep = "")

    } else {

        nloc <- 0
        loc <- 0
	lname <- NULL

    }

    X.sc <- des$X.sc
    nsc <- ncol(X.sc)
    if(phi) scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE) * X.sc))
    else scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], n, nsc, byrow=TRUE) * X.sc)
    if(nsc == 1) scname <- "scale"
    else scname <- paste("phi", 0:(nsc - 1), sep="")

    if(!is.null(des$X.sh)) {

        X.sh <- des$X.sh
        nsh <- ncol(X.sh)
        shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], n, nsh, byrow=TRUE) * X.sh)
	if(nsh == 1) shname <- "shape"
	else shname <- paste("xi", 0:(nsh - 1), sep = "")

    } else {

        nsh <- 0
        shape <- 0
	shname <- NULL

    }

    if(!is.null(blocks)) {

      blocks$location <- rowSums(matrix(p[1:nloc], blocks$nBlocks, nloc, byrow=TRUE) * blocks$designs$X.loc)
      if(phi) blocks$scale <- exp(rowSums(matrix(p[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE) * blocks$designs$X.sc))
      else blocks$scale <- rowSums(matrix(p[(nloc+1):(nloc+nsc)], blocks$nBlocks, nsc, byrow=TRUE) * blocks$designs$X.sc)
      blocks$shape <- rowSums(matrix(p[(nloc+nsc+1):(nloc+nsc+nsh)], blocks$nBlocks, nsh, byrow=TRUE) * blocks$designs$X.sh)

    }
    
    if(missing(span)) res <- levd(x=x, threshold=u, location=loc, scale=scale, shape=shape, type=type, npy=npy, infval=1e10, weights=w, blocks=blocks)
    else res <- levd(x=x, threshold=u, location=loc, scale=scale, shape=shape, type=type, span=span, npy=npy, infval=1e10, weights=w, blocks=blocks)

    # newp <- cbind(loc, scale, shape)
    # colnames(newp) <- c(lname, scname, shname)

    # res2 <- do.call(priorFun, c(list(x = newp), priorParams))
    res2 <- do.call(priorFun, c(list(x = p), priorParams))

    return(res + res2)

} # end of 'oevdgen' function.

shapePriorBeta <- function(x, p, q) {
    ##
    ## Function to calculate the GMLE prior for the shape parameter.
    ##

    # if(!is.element("shape", colnames(x))) stop("shapePriorBeta: invalid model for this prior penalty function.") 
    # xi <- 0.5 + x[,"shape"]

    xi <- x[ length(x) ]

    # res <- -sum(dbeta(x=xi, shape1=q, shape2=p, log=TRUE), na.rm=TRUE)

    res <- - dbeta(x = xi, shape1 = q, shape2 = p, log = TRUE)
    if(!is.finite(res)) res <- sign(res) * 1e16

    return(res)

} # end of 'shapePriorBeta' function.

initializer <- function(x, model, threshold, ...) {

    UseMethod("initializer", x)

} # end of 'initializer' function.

initializer.lmoments <- function(x, model, threshold, ..., npy, blocks=NULL) {

    if(is.element(model, c("gev", "weibull", "frechet", "gumbel"))) {

	lambda <- Lmoments(x)
	tau3 <- lambda[3] / lambda[2]
        co <- 2 / (3 + tau3) - log(2) / log(3)
        kappa <- 7.8590 * co + 2.9554 * co^2
        g <- gamma(1 + kappa)
        sigma <- (lambda[2] * kappa)/((1 - 2^(-kappa)) * g)
        mu <- lambda[1] - (sigma / kappa) * (1 - g)
        xi <- -kappa
	res <- c(mu, sigma, xi)
	names(res) <- c("location", "scale", "shape")

    } else if(is.element(model, c("gp", "exponential", "beta", "pareto", "pp"))) {

	n <- length(x)
	u <- threshold
	if(length(u)==1) u <- rep(u, n)
	id <- x > u
	lambda <- Lmoments(x[id] - u[id])
        tau2 <- lambda[2] / lambda[1]
        sigma <- lambda[1] * (1 / tau2 - 1)
        kappa <- 1 / tau2 - 2
        xi <- -kappa

        if(model=="pp") {

          if(is.null(blocks)) {

            rate <- npy * sum(id)/n

          } else {

            rate <- sum(id)/(blocks$nBlocks * mean(blocks$weights))

          }

          sigma <- exp(log(sigma) + xi*log(rate))
          mu <- mean(threshold) - (sigma/xi)*(rate^(-xi) - 1)
          res <- c(mu, sigma, xi, rate)
          names(res) <- c("location", "scale", "shape", "rate")

        } else if(model != "exponential") {

	    res <- c(sigma, xi)
	    names(res) <- c("scale", "shape")

	} else {

	    res <- sigma
	    names(res) <- c("scale")

	}

    } # end of 'which models to estimate' stmts.

    return(res)

} # end of 'initializer.lmoments' function.

initializer.moms <- function(x, model, threshold, ..., npy, blocks=NULL) {
    if(is.element(model, c("gev","gumbel","frechet","weibull"))) {
	m <- mean(x)
	s <- sqrt(6*var(x))/pi
	mu <- m - 0.57722*s
	sigma <- log(s)
	if(model=="weibull") xi <- -1e-8
	else if(is.element(model,c("gev","frechet"))) xi <- 1e-8
	else xi <- 0
	res <- c(mu, sigma, xi)
	names(res) <- c("location", "scale", "shape")
    } else if(is.element(model, c("pp", "gp", "exponential", "beta", "pareto"))) {
	n <- length(x)
	u <- threshold
	if(length(u)==1) u <- rep(u,n)
	id <- x > u
	x.u <- x[id] - u[id]
	m <- mean(x.u)
	sigma <- sqrt(var(x.u))
	if(model=="beta") xi <- -1e-8
	else if(is.element(model, c("pp","gp","beta"))) xi <- 1e-8
	else xi <- 0
	if(model=="pp") {
          if(is.null(blocks)) {
	    rate <- npy * sum(id)/n
          } else {
            rate <- sum(id)/(blocks$nBlocks * mean(blocks$weights))
          }
          sigma <- exp(log(sigma) + xi*rate)
          mu <- mean(u) - (sigma/xi)*(rate^(-xi) - 1)
          res <- c(mu, sigma, xi)
          names(res) <- c("location", "scale", "shape")
	} else if(model != "exponential") {
	    res <- c(sigma, xi)
	    names(res) <- c("scale", "shape")
	} else {
	    res <- sigma
	    names(res) <- "scale"
	}
    } else stop("initializer: invalid model argument.")
    return(res)
} # end of 'initializer.moms' function.

initializer.mle <- function(x, model, threshold, ..., data=NULL, u, u.fun, loc.fun, sc.fun, sh.fun, use.phi, type, span, time.units, period.basis, blocks=NULL) {

    attributes(x) <- NULL

    if(is.null(data)) {
	if(missing(span)) res <- fevd(x=x, threshold=u, use.phi=use.phi, type=type, time.units=time.units, period.basis=period.basis, blocks=blocks)
	else res <- fevd(x=x, threshold=u, use.phi=use.phi, type=type, span=span, time.units=time.units, period.basis=period.basis, blocks=blocks)
    } else {
	if(missing(span)) res <- fevd(x=x, data=data, threshold=u, threshold.fun=u.fun, location.fun=loc.fun, scale.fun=sc.fun, shape.fun=sh.fun,
				    use.phi=use.phi, type=type, time.units=time.units, period.basis=period.basis, blocks=blocks)
        else res <- fevd(x=x, data=data, threshold=u, threshold.fun=u.fun, location.fun=loc.fun, scale.fun=sc.fun, shape.fun=sh.fun,
				    use.phi=use.phi, type=type, span=span, time.units=time.units, period.basis=period.basis, blocks=blocks)
    }
    return(res$results$par)
} # end of 'initializer.bayesian' function.

will.accept <- function( ll.i, prior.i, ll.star, prior.star, log=TRUE) {
    ##
    ## Function to determine whether to accept a proposed parameter value or not.
    ##
    ## Arguments:
    ##
    ## 'll.i' numeric giving the log-likelihood value for the i-th step parameter combination.
    ## 'prior.i' numeric giving the (log) prior probability for the i-th step parameter combination.
    ## 'll.star', 'prior.star' are the log-likelihood and (log) prior values for the proposed parameter vector.
    ##
    ## Details: Calculates the ratio R=(L(theta*;data)p(theta*))/(L(theta.i;data)p(theta.i)), and determines
    ##   whether to accept or reject the proposed parameters with probability alpha = min(1, R).
    ##
    ## Value: list with components:
    ## 'accept' a logical saying whether to accept the proposed parameters (if TRUE) or not (if FALSE),
    ## 'ratio' numeric giving the value of R (see Details above),
    ## 'alpha' numeric giving the value of 'alpha' (see Details above),
    ## 'll.star' numeric giving the log-likelihood of the proposed parameter combination,
    ## 'prior.star' the prior probability for the proposed parameter combination.
    ##
    out <- list()
    ratio <- (ll.star + prior.star) - (ll.i + prior.i)
    if(!log) ratio <- exp( ratio)
    if( is.na( ratio)) {
         out$accept <- FALSE
         out$ratio <- out$alpha <- NA
         return( out)
    }

    if(!log) alpha <- min( 1, ratio)
    else alpha <- min(0, ratio)

    if(!log) {
 	if( runif( 1) < alpha) accept <- TRUE
 	else accept <- FALSE
    } else {
 	if(log(runif(1)) < alpha) accept <- TRUE
	else accept <- FALSE
    }
    out$accept <- accept
    out$ratio <- ratio
    out$alpha <- alpha
    return( out)
} # end of 'will.accept' function.

fevdPriorDefault <- function(theta, m, v, log=TRUE) {

    np <- length(theta)
    dfun <- function(th) dnorm(th[1], mean=th[2], sd=th[3], log=log)

    if(length(m)==1) m <- rep(m, np)
    else if(length(m) != np) stop("fevdPriorDefault: m must have length 1 or same as number of parameters.")

    if(length(v)==1) v <- rep(v, np)
    else if(length(v) != np) stop("fevdPriorDefault: v must have length 1 or same as number of parameters.")

    th <- cbind(theta, m, v)
    res <- apply(th, 1, dfun)
    if(log) return(sum(res))
    else return(prod(res))

} # end of 'fevdPriorDefault' function.

fevdProposalDefault <- function(p, ind, ...) {

    a <- list(...)

    if(!is.null(a$mean)) {

	#if(length(a$mean) > 1) m <- a$mean[ind]
	#else m <- a$mean
	m <- a$mean

    } else m <- 0

    if(!is.null(a$sd)) {

	# if(length(a$sd) > 1) s <- a$sd[ind]
	# else s <- a$sd 
	s <- a$sd

    } else s <- 0.1
    
    return(p + rnorm(length(p), mean=m, sd=s))

} # end of 'fevdProposalDefault' function.

fpois <- function(x, na.action = na.fail, ...) {

    UseMethod("fpois", x)

} # end of 'fpois' function.

fpois.default <- function(x, na.action = na.fail, ...) {

    theCall <- match.call()
    dname <- deparse(substitute(x))

    x <- na.action(x)

    n <- length(x)
    m <- mean(x)
    v <- var(x)

    PARAMETER <- c(m, v, n - 1)
    names(PARAMETER) <- c("mean", "variance", "degrees of freedom")

    CRITVAL <- (n - 1) * v / m
    names(CRITVAL) <- "Chi-square(n - 1)"
    PVAL <- pchisq(CRITVAL, df = n - 1, lower.tail = FALSE)

    structure(list(statistic = CRITVAL, parameter = PARAMETER, 
        alternative = "greater", p.value = PVAL, method = "Test for Equality of (Poisson) Mean and Variance", 
        data.name = dname, call = theCall), class = "htest")

} # end of 'fpois.default' function.

fpois.data.frame <- function(x, na.action = na.fail, ..., which.col = 1) {

    theCall <- match.call()
    dname <- deparse(substitute(x))

    if(length(which.col) > 1) stop("fpois: invalid which.col argument.  Length must be one.  To include covariates, use glm with poisson() family link.")
    x <- x[, which.col]
    x <- na.action(x)

    n <- length(x)
    m <- mean(x)
    v <- var(x)

    PARAMETER <- c(m, v, n - 1)
    names(PARAMETER) <- c("mean", "variance", "degrees of freedom")

    CRITVAL <- (n - 1) * v / m
    names(CRITVAL) <- "Chi-square(n - 1)"
    PVAL <- pchisq(CRITVAL, df = n - 1, lower.tail = FALSE)

    structure(list(statistic = CRITVAL, parameter = PARAMETER,
        alternative = "greater", p.value = PVAL, method = "Test for Equality of (Poisson) Mean and Variance",
        data.name = dname, call = theCall), class = "htest")

} # end of 'fpois.data.frame' function.

fpois.matrix <- function(x, na.action = na.fail, ..., which.col = 1) {


    theCall <- match.call()
    dname <- deparse(substitute(x))

    if(length(which.col) > 1) stop("fpois: invalid which.col argument.  Length must be one.  To include covariates, use glm with poisson() family link.")
    x <- x[, which.col]
    x <- na.action(x)

    n <- length(x)
    m <- mean(x)
    v <- var(x)

    PARAMETER <- c(m, v, n - 1)
    names(PARAMETER) <- c("mean", "variance", "degrees of freedom")

    CRITVAL <- (n - 1) * v / m
    names(CRITVAL) <- "Chi-square(n - 1)"
    PVAL <- pchisq(CRITVAL, df = n - 1, lower.tail = FALSE)

    structure(list(statistic = CRITVAL, parameter = PARAMETER,
        alternative = "greater", p.value = PVAL, method = "Test for Equality of (Poisson) Mean and Variance",
        data.name = dname, call = theCall), class = "htest")

} # end of 'fpois.matrix' function.

fpois.list <- function(x, na.action = na.fail, ..., which.component = 1) {


    theCall <- match.call()
    dname <- deparse(substitute(x))

    if(length(which.component) > 1) stop("fpois: invalid which.component argument.  Length must be one.  To include covariates, use glm with poisson() family link.")
    x <- x[[which.component]]
    x <- na.action(x)

    n <- length(x)
    m <- mean(x)
    v <- var(x)

    PARAMETER <- c(m, v, n - 1)
    names(PARAMETER) <- c("mean", "variance", "degrees of freedom")

    CRITVAL <- (n - 1) * v / m
    names(CRITVAL) <- "Chi-square(n - 1)"
    PVAL <- pchisq(CRITVAL, df = n - 1, lower.tail = FALSE)

    structure(list(statistic = CRITVAL, parameter = PARAMETER,
        alternative = "greater", p.value = PVAL, method = "Test for Equality of (Poisson) Mean and Variance",
        data.name = dname, call = theCall), class = "htest")

} # end of 'fpois.list' function.
