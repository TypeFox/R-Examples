bayesQR.single <- function(formula=NULL, data=NULL, quantile=0.5, alasso=FALSE, ndraw=NULL, keep=1, prior=NULL){

  # Create function for error message
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }

  # Extract input data
	if (is.null(formula)){
    pandterm("Formula has to be specified")
	}
	mf <- model.frame(formula=formula, data=data)
	X <- model.matrix(attr(mf,"terms"), data=mf)
	y <- model.response(mf)
	n <- length(y)
	nvar <- ncol(X)
	names <- colnames(X)

	# Check specified quantile
  if ((quantile<=0)|(quantile>=1)) {
    pandterm("Quantiles should be between zero and one")
  }

  # Check specified number of mcmc draws
  if (is.null(ndraw)) {
    pandterm("Number of mcmc draws (ndraw) has to be specified")
  }

  # Check specified number of mcmc draws
  if (keep>ndraw) {
    pandterm("'keep' cannot be larger than 'ndraw'")
  }

  # Determine method based on depvar
	dv <- unique(y)
	QRc <- FALSE; QRc.AL <- FALSE; QRb <- FALSE; QRb.AL <- FALSE;
  if (length(dv)>2){
	  if (alasso) {
		  QRc.AL <- TRUE
		} else {
		  QRc <- TRUE
		}
	} else if (length(dv)==2){
	  if (!identical(as.numeric(sort(dv)),c(0,1))){
		  pandterm("The dependent variable should be coded as vector of zeros and ones")
		}
	  if (alasso) {
		  QRb.AL <- TRUE
		} else {
		  QRb <- TRUE
		}
	} else {
	  pandterm("Something is wrong with the dependent variable")
	}

	# QRc 
	#================================
	if (QRc) {
	  # If missing, set default prior 
    if (is.null(prior)) {
		  prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
    } else if (class(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  } else if (prior$method!="QRc") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  }

    # Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    ndraw <- as.integer(ndraw)
    keep <- as.integer(keep)
    y <- as.double(y)
    quantile <- as.double(quantile)
    X <- as.double(X)
    beta0 <- as.double(prior$beta0)
    V0i <- as.double(chol2inv(chol(prior$V0)))
		shape0 <- as.double(prior$shape0)
		scale0 <- as.double(prior$scale0)
    betadraw <- double(nvar*ndraw/keep)
    sigmadraw <- double(ndraw/keep)

    # Call Fortran routine
    fn_val <- .Fortran("QRc_mcmc", n, nvar, ndraw, keep, y, quantile, X, beta0, V0i, shape0, scale0, betadraw, sigmadraw)

    # Return bayesQR object
		out <- list(method="QRc",
		            quantile=quantile,
								names=names,
                betadraw=matrix(fn_val[[12]], nrow=ndraw/keep,ncol=nvar),
                sigmadraw=fn_val[[13]]
					  		)

	# QRc.AL
	#================================
	} else if (QRc.AL) {
	  # If missing, set default prior 
    if (is.null(prior)) {
		  prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
    } else if (class(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  } else if (prior$method!="QRc.AL") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  }

    # Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    ndraw <- as.integer(ndraw)
    keep <- as.integer(keep)
    y <- as.double(y)
    quantile <- as.double(quantile)
    x <- as.double(X)
    a <- as.double(prior$a)
    b <- as.double(prior$b)
    c <- as.double(prior$c)
    d <- as.double(prior$d)
    betadraw <- double(nvar*ndraw/keep)
    sigmadraw <- double(ndraw/keep)

    # Call Fortran routine
    fn_val <- .Fortran("QRc_AL_mcmc",n, nvar, ndraw, keep, y, quantile, x, a, b, c, d, betadraw, sigmadraw) 

		# Return bayesQR object
		out <- list(method="QRc.AL",
	              quantile=quantile,
								names=names,
	              betadraw=matrix(fn_val[[12]], nrow=ndraw/keep, ncol=nvar),
	              sigmadraw=fn_val[[13]]
	              )
	
	# QRb
	#================================
	} else if (QRb) {
	  # If missing, set default prior 
    if (is.null(prior)) {
		  prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
    } else if (class(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  } else if (prior$method!="QRb") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  }

    # Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    ndraw <- as.integer(ndraw)
    keep <- as.integer(keep)
    y <- as.integer(y)
    quantile <- as.double(quantile)
    X <- as.double(X)
    beta0 <- as.double(prior$beta0)
    V0i <- as.double(chol2inv(chol(prior$V0)))
    betadraw <- double(nvar*ndraw/keep)
    sigmadraw <- double(ndraw/keep)

    # Call Fortran routine
    fn_val <- .Fortran("QRb_mcmc", n, nvar, ndraw, keep, y, quantile, X, beta0, V0i, betadraw)

    # Return bayesQR object
		out <- list(method="QRb",
		            quantile=quantile, 
								names=names,
                betadraw=matrix(fn_val[[10]], nrow=ndraw/keep,ncol=nvar)
					  		)

	# QRb.AL
	#================================
	} else if (QRb.AL) {
	  # If missing, set default prior 
    if (is.null(prior)) {
		  prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
    } else if (class(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  } else if (prior$method!="QRb.AL") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
	  }
	
    # Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    ndraw <- as.integer(ndraw)
    keep <- as.integer(keep)
    y <- as.integer(y)
    quantile <- as.double(quantile)
    x <- as.double(X)
    c <- as.double(prior$c)
    d <- as.double(prior$d)
    betadraw <- double(nvar*ndraw/keep)

    # Call Fortran routine
    fn_val <- .Fortran("QRb_AL_mcmc",n, nvar, ndraw, keep, y, quantile, x, c, d, betadraw) 

		# Return bayesQR object
		out <- list(method="QRb.AL",
	              quantile=quantile,
								names=names,
	              betadraw=matrix(fn_val[[10]], nrow=ndraw/keep, ncol=nvar)
	              )
	
	} else {
	  pandterm("Something is wrong with the dependent variable")
	}
	return(out)
} 
