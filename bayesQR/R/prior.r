prior <- function(formula=NULL, data=NULL, alasso=FALSE, ...){

  # Create function for error message
  pandterm <- function(message) {
    stop(message, call. = FALSE)
  }

  # Extract input
	if (is.null(formula)){
    pandterm("Formula has to be specified")
	}
	mf <- model.frame(formula=formula, data=data)
	nvar <- ncol(model.matrix(attr(mf,"terms"), data=mf))
	prior <- list(...)

  # Determine method based on depvar
	dv <- unique(model.response(mf))
	QRc <- FALSE; QRc.AL <- FALSE; QRb <- FALSE; QRb.AL <- FALSE;
  if (length(dv)>2){
	  if (alasso) {
		  QRc.AL <- TRUE
		} else {
		  QRc <- TRUE
		}
	} else if (length(dv)==2){
	  if (alasso) {
		  QRb.AL <- TRUE
		} else {
		  QRb <- TRUE
		}
	} else {
    pandterm("Dependent variable is constant")	
	}

  # Create correct prior object for given method
	if (QRc) {
	#================================
	  # If missing, set default prior 
    if (length(prior)==0) {
      prior$beta0 <- rep(0, nvar)
      prior$V0 <- 100 * diag(nvar)
			prior$shape0 <- 0.01
			prior$scale0 <- 0.01
		# Else check provided prior
    } else {
		  if(!all(names(prior) %in% c("beta0","V0","shape0","scale0"))){
		    pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
			} else {
        if (is.null(prior$beta0)) {
          prior$beta0 = rep(0, nvar)
        } else {
          prior$beta0 = prior$beta0
        }
        if (is.null(prior$V0)) {
          prior$V0 = 100 * diag(nvar)
        } else {
          prior$V0 = prior$V0
        }
        if (is.null(prior$shape0)) {
          prior$shape0 = 0.01 
        } else {
          prior$shape0 = prior$shape0
        }
        if (is.null(prior$scale0)) {
    	    prior$scale0 = 0.01
        } else {
          prior$scale0 = prior$scale0
        }
        if (ncol(prior$V0) != nrow(prior$V0) || ncol(prior$V0) != nvar || nrow(prior$V0) != nvar) {
          pandterm("Bad dimensions for V0")
        }
        if (length(prior$beta0) != nvar) {
          pandterm("Wrong lenghth for beta0")
        }
  		}
  	}
		prior$method <- "QRc"

	# QRc.AL
	#================================
	} else if (QRc.AL) {
	  # If missing, set default prior 
    if (length(prior)==0) {
      prior$a = 0.01
      prior$b = 0.01
      prior$c = 0.01
      prior$d = 0.01
		# Else check provided prior
    } else {
		  if(!all(names(prior) %in% c("sigma_shape","sigma_scale","etasq_shape","etasq_scale"))){
		    pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
      } else {
        if (is.null(prior$sigma_shape)) {
  	      prior$a = 0.01
        } else {
          prior$a = prior$sigma_shape
        }
        if (is.null(prior$sigma_scale)) {
  	      prior$b = 0.01
        } else {
          prior$b = prior$sigma_scale
        }
        if (is.null(prior$etasq_shape)) {
  	      prior$c = 0.01
        } else {
          prior$c = prior$etasq_shape
        }
        if (is.null(prior$etasq_scale)) {
  	      prior$d = 1/0.01
        } else {
          prior$d = 1/prior$etasq_scale
        }
      }
		}
		prior$method <- "QRc.AL"

	# QRb
	#================================
	} else if (QRb) {
	  # If missing, set default prior 
    if (length(prior)==0) {
      prior$beta0 <- rep(0, nvar)
      prior$V0 <- 100 * diag(nvar)
		# Else check provided prior
    } else {
		  if(!all(names(prior) %in% c("beta0","V0"))){
		    pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
			} else {
        if (is.null(prior$beta0)) {
          prior$beta0 = rep(0, nvar)
        } else {
          prior$beta0 = prior$beta0
        }
        if (is.null(prior$V0)) {
          prior$V0 = 100 * diag(nvar)
        } else {
          prior$V0 = prior$V0
        }
        if (ncol(prior$V0) != nrow(prior$V0) || ncol(prior$V0) != nvar || nrow(prior$V0) != nvar) {
          pandterm("Bad dimensions for V0")
        }
        if (length(prior$beta0) != nvar) {
          pandterm("Wrong lenghth for beta0")
        }
			}
		}
		prior$method <- "QRb"

	# QRb.AL
	#================================
	} else if (QRb.AL) {
	  # If missing, set default prior 
    if (is.null(prior)) {
      c = 0.01
      d = 0.01
		# Else check provided prior
    } else {
		  if(!all(names(prior) %in% c("lambdasq_shape","lambdasq_scale"))){
		    pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
      } else {
        if (is.null(prior$lambdasq_shape)) {
  	      prior$c = 0.01
        } else {
          prior$c = prior$lambdasq_shape
        }
        if (is.null(prior$lambdasq_scale)) {
  	      prior$d = 1/0.01
        } else {
          prior$d = 1/prior$lambdasq_scale
        }
      }
		}
	  prior$method <- "QRb.AL"
	}

	# Define class of object and return
	class(prior) <- "bayesQR.prior"
	return(prior)
}
