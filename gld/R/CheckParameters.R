.gl.parameter.tidy <- function(lambda1,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL) 
{
# parameter labels in .gl.parameter.tidy are not correct for GPD parameterisation, but the rest should work
# Don't allow characters in lambda5 - common error with parameterisation stuff
if(is.character(lambda5)) {stop(paste("lambda5=",lambda5,"It should be a number between -1 and 1"))}
# Don't allow numbers in parameterisation - included as a warning here, so the main one is a stop.
if(!is.character(param)) {warning(paste("param=",param,"It shouldn't be a number, it should be a string describing the parameterisation"))}
if(is.null(lambda1)) { stop("No values provided for lambda parameters, argument lambda1 is NULL") }
if(length(lambda1) > 1) #using a vector for the parameters.  
	# Check that there aren't values in the individual lambda arguments
	{
	if (!(is.null(lambda2) & is.null(lambda3)& is.null(lambda4) & is.null(lambda5)) ) 
		{ stop("Call includes vector version of the lambda parameters as well as the \nscalar version") }
	if ((length(lambda1) < 4) | (length(lambda1) > 5 ) )  
		{ stop(paste("argument lambda1 has length", length(lambda1),"\nThis should be 1 (lambda parameters given as seperate arguments), 4 (vector argument \n for RS, FKML or GPD parameterisation) or 5 (vector argument for fm5 parameterisation")) }
	if (length(lambda1)== 5)
		{ if (param != "fm5") { 
			stop(paste("argument lambda1 has length",length(lambda1),"which is not valid for the",param,"\nparameterisation")) 
			}
		# else --- fm5, in vector form, ready for gl.check.lambda 
		}
	if (length(lambda1)== 4)
		{ if (param == "fm5" ) 
			{ stop(paste("argument lambda1 has length 4, which is not valid for the fm5 \nparameterisation")) }
		# else --- 4 parameter versions in vector form, ready for gl.check.lambda 
		}
	}
else { # single parameter arguments - check they are there, then collect them together.  Ideally I would have a special section here to deal with GPD
	if (is.null(lambda2)) { stop("No value for lambda2") }
	if (is.null(lambda3)) { stop("No value for lambda3") }
	if (is.null(lambda4)) { stop("No value for lambda4") }
	if ((is.null(lambda5)) & param=="fm5" ) { stop("No value for lambda5") }
	if (!(is.null(lambda5)) & param!="fm5") { stop(paste("lambda5=",lambda5," but there is no lambda 5 for the\n",param,"parameterisation")) }
	if (param != "fm5") { # A 4 parameter version
		lambda1 <- c(lambda1,lambda2,lambda3,lambda4)
		}
	else { # fm5
		lambda1 <- c(lambda1,lambda2,lambda3,lambda4,lambda5)
		}
	}
# There is now an error if there is the wrong number of parameters, and 
# lambda1 returned as a vector with 4 or 5 elements
# as.double is needed to remove data.frame attributes if lambda1 was
# extracted from a data.frame
as.double(lambda1)
}

gl.check.lambda <- function (lambdas, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, param = "fkml", lambda5 = NULL, vect = FALSE)
{
    if (vect) {
        if (!is.null(lambda3)) {
            warning("lambda3 should be null because you claim the parameters are in a vector")
        }
    }
    else {
        lambdas <- .gl.parameter.tidy(lambdas, lambda2, lambda3, 
            lambda4, param, lambda5)

    }
    if (param == "fm5") {
        lambda5 = lambdas[5]
    }
    lambda4 = lambdas[4]
    lambda3 = lambdas[3]
    lambda2 = lambdas[2]
    lambda1 = lambdas[1]
    param <- switch(param, freimer = , frm = , FMKL = , FKML = , 
        fmkl = , fkml = {
            if (lambda2 <= 0) {
                return(FALSE)
            } else {
                return(TRUE)
            }
        }, ramberg = , ram = , RS = , rs = {
            if (lambda3 * lambda4 > 0) {
                if ((lambda3 > 0) & (lambda4 > 0)) {
                  if (lambda2 <= 0) {
                    ret <- FALSE
                  } else {
                    ret <- TRUE
                  }
                }
                if ((lambda3 < 0) & (lambda4 < 0)) {
                  if (lambda2 >= 0) {
                    ret <- FALSE
                  } else {
                    ret <- TRUE
                  }
                }
            } else {
                if (lambda2 >= 0) {
                  return(FALSE)
                }
                if ((lambda3 > 0) & (lambda3 < 1) & (lambda4 < 
                  0)) {
                  return(FALSE)
                }
                if ((lambda4 > 0) & (lambda4 < 1) & (lambda3 < 
                  0)) {
                  return(FALSE)
                }
                lc <- lambda3
                ld <- lambda4
                if ((lambda3 > -1) & (lambda3 < 0) & (lambda4 > 
                  1)) {
                  if (((1 - lc)^(1 - lc) * (ld - 1)^(ld - 1))/((ld - 
                    lc)^(ld - lc)) > -lc/ld) {
                    return(FALSE)
                  } else {
                    return(TRUE)
                  }
                }
                if ((lambda4 > -1) & (lambda4 < 0) & (lambda3 > 
                  1)) {
                  if (((1 - ld)^(1 - ld) * (lc - 1)^(lc - 1))/((lc - 
                    ld)^(lc - ld)) > -ld/lc) {
                    return(FALSE)
                  } else {
                    return(TRUE)
                  }
                }
                if (lambda3 == 0) {
                  if (lambda4 > 0) {
                    if (lambda2 < 0) {
                      return(FALSE)
                    }
                    ret <- TRUE
                  }
                  if (lambda4 == 0) {
                    warning("RS parameterisation: lambda3 and lambda4 zero gives a point mass at lambda1")
                  }
                  if (lambda4 <0) {
                    if (lambda2 > 0) {
                      return(FALSE)
                    }
                  ret <- TRUE
                  }
                }
                if (lambda4 == 0) {
                  if (lambda3 > 0) {
                    if (lambda2 < 0) {
                      return(FALSE)
                    }
                    ret <- TRUE
                  }
                  if (lambda3 == 0) {
                    warning("RS parameterisation: lambda3 and lambda4 zero gives a point mass at lambda1")
                  }
                  if (lambda3 <0) {
                    if (lambda2 > 0) {
                      return(FALSE)
                    }
                    ret <- TRUE
                  }
                }
                if (is.null(ret)) {warning("RS param return not set: please email maintainer with example")
                ret <- TRUE}
            }
        }, fm5 = {
            lambda5 <- lambdas[5]
            if (lambda2 <= 0) {
                ret <- FALSE
            } else {
                if ((lambda5 >= -1) | (lambda5 <= 1)) {
                  ret <- TRUE
                } else {
                  ret <- FALSE
                }
            }
        }, vsk=, VSK =, gpd=, GPD= {
		if (lambda2 <= 0) {
			ret <- FALSE
			warning("Negative or zero beta")
		} else {  ## delta check
			if ((lambda3 < 0)|(lambda3 > 1)) { ret <- FALSE 
			} else 	{
				ret <- TRUE 
			} 
		}	
	}, stop("Error when checking validity of parameters.\n Parameterisation must be fmkl, rs, gpd, vsk or fm5")) 
ret
}
