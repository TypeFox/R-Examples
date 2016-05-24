qqgl <- function (y=NULL,lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fkml",lambda5=NULL,abline=TRUE,lambda.pars1=NULL,lambda.pars2=NULL,
  param2="fkml",points.for.2.param.sets=4000,...)
{
# standard parameter fixin' - copied directly from dgl, but we want the 
# warnings to happen in this function.
# Tidy the parameters so gl.check.lambda will work
if(is.null(lambda2)&&length(lambda1)==1){ #lambda2 is NULL and length(lambda1)=1 (see
  	# default above, so parameters must be set using lambda.pars1 
	# and possibly lambda.pars2
  	lambdas <- .gl.parameter.tidy(lambda.pars1)
  	} else{ #lambda2 is not NULL, or lambda1 is a vector longer than 1, so parameters are defined old-school - tidy them the old way
	# If you are doing it this way, you shouldn't use lambda.pars1
	if(!is.null(lambda.pars1)) {stop("Don\'t use variables lambda1, lambda2 etc. in the same function call as lambda.pars1")}
	lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
	}
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
if(is.null(lambda.pars2)) { # No second set of parameters - compare to a dataset
	n <- length(y)
	}
else { # A second set of parameters - compare two distributions
	if(!is.null(y)){stop("Can\'t produce a QQ plot comparing 2 sets of parameter values and also a dataset")}
	lambdas2 <- .gl.parameter.tidy(lambda.pars2,param=param2)
	if(!gl.check.lambda(lambdas2,param=param2,vect=TRUE)) {
        	stop(paste("The second set of parameter values", lambdas2,
		"\ndo not produce a proper distribution with the",param2,
		"parameterisation - see \ndocumentation for gl.check.lambda"))
	}
	n <- points.for.2.param.sets
	}
u <- seq(from = 1/(n + 1), by = 1/(n + 1), length = n)
    # change this to use ideal depths
if(is.null(lambda.pars2)) { # No second set of parameters - compare to a dataset
	Theoretical.Quantiles <- qgl(u, lambda1=lambdas, param=param)
	Data <- y
	ret <- qqplot(Theoretical.Quantiles,Data,...)
	}
else	{ # A second set of parameters - compare two distributions
	Quantiles1 <- qgl(u, lambda1=lambdas, param=param)
	Quantiles2 <- qgl(u, lambda1=lambdas2, param=param2)
	ret <- qqplot(Quantiles1,Quantiles2,...)
	}
if(abline) { 
	abline(0,1)
	}
invisible(ret)
}
