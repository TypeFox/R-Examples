rgl <- function(n,lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,param="fkml",lambda5=NULL)
{
# Check the parameters
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the values are OK
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,"\ndo not produce a proper distribution for the",param,"parameterisation \n - see documentation for gl.check.lambda"))
	}
# Produce the uniform data
p <- runif(n)
# convert to gl
res <- qgl(p,lambda1=lambdas,param=param)
res
}
