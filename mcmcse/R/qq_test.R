qqTest <- function(x, covmat, g = NULL)
{
	chain <- as.matrix(x)
	if(!is.matrix(x) && !is.data.frame(x))
	  stop("'x' must be a matrix or data frame.")

	if (is.function(g)) 
	{
	  chain <- apply(x, 1, g)

	  if(is.vector(chain))
	  {
	    chain <- as.matrix(chain)
	  }else
	  {
	    chain <- t(chain)
	  }
	}

	mu <- colMeans(chain)
	n <- dim(chain)[1]
	p <- dim(chain)[2]
	decomp  <- svd(covmat)
	inv.root <- decomp$v %*% diag( (decomp$d^(-1/2)), p) %*% t(decomp$u)
	qqnorm(inv.root%*%mu)
	qqline(inv.root%*%mu)
}
