##### Estimation for Regularized Logit Multinomial Regression  ######

## currently nothing more than a call to dmr
mnlm <- function(cl, covars, counts, mu=NULL, bins=NULL, verb=0, ...)
{
  out <- dmr(cl, covars, counts, mu, bins, verb, ...)
  class(out) <- c("dmr", class(out))
  return(out)
}  

## sufficient reduction projection
srproj <- function(obj, counts, dir=1:K, ...)
{
	if(inherits(obj,"dmr")) obj <- coef(obj, ...)
	if(!inherits(obj,"dmrcoef"))
		stop("obj must be of class dmr or dmrcoef.")

	# row totals
	m <- rowSums(counts)

    # remove the intercept
    obj <- obj[-1,,drop=FALSE]
    # use specified directions
    K <- nrow(obj)
    obj <- obj[dir,,drop=FALSE]

    # project
    z <- tcrossprod(counts,obj)
    # normalize
    z <- as.matrix(z)/(m + 1*(m==0))
    colnames(z) <- rownames(obj)
    rownames(z) <- rownames(counts)
    cbind(z,m=m)
}

