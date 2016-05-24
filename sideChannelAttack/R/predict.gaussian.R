p <-
function (x_, mean_, cov_, detCov_) 
{
    N = length(x_)
    if (detCov_ == 0) 
        stop("the determinant of the covariance is zero")
    return(exp(-0.5 * (x_ - mean_) %*% ginv(cov_) %*% matrix(c(x_ - 
        mean_), ncol = 1))/sqrt(abs(((2 * pi)^N) * detCov_)))
}
keyOfList <-
function (x_)
{
    keys_ = c()
    n = 1
    for(i in x_)
    {
        if(!is.null(i))
            keys_ = c(keys_,n)
        n=n+1
    }
    return (keys_)
}
predict.gaussian <-
function (object, newdata, ...) 
{
	if(is.matrix(newdata))
		nbreDExempleMax=dim(newdata)[1]
	else
		nbreDExempleMax=1
	res=c()
	for(nbreDExemple in 1:nbreDExempleMax)
	{
		maximum = -1
		maximumIndice = -1
		for (supposonsLaclCle in keyOfList(object$mean)) {
			if(is.matrix(newdata))
			{
				proba <- p(newdata[nbreDExemple,], (object$mean)[[supposonsLaclCle]], 
					(object$cov)[[supposonsLaclCle]], (object$detCov)[[supposonsLaclCle]])
			}
			else
			{
				proba <- p(newdata, (object$mean)[[supposonsLaclCle]], 
						   (object$cov)[[supposonsLaclCle]], (object$detCov)[[supposonsLaclCle]])
			}
			if (maximum < proba) {
				maximumIndice = supposonsLaclCle
				maximum = proba
			}
		}
		res=c(res,maximumIndice - 1)
	}
    return (factor(res))
}
