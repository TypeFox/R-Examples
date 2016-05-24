dmvnorm <- 
function (x, mean, sigma, log = FALSE, trustme = FALSE) 
{		
## Multivariate normal density taken from package mvtnorm
		if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (!is.null(dim(mean))) 
        dim(mean) <- NULL
 
  	dec <- chol(sigma)
   
    tmp <- forwardsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * length(mean) * 
        log(2 * pi) - 0.5 * rss
    names(logretval) <- rownames(x)
    if (log) 
        return(logretval)
    exp(logretval)
}
