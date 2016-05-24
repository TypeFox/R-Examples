ldb <- 
    function(means,
        covs,
	covstruct = c("unstructured", "scaledIdentity", "diagonal", "identity"), 
	noise = 10)
{
    covstruct <- match.arg(covstruct)
    if(!is.list(means)) stop("means is not a list")
    else M <- means
    if(length(M)==1) M[[2L]] <- M[[1L]]
    if(is.list(covs)) covs <- as.matrix(covs[[1L]])
    if(covstruct == "scaledIdentity")
	covs <- diag(mean(diag(covs)), nrow=dim(covs)[1L], ncol=dim(covs)[2L])
    else if(covstruct == "diagonal")
	covs <- diag(diag(covs), nrow=dim(covs)[1L], ncol=dim(covs)[2L])
    else if(covstruct == "identity")
	covs <- diag(1, nrow=dim(covs)[1L], ncol=dim(covs)[2L])
    
    if(qr(covs)$rank != max(dim(covs)))
	stop("Covariance matrix not full rank.")
    offset <- ifelse( M[[2L]] == M[[1L]], .Machine$double.eps, 0)
    b <- ((M[[2L]] - M[[1L]]) + offset) %*% qr.solve(covs)
    enorm <- sqrt(sum(b^2))
    res <- NULL
    res$noise <- noise
    res$coeffs <- as.vector(b/enorm)
    res$bias <- as.vector( (-.5 * b %*% (M[[2L]] + M[[1L]])) /enorm)
    class(res) <- c("glcStruct", "list")
    res
}

ldb.p.correct <- function(means, covs, noise = 0)
{
    if(!is.list(means)) stop("means is not a list")
    else M <- means
    if(length(M)==1) M[[2L]] <- M[[1L]]
    if(is.list(covs) && length(covs) > 1)
        K <- as.matrix(Reduce("+", covs)/length(covs))
    else K <- as.matrix(covs)
    g <- length(M[[1L]])
    if(g < 2) return(.5)
    else if(g > 2)
        stop("A number of categories greater than 2 is not currently supported.")
    
    Kinv = qr.solve(K)
    b <- t(M[[2L]] - M[[1L]]) %*% Kinv
    muhx = b %*% M[[1L]] + 
        .5*t(M[[1L]]) %*% Kinv %*% M[[1L]] -
        .5*t(M[[2L]]) %*% Kinv %*% M[[2L]]
    varhx = b %*% (M[[2L]] - M[[1L]])
    # normalize
    enorm <- sqrt(sum(b^2))
    muhx <- muhx / enorm
    varhx <- varhx / (enorm^2)
    #obtain the z score from mu and var
    z <- -muhx/sqrt(varhx + noise^2)
    pc <- as.vector(pnorm(z, mean=0, sd=1))
    pc
}
