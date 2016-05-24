qdb <- 
    function(means,
        covs,
	pnoise = 10, 
	cnoise = 100,
	sphere = FALSE)
{
    if(!is.list(means)) stop("means is not a list")
    if(sphere){
        dimen <- length(means[[1L]])
	A <- diag(1, nrow=dimen, ncol=dimen)
	b <- -2*means[[1L]]
	c0 <- sum(means[[1L]]^2) - sum((means[[2L]] - means[[1L]])^2)/4
    } else {
        if(!is.list(covs)) stop ("covs is not a list")
        if(length(covs) < 2) stop("length(covs) is less than 2")
        Kinv <- lapply(covs, qr.solve)
        A <- -.5*( Kinv[[2L]] - Kinv[[1L]])
        b <- t(means[[2L]] %*% Kinv[[2L]] - means[[1L]] %*% Kinv[[1L]])
        c0 <- .5*(-log(det(covs[[2L]])/det(covs[[1L]])) - 
            means[[2L]] %*% Kinv[[2L]] %*% means[[2L]] + 
            means[[1L]] %*% Kinv[[1L]] %*% means[[1L]])
    }
    res <- NULL
    res$pnoise <- pnoise
    res$cnoise <- cnoise
    res$coeffs <- c(diag(A), 2*A[upper.tri(A)], b)
    res$bias <- as.vector(c0)
    class(res) <- c("gqcStruct", "list")
    res
}

qdb.p.correct <- function(x, qdb, refpts = colMeans(x))
{
    x <- as.matrix(x)
    dimen <- ncol(x)
    N <- nrow(x)
    if(length(refpts) != dimen) 
        stop("length(refpts) and ncol(x) are different")
    if(inherits(qdb, "gqcStruct")) 
        qdb <- unlist(qdb[c("coeffs","bias")])
    if(dimen < 2 | dimen > 3) stop("Unsupported number of dimenensions.")
    if((sum(1:dimen) + dimen) != length(qdb)-1) 
        stop("Inappropriate number of coeffs")
    idx <- 1:dimen
    comb <- cbind(rbind(idx,idx),combn(idx,2))
    refs <- -sign(c(refpts[comb[1,]] * refpts[comb[2,]], refpts, 1) %*% qdb)
    hvec <- cbind(x[,comb[1,]] * x[,comb[2,]], x, 1) %*% qdb
    hvec <- as.vector(refs) * as.vector(hvec)
    pc <- length(hvec[hvec < 0])/N
    pc
}