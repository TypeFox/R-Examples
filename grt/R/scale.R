scale <- function(x, ...) UseMethod("scale")

scale.glc  <- function(x, initdb = FALSE, zlimit = Inf, ...)
{
    if (!inherits(x, "glc")) stop("object not of class \"glc\"")
    data <- x$model
    Terms <- x$terms
    X <- model.matrix(delete.response(Terms), data)
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    if(initdb) params <- x$initpar
    else params <- x$par
    #calculate z-scores
    z_coefs <- c(params$coeffs, params$bias) / params$noise
    z <- drop( cbind(X,1) %*% as.matrix(z_coefs) )
    z[z >  zlimit] <- zlimit
    z[z < -zlimit] <- -zlimit
    z
}

scale.gqc <- function(x, initdb = FALSE, zlimit = Inf, ...)
{
    if (!inherits(x, "gqc")) stop("object not of class \"gqc\"")
    data <- x$model
    Terms <- x$terms
    X <- model.matrix(delete.response(Terms), data)
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    if(initdb) params <- x$initpar
    else params <- x$par
    #calculate z-scores
    dimen <- ncol(X)
    coeffs <- params$coeffs
    ncoeff <- length(coeffs)
    A <- diag(coeffs[1:dimen], nrow=dimen, ncol=dimen)
    A[lower.tri(A)] <- A[upper.tri(A)] <- (.5 * coeffs[(dimen+1):(ncoeff-dimen)])
    b <- coeffs[(ncoeff-dimen+1):(ncoeff)]
    c0 <- as.vector(params$bias)
    pnoise <- params$pnoise
    cnoise <- params$cnoise
    comb <- cbind(rbind(1:dimen,1:dimen),combn(1:dimen,2))
    tmp1 <- pnoise*sum(diag(A))
    meanhxs <- tmp1 + cbind(X[,comb[1,]] * X[,comb[2,]], X, 1) %*% c(coeffs, c0)
    bvec <- matrix(rep(b,nrow(X)), ncol=nrow(X), nrow=dimen, byrow=FALSE)
    varhxs <- 2*tmp1*tmp1 + pnoise*colSums((bvec + 2*A %*% t(X))^2)
    z <- drop( meanhxs / sqrt(varhxs+cnoise) )
    z[z >  zlimit] <- zlimit
    z[z < -zlimit] <- -zlimit
    z
}