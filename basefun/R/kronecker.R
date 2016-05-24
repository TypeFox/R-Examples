
### row tensor
.Rt <- function(X1, X2 = NULL) {

   if (is.null(X2)) X2 <- X1
   stopifnot(nrow(X1) == nrow(X2))
   p1 <- ncol(X1)
   p2 <- ncol(X2)
   ip1 <- 1:p1
   ip2 <- 1:p2
   ret <- X1[, rep(ip1, rep(p2, p1)), drop = FALSE] * 
          X2[, rep(ip2, p1), drop = FALSE]
   return(ret)
}

# H-transform of an array A by a matrix X
.H <- function(X, A, d, G = TRUE, crossp = G) {
    stopifnot(nrow(A) == d[1])
    if (G) {
        ### row tensor
        X <- .Rt(X)
    }
    if (crossp) {
        XA <- crossprod(X, A)
    } else {
        XA <- X %*% A
    }
    array(XA, c(nrow(XA), d[-1]))
}

# Rotation of an array A
.Rotate <- function(A) {
    d <- 1:length(dim(A))
    aperm(A, c(d[-1], d[1]))
}

# Rotated H-transform of an array A by a matrix X
.RH <- function(X, A, d = dim(A), G = TRUE, crossp = TRUE) {
    if (length(dim(A)) > 2)
        A <- matrix(A, nrow = d[1])
    .Rotate(.H(X, A, d = d, G = G, crossp = crossp))
}

### X = X3 %x% X2 %x% X1
.cXb <- function(..., beta) {
    X <- list(...)
    ret <- .RH(X[[1]], A = beta, G = FALSE, crossp = FALSE)  
    for (i in 2:length(X)) 
        ret <- .RH(X[[i]], A = ret, G = FALSE, crossp = FALSE)
    ret
}

