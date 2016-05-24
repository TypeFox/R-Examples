mypoly <-
  function (x, ..., degree = 1, coefs = NULL, raw = FALSE) 
{
  dots <- list(...)
  if (nd <- length(dots)) {
    if (nd == 1 && length(dots[[1L]]) == 1L) 
      degree <- dots[[1L]]
    else return(polym(x, ..., degree = degree, raw = raw))
  }
  if (is.matrix(x)) {
    m <- unclass(as.data.frame(cbind(x, ...)))
    return(do.call("polym", c(m, degree = degree, raw = raw)))
  }
  if (degree < 1) 
    stop("'degree' must be at least 1")
  if (any(is.na(x))) 
    stop("missing values are not allowed in 'poly'")
  n <- degree + 1
  if (raw) {
    Z <- outer(x, 1L:degree, "^")
    colnames(Z) <- 1L:degree
    attr(Z, "degree") <- 1L:degree
    class(Z) <- c("poly", "matrix")
    return(Z)
  }
  if (is.null(coefs)) {
    if (degree >= length(unique(x))) 
      stop("'degree' must be less than number of unique points")
    xbar <- mean(x)
    x <- x - xbar
    X <- outer(x, seq_len(n) - 1, "^")
    QR <- qr(X)
browser()
    if (QR$rank < degree) 
      stop("'degree' must be less than number of unique points")
    z <- QR$qr
    z <- z * (row(z) == col(z))
    raw <- qr.qy(QR, z)
    norm2 <- colSums(raw^2)
    alpha <- (colSums(x * raw^2)/norm2 + xbar)[1L:degree]
    Z <- raw/rep(sqrt(norm2), each = length(x))
    colnames(Z) <- 1L:n - 1L
    Z <- Z[, -1, drop = FALSE]
    attr(Z, "degree") <- 1L:degree
    attr(Z, "coefs") <- list(alpha = alpha, norm2 = c(1, 
                                              norm2))
    class(Z) <- c("poly", "matrix")
  }
  else {
    alpha <- coefs$alpha
    norm2 <- coefs$norm2
    Z <- matrix(, length(x), n)
    Z[, 1] <- 1
    Z[, 2] <- x - alpha[1L]
    if (degree > 1) 
      for (i in 2:degree) Z[, i + 1] <- (x - alpha[i]) * 
        Z[, i] - (norm2[i + 1]/norm2[i]) * Z[, i - 1]
    Z <- Z/rep(sqrt(norm2[-1L]), each = length(x))
    colnames(Z) <- 0:degree
    Z <- Z[, -1, drop = FALSE]
    attr(Z, "degree") <- 1L:degree
    attr(Z, "coefs") <- list(alpha = alpha, norm2 = norm2)
    class(Z) <- c("poly", "matrix")
  }
  Z
}
