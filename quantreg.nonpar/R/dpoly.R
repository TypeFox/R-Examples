dpoly <-
function (x, ..., degree = 1, coefs = NULL, raw = FALSE) 
{
  dots <- list(...)
  if (nd <- length(dots)) {
    if (nd == 1 && length(dots[[1]]) == 1) 
      degree <- dots[[1]]
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
    if (degree >= length(unique(x))) 
      stop("'degree' must be less than number of unique points")
    Z <- outer(x, 1:degree, "^")
    colnames(Z) <- 1:degree
    attr(Z, "degree") <- 1:degree
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
    if (QR$rank < degree) 
      stop("'degree' must be less than number of unique points")
    z <- QR$qr
    z <- z * (row(z) == col(z))
    raw <- qr.qy(QR, z)
    norm2 <- colSums(raw^2)
    alpha <- (colSums(x * raw^2)/norm2 + xbar)[1:degree]
    Z1 <- raw/rep(sqrt(norm2), each = length(x))
    norm2 <- c(1, norm2)
    x <- x + xbar
    if (degree > 1) 
      Z <- matrix(, length(x), n)
    Z[, 1] <- 0
    Z[, 2] <- 1
    for (i in 2:degree) Z[, i + 1] <- Z1[, i] + (x - alpha[i]) * 
      Z[, i] - (norm2[i + 1]/norm2[i]) * Z[, i - 1]
    Z <- Z/rep(sqrt(norm2[-1]), each = length(x))
    colnames(Z) <- 1:n - 1
    Z <- Z[, -1, drop = FALSE]
    attr(Z, "degree") <- 1:degree
    attr(Z, "coefs") <- list(alpha = alpha, norm2 = norm2)
    class(Z) <- c("poly", "matrix")
  }
  else {
    alpha <- coefs$alpha
    norm2 <- coefs$norm2
    Z1 <- matrix(, length(x), n)
    Z1[, 1] <- 1
    Z1[, 2] <- x - alpha[1]
    if (degree > 1) 
      for (i in 2:degree) Z1[, i + 1] <- (x - alpha[i]) * 
      Z1[, i] - (norm2[i + 1]/norm2[i]) * Z1[, i - 1]
    Z <- matrix(, length(x), n)
    Z[, 1] <- 0
    Z[, 2] <- 1
    if (degree > 1) 
      for (i in 2:degree) Z[, i + 1] <- Z1[, i] + (x - alpha[i]) * 
      Z[, i] - (norm2[i + 1]/norm2[i]) * Z[, i - 1]
    Z <- Z/rep(sqrt(norm2[-1]), each = length(x))
    colnames(Z) <- 0:degree
    Z <- Z[, -1, drop = FALSE]
    attr(Z, "degree") <- 1:degree
    attr(Z, "coefs") <- list(alpha = alpha, norm2 = norm2)
    class(Z) <- c("poly", "matrix")
  }
  Z
}
