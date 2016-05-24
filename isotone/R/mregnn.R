## z = Xb unrestricted
mregnn <- function (x, y, a) {
  k <- qr.Q(qr(x))
  u <- drop (crossprod(k, y))
  v <- -crossprod (k, t(a))
  lb <- nnls(v, u)$x
  xb <- drop(k %*% (u - v %*% lb))
  return (list(xb = xb, lb = lb, f = sum((y - xb) ^ 2)))
}

## z monotone restricted
mregnnM <- function (x, y) {
  k <- qr.Q(qr(x))
  u <- drop (crossprod(k, y))
  v <- -t(diff(k))
  lb <- nnls(v, u)$x
  xb <- drop(k %*% (u - v%*% lb))
  return (list(xb = xb, lb = lb, f = sum((y - xb) ^ 2)))
}

## z positive
mregnnP <- function (x, y) {
  k <- qr.Q(qr(x))
  u <- drop (crossprod(k, y))
  v <- -t(k)
  lb <- nnls(v, u)$x
  xb <- drop(k %*% (u - v%*% lb))
  return (list(xb = xb,lb = lb, f = sum((y - xb) ^ 2)))
}
