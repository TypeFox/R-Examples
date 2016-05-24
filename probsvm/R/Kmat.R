Kmat <- function(x, y, kernel, kparam = NULL) {
if (kernel == "linear") {
obj <- x %*% t(y)
} else if (kernel == "polynomial") {
obj <- (x %*% t(y) + 1)^kparam
} else if (kernel == "radial") {
normx <- drop((x^2) %*% rep(1, ncol(x)))
normy <- drop((y^2) %*% rep(1, ncol(y)))
temp <- x %*% t(y)
temp <- (-2 * temp + normx) + outer(rep(1, nrow(x)), normy, "*")
obj <- exp(-temp * kparam)
} else obj <- NULL
obj
}
