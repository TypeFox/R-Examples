`scaleChi` <-
function(Y, Kn, R0, eps = 0.000001)
  {
    R <- as.matrix(rowSums(Y))
    R0 <- as.matrix(R0)
    Keps <- Kn * (Kn > eps) + eps * (Kn < eps)
    Yr <- diag(1 / R) %*% Y %*% diag(1 / Keps) - 1
    Yr <- diag(sqrt(R0)) %*% Yr %*% diag(sqrt(Kn))
    Yr
  }

`scaleLin` <-
function(X, mean, sd, r, eps = 0.00000001)
  {
    ## [Xs] = scale_lin(X, mean, sd, r)
    ## centers and standardizes wrt to a prefined mean and
    ## standard deviation used in crossval_chi_lin
    X <- sweep(X, 2, mean)
    X <- sweep(X, 2, sd + eps, "/")
    r <- as.matrix(r)
    Xs <- diag(sqrt(r)) %*% X
    Xs
}
