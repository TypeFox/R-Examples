###########################################
## Functions to estimate the parameters  ##
###########################################

## Estimate F
fj_pl <- function(x, y, pij, sigma2, lambda,
                  B = getbasismatrix(x,
                                     create.bspline.basis(range(x),
                                                          min(N %/% 4, 40))),
                  R = getbasispenalty(create.bspline.basis(range(x),
                                                           min(N %/% 4, 40))), ...) {
    N <- length(y)
    w_j <- (pij/sigma2) / (sum(pij/sigma2) / N)
    hat_matrix <- B %*% solve(t(B) %*% diag(w_j) %*% B +
                              2*exp(lambda)*R) %*%
                                  t(B) %*% diag(w_j)
    list(f_hat = hat_matrix %*% y,
         H = hat_matrix,
         trace = sum(pij * diag(hat_matrix)))
}
