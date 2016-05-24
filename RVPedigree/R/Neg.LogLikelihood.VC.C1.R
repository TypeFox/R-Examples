Neg.LogLikelihood.VC.C1 <- function(h, eps, U, S) # f 4.2 articl
{
    B <- h * (S - 1) + 1
    return(sum(log(B)) + t(eps) %*% U %*% diag(1/B) %*% t(U) %*% eps)
}
