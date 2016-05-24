MVTMLE0 <- function(X, nu = 1, prewhitened=FALSE, eps = 1e-06, maxiter = 100)
    {
    .Call( "MVTMLE0", X, nu, prewhitened, eps, maxiter, PACKAGE = "fastM")
    }
