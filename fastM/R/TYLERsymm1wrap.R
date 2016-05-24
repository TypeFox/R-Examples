TYLERsymm1 <- function(X, eps = 1e-06, maxiter = 100)
    {
    .Call( "Tylersymm1", X, eps, maxiter, PACKAGE = "fastM")
    }


