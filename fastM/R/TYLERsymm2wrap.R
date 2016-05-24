TYLERsymm2 <- function(X, eps = 1e-06, maxiter = 100)
    {
    .Call( "Tylersymm2", X, eps, maxiter, PACKAGE = "fastM")
    }


