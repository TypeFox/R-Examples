ipfp.control <-
function (tol = 1e-06, maxit = 200) 
{
    if (!is.numeric(tol) || tol <= 0) 
        stop("value of the IPFP 'tol' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of the IPFP iterations must be > 0")
    list(tol = tol, maxit = maxit)
}

