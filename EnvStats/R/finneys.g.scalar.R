finneys.g.scalar <-
function (m, z, n.terms.inc = 10, iter = 100, tol = .Machine$double.eps) 
{
    if (any(length.list(m, z, n.terms.inc, iter, tol) > 1)) 
        stop("'m', 'z', 'n.terms.inc', 'iter', and 'tol' must be scalars")
    if (m < 1 || m != trunc(m)) 
        stop("'m' must be an integer greater than or equal to 1")
    if (n.terms.inc < 3 || n.terms.inc != trunc(n.terms.inc)) 
        stop("'n.terms.inc' must be an integer greater than or equal to 3")
    if (iter < 1 || iter != trunc(iter)) 
        stop("'iter' must be an integer greater than or equal to 1")
    if (tol < .Machine$double.eps) 
        stop("'tol' must be at least as large as .Machine$double.eps")
    for (i in 1:iter) {
        vec <- finneys.g.scalar.sub(m, z, n.terms = n.terms.inc * 
            i)
        if (abs(vec[n.terms.inc * i]) <= tol) 
            return(sum(vec))
    }
    stop(paste("Failure to converge with n.terms.inc = ", n.terms.inc, 
        "and iter = ", iter))
}
