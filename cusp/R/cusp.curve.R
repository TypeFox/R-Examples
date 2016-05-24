`cusp.curve` <-
function (alpha, beta) 
{
    roots <- polyroot(c(beta + alpha^2, 2 * alpha * beta, beta^2 - 
        3, -2 * alpha, -2 * beta, 0, 1))
    real <- abs(Im(roots)) < .Machine$double.eps^0.5
    if (all(real)) 
        sort(Re(roots))
    else rep(Re(roots[real]), 1)
}

