`cusp.extrema` <-
function (alpha, beta) 
{
    roots <- polyroot(c(alpha, beta, 0, -1))
    real <- abs(Im(roots)) < .Machine$double.eps^0.5
    if (all(real)) 
        sort(Re(roots))
    else rep(Re(roots[real]), 3)
}

