`cusp.bifset` <-
function (beta) 
{
    alpha = 2 * (beta/3)^(3/2)
    cbind(beta = beta, alpha.l = -alpha, alpha.u = alpha)
}

