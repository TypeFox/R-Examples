get_empirical_variances <-
function (sigma2, betahat, bin = 10, rescaleconst = NULL) 
{
    wipeout = function(betahat2, bin) {
        n = length(betahat2)
        indicator = rep(0, n)
        return(.C("wipeout", as.double(betahat2), as.integer(indicator), 
            as.integer(n), as.integer(bin))[[2]])
    }
    increasingreg = function(y) {
        n = length(y)
        return(.C("increasingreg", as.double(y), as.integer(n))[[1]])
    }
    getrescaleconst = function(bin, iterN = 50000) {
        if (bin == 10) 
            rescaleconst = 0.689
        else {
            for (i in 1:50000) stats[i] = max(rnorm(10)^2)
            rescaleconst = (bin - mean(stats))/(bin - 1)
        }
        return(rescaleconst)
    }
    if (is.null(rescaleconst)) 
        rescaleconst = getrescaleconst(bin)
    n = length(sigma2)
    neworder = order(sigma2)
    betahat2 = betahat[neworder]^2/rescaleconst
    adjustctl = as.logical(wipeout(betahat2, bin))
    yc = betahat2[adjustctl]
    ychat = increasingreg(yc)
    yhat = rep(-1, n)
    yhat[adjustctl] = ychat
    holes = which(!adjustctl)
    n1 = sum(!adjustctl)
    if (holes[1] == 1) 
        yhat[1] = yhat[2]
    else yhat[holes[1]] = yhat[holes[1] - 1]
    for (i in 2:n1) yhat[holes[i]] = yhat[holes[i] - 1]
    varbetahat = 1:n
    varbetahat[neworder] = yhat
    adjustctl[neworder] = adjustctl
    return(varbetahat)
}
