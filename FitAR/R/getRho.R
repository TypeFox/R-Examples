`getRho` <-
function(ans){
#rho test-statistic
    phiHat <- ans$phiHat
    res <- ans$res
    n <- length(res)
    n*(phiHat-1)
    }

