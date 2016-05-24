partocov <-
function(par, covt){
#  partocov()
    n <- ncol(par)
    if (ncol(covt) != n) {
        stop("Par and covt must be same size in partocov():\n")
    }
    covg <- matrix(0, n, n)
    for (i in 1:n) {
        covg[i, i] <- par[i, i] * covt[i, i]
    }
    for (j in 1:n) {
        for (i in 1:n) {
            if (i != j) {
                able <- covg[i, i] * covg[j, j]
                if (able > 0) {
                  covg[i, j] <- par[i, j] * sqrt(covg[i, i] * 
                    covg[j, j])
                }
                else {
                  covg[i, j] <- 0
                }
            }
        }
    }
    return(covg)
}
