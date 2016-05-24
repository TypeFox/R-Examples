covtopar <-
function(covg, covt){
#  covtopar() - var/cov component to fraction and correlation ( eg hsq and rg)
#             - same effect cross trait covariance cases
    n <- ncol(covg)
    if (ncol(covt) != n) {
        stop("Covg and Covt must be same size in covtopar():\n")
    }
    fract <- rep(0,n)
    corre <- matrix(0, n, n)
    for (j in 1:n) {
        for (i in 1:n) {
            if (i == j) {
                fract[i] <- covg[i, i]/covt[i, i]
                if (fract[i] < 0) {
                  fract[i] <- 0
                }
                corre[i,i] <- 1
            }
            else if (i != j) {
                able <- covg[i, i] * covg[j, j]
                if (able > 0) {
                  corre[i, j] <- covg[i, j]/sqrt(able)
                }
                else {
                  corre[i, j] <- 0
                }
            }
        }
    }
    return(list(corre=corre, fract=fract))
}
