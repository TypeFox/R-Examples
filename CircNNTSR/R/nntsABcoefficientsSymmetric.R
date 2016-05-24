nntsABcoefficientsSymmetric<-
function (cpars = c(0, 0), M = 0) 
{
    if (M == 0) 
        return("Uniform case: No coefficients AB")
    size <- length(cpars)
    if ((size%%2) == 1) 
        return("Length of cpar must be an even number")
    if (size < 2 * M) {
        temp <- size/2 + 1
        cparscorr <- c(cpars[1:(size/2)], array(0, (M - temp + 
            1)), cpars[temp:size], array(0, (M - temp + 1)))
        cpars <- cparscorr
        cat("Warning: Missing parameters set to 0\n")
    }
    if (M > 0) {
        if (1/(2 * pi) - sum(cpars[1:M]) < 0) 
            return("sum of componentes greater than condition")
        else {
            aux <- complex(M + 1)
            aux[1] <- sqrt(1/(2 * pi) - sum(cpars[1:M]))
            for (k in 1:M) {
                aux[k + 1] <- sqrt(cpars[k]) * exp((0+1i) * cpars[(k + 
                  M)])
            }
            aux2 <- complex(M)
            for (j in 1:M) {
                for (k in 1:(M - j + 1)) {
                  aux2[j] <- aux2[j] + 2 * aux[(k + j)] * Conj(aux[k])
                }
            }
            ab <- rep(0, 2 * M)
            for (k in 1:M) {
                ab[k] = Re(aux2[k])
                ab[(k + M)] = -Im(aux2[k])
            }
        }
    }
    return(ab)
}
