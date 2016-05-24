nntsABcoefficients <-
function (cpars = 1/sqrt(2*pi), M = 0) 
{
    if (M == 0) 
        return("Uniform case: No coefficients AB")
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")

    if (M > 0) {
aux <- complex(M)
        for (j in 1:M) {
            for (k in 1:(M - j + 1)) {
            aux[j] <- aux[j] + 2 * cpars[(k + j)] * Conj(cpars[k])
            }
        }
        ab <- rep(0, 2 * M)
        for (k in 1:M) {
            ab[k] = Re(aux[k])
            ab[(k + M)] = -Im(aux[k])
        }
    }
    return(ab)
}

