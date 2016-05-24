checkgenZ <-
function (genZ, tol = 1e-04) 
{
    if (ncol(genZ)%%3 != 0) {
        stop("The number of columns in genZ must be a multiple of 3")
    }
    ok <- 1
    for (i in 1:(nrow(genZ))) {
        for (l in 1:(ncol(genZ)/3)) {
            s <- sum(genZ[i, (3 * l - 2):(3 * l)])
            if (abs(s - 1) > tol) {
                warning("Indiv", i, "; Locus", l, "; sum =", 
                  s, ", should be 1")
                ok <- 0
            }
        }
    }
    if (ok != 1) {
        stop()
    }
}
