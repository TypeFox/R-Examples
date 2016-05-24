HWf <- function (X) 
{
    if (length(X) != 3 | any(X < 0)) 
        stop("HWf: X is not a 3 by 1 non-negative count vector")
    if (any(!is.wholenumber(X))) {
        warning("Genotype counts are not integers, counts will be rounded.")
        X <- round(X, digits = 0)
    }
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
     # reorganize counts in standard AA, AB, BB format.
    X <- c(min(Xhom), Xhet, max(Xhom))
    if(maf(X)==0) {
        warning("Monomorphic marker, f set to zero.")
        f <- 0
    } else {
       nA <- 2 * X[1] + X[2]
       nB <- 2 * X[3] + X[2]
       f <- (4 * X[1] * X[3] - (X[2])^2)/(nA * nB)
       names(f) <- NULL
   }
    return(f)
}


