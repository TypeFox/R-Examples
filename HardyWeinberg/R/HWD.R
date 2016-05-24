HWD <- function(X) {
      if (length(X) != 3 | any(X < 0)) 
        stop("HWD: X is not a 3 by 1 non-negative count vector")
    if (any(!is.wholenumber(X))) {
        warning("HWD: Genotype counts are not integers, counts will be rounded.")
        X <- round(X, digits = 0)
    }
    n <- sum(X)   
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    X <- c(min(Xhom), Xhet, max(Xhom))
    nA <- 2 * X[1] + X[2]
    nB <- 2 * X[3] + X[2]
    pAA <- X[1]/n
    pA <- nA/(2*n)
    D <- pAA - pA^2
    names(D) <- NULL
    return(D)
}
