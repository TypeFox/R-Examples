"ffFullMatrix" <-
function (X, x, maxInt, blk = NULL) 
{
    if (!is.data.frame(X)) 
        X <- as.data.frame(X)
    Z <- data.frame(one = rep(1, nrow(X)))
    k <- 0
    if (!is.null(blk)) {
        blk <- as.matrix(blk)
        if (nrow(blk) != nrow(X)) 
            stop("Matrix and block should have the same number of rows.")
        k <- ncol(blk)
        Z <- cbind(Z, blk)
        names(Z)[-1] <- paste("bk", seq(k), sep = "")
    }
    ord <- min(length(x), maxInt)
    nT <- rep(0, ord)
    for (i in seq(ord)) {
        tt <- subsets(length(x), i, x)
        if (is.null(dim(tt))) 
            tt <- matrix(tt, nrow = 1)
        for (j in 1:nrow(tt)) {
            nT[i] <- nT[i] + 1
            Z <- cbind(Z, eval(parse(text = (paste("X[", tt[j,  
                ],  "]", collapse = "*", sep = "")))))
            names(Z)[ncol(Z)] <- paste("x", tt[j, ], collapse = "*", 
                sep = "")
        }
    }
    nT <- c(k, nT)
    if (ord > 1) 
        names(nT) <- c("blk", "main", paste("int", 2:ord, sep = "."))
    else names(nT) <- c("blk", "main")
    list(Xa = as.matrix(Z), x = x, maxInt = ord, nTerms = nT)
}
