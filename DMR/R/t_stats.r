t_stats <- function(model, M, ind1, ind2, sigma_sq, z){
    if (ind1 != ind2){
        n <- nrow(M)
        Tmat <- matrix(0, n, n)
        ind <- t(combn(1:n, 2))
        x <- as.matrix(M[ind[, 1],] - M[ind[, 2], ])
        if (ncol(x) == 1) x <- t(x)
        Tmat[ind] <- ((x%*%z)^2)/(sigma_sq*(apply(x, 1, function(y) t(y)%*%y)))
        Tmat <- cbind(0, Tmat)
        t_st <- c(0, (as.numeric(summary(model)$coeff[ind1:ind2,3]))^2)
        Tmat <- rbind(t_st, Tmat)
    } else {
        Tmat <- matrix(c(0, 0, (as.numeric(summary(model)$coeff[ind1,3]))^2,0), 2, 2)
    }
    return(Tmat)
}
