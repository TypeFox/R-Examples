"permutations" <- function(n) { 
    ##function permutations from package e1071 
    if (n == 1) 
        return(matrix(1)) 
    else if (n < 2) 
        stop("n must be a positive integer") 
    z <- matrix(1) 
    for (i in 2:n) { 
        x <- cbind(z, i) 
        a <- c(1:i, 1:(i - 1)) 
        z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x)) 
        z[1:nrow(x), ] <- x 
        for (j in 2:i - 1) { 
            z[j * nrow(x) + 1:nrow(x), ] <- x[, a[1:i + j]] 
        } 
    } 
    dimnames(z) <- NULL 
    z 
} 
