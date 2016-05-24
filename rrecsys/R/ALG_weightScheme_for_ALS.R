# weight function ####

weightScheme <- function(x, scheme, delta) {
    
    # uniform
    if (scheme == "uni") {
        W <- matrix(delta, nrow = nrow(x), ncol = ncol(x))
        W[x == 1] <- 1
        
        return(W)
        
    }
    
    # user oriented
    if (scheme == "uo") {
        W <- matrix(0, nrow = nrow(x), ncol = ncol(x))
        s <- apply(x, 1, sum)
        row_max_rated <- max(s)
        for (i in 1:nrow(x)) {
            # if (s[i]!=0){
            W[i, ] <- 1/(1 + exp(-s[i]/10)) + delta - 0.5
            # W[i,] <- delta + .6 * s[i]/row_max_rated
        }
        W[x == 1] <- 1
        
        return(W)
        
    }
    
    # item oriented
    if (scheme == "io") {
        W <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
        s <- (nrow(x) - apply(x, 2, sum))
        col_max_rated <- nrow(x)
        for (j in 1:ncol(x)) W[, j] <- delta + 0.6 * s[j]/col_max_rated
        W[x == 1] <- 1
        
        return(W)
        
    }
    
    # combined scheme
    regterm <- 0.3
    W <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    s <- apply(x, 1, sum)
    for (i in 1:nrow(x)) {
        if (s[i] != 0) {
            W[i, ] <- (1 - regterm) * (1/(1 + 1/s[i]) - 0.49)
            # W[i,] <- 1/s[i] + 0.01
        }
    }
    s <- (nrow(x) - apply(x, 2, sum))/nrow(x)
    m <- nrow(x)
    for (j in 1:ncol(x)) W[, j] <- W[, j] + regterm * 0.8 * s[j]
    W[x == 1] <- 1
    return(W)
} 
