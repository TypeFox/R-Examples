`Pmatch` <-
function(prob, k = c(1, 0, 0), theta = 0) {
    if (is.vector(prob) == TRUE) {
        if (length(prob) %% 2 == 0) 
	    prob <- matrix(prob, nrow = 2)
        else 
	    stop("'prob' must be a vector of odd length or a matrix with 2 rows")
    }
    prob <- as.matrix(prob)
    if (nrow(prob) != 2) 
        stop("'prob' must be a vector of odd length or a matrix with 2 rows")
    if (!is.numeric(prob) || any(is.na(prob)) || any(prob < 0) || any(prob >= 1))
        stop("all entries of 'prob' must be in the interval [0, 1)")
    if (any(apply(prob, 2, sum) > 1)) 
        stop("sum of allele proportions at a locus is larger than 1")
    if (any(apply(prob, 2, sum) == 0))
        stop("zero column in 'prob'")      
     if (!is.numeric(k) || any(is.na(k)) || any(k < 0) || any(k > 1)) 
        stop("all entries of 'k' must be numbers from the interval [0, 1]")
    if (length(k) != 3)
        stop("'k' should be a vector containing three kinship coefficients")
    if (sum(k) != 1) 
        stop("sum of kinship coefficients is different from 1")
    n <- ncol(prob) # number of loci
    dimnames(prob) <- list(c('allele 1', 'allele 2'), paste("locus ", 1:n, sep=""))
    M <- matrix(NA, nrow = 1, ncol = n)
    dimnames(M) <- list(NULL, paste("locus ", 1:n, sep=""))
    denom <- (1 + theta)*(1 + 2*theta)
    for (i in 1:n) {
        if (min(prob[, i]) == 0) {
	    p <- max(prob[1, i], prob[2, i]) 
	    M[1, i] <- k[1]*(2*theta + (1 - theta)*p)* 
	    (3*theta + (1 - theta)*p)/denom +
	    k[2]*(2*theta + (1 - theta)*p)/(1 + theta) +
	    k[3]
	}
        else 
            M[1, i] <- k[1]*2*(theta + (1 - theta)*prob[1, i])*
	    (theta + (1 - theta)*prob[2, i])/denom + 
	    k[2]*(2*theta + (1 - theta)*(prob[1, i] + prob[2, i]))/(2*(1 + theta)) + 
	    k[3]	
    }
    match.prod <- apply(M, 1, prod)
    z <- list(prob = prob, match = M, total_match = match.prod)
    z
}

