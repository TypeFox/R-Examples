`Pevid.rel` <-
function(alleles, prob, x, u = NULL, k = c(1, 0, 0), S = NULL) {
    require(combinat)
    if (!is.numeric(x) || is.na(x) || x < 0) 
        stop("'x' must be a nonnegative integer")
    x0 <- x
    x <- round(x)
    if (x != x0) 
        warning("'x' has been rounded to integer") 
    if (!is.vector(alleles)) 
        stop("'alleles' must be a vector")
    if (!is.vector(prob)) 
        stop("'prob' must be a vector")
    n_a <- length(alleles)
    if (length(union(alleles, alleles)) < n_a)
        stop("'alleles' must be a vector of distinct elements")
    if (n_a != length(prob)) 
        stop("'alleles' and 'prob' must have the same length")
    if (!is.numeric(prob) || any(is.na(prob)) || any(prob <= 0) || any(prob >= 1)) 
        stop("all entries of 'prob' must be numbers between 0 and 1")
    if (sum(prob) > 1) 
        stop("sum of allele proportions is larger than 1")
    u <- union(u, u)
    n_u <- length(u)
    if (n_u > 0) {  
        if (setequal(intersect(u, alleles), u) == FALSE)
            stop("'u' must contain only elements from 'alleles'")
    }
    if (!is.null(S)) {
        if (is.genotype(S) == FALSE)
	    S <- as.genotype(S)
	if (any(is.na(S)) || length(S) != 1)
	    stop("entries of 'S' are not correct")
        if (x < 1)
            stop("'x' must be at least 1")
    }
    if (is.null(S) && x < 2)
        stop("'x' must be at least 2")
    if (!is.numeric(k) || any(is.na(k)) || any(k < 0) || any(k > 1)) 
        stop("all entries of 'k' must be numbers from the interval [0, 1]")
    if (length(k) != 3)
        stop("'k' should be a vector containing three kinship coefficients")
    if (sum(k) != 1) 
        stop("sum of kinship coefficients is different from 1")
    min.x <- trunc(n_u / 2) + n_u %% 2
    if (x > 0 && x < min.x)
        stop("not enough unknown contributors")
    ind <- NULL
    if (n_u > 0) {
        ind <- rep(NA, n_u)
        for (j in 1:n_u) 
            ind[j] <- which(alleles == u[j])
    }
    L <- function(r, uni = ind, p = prob) {
        nn <- length(uni)
	T <- rep(0, nn + 1)
        T[1] <- sum(p)^r
        if (nn > 1) {
            for (skip in 1:(nn - 1)) {
                K <- combn(uni, skip)
                Kcol <- ncol(K)
                t <- rep(0, Kcol)
                for (j in 1:Kcol)
                    t[j] <- sum(p[ -K[,j] ])^r
                T[skip + 1] <- (-1)^(skip) * sum(t)
            }
        }
        if (nn > 0)
            T[nn + 1] <- (-1)^(nn) * sum(p[ -uni ])^r
        return(sum(T))
    }
    if (is.null(S)) 
        res <- k[1]*L(r = 2*x) + k[2]*L(r = 2*x - 1) + k[3]*L(r = 2*x - 2) 	
    else if(heterozygote(S)) {                                            			
        I1u <- is.element(allele.names(S)[1], u) 	
        I2u <- is.element(allele.names(S)[2], u)
        I1d <- is.element(allele.names(S)[1], setdiff(alleles, u))
        I2d <- is.element(allele.names(S)[2], setdiff(alleles, u))
        ws1 <- which(alleles==allele.names(S)[1]) 	
        ws2 <- which(alleles==allele.names(S)[2])
        res <- k[1]*L(r = 2*x, uni = ind) + 
            k[2]*I1u*L(r = 2*x - 1, uni = setdiff(ind, ws1))/2 +
            k[2]*I1d*L(r = 2*x - 1)/2 +
	    k[2]*I2u*L(r = 2*x - 1, uni = setdiff(ind, ws2))/2 +
            k[2]*I2d*L(r = 2*x - 1)/2 +
	    k[3]*I1d*I2d*L(r = 2*x - 2) +
	    k[3]*(I1u*I2u - I1d*I2d)*
	    L(r = 2*x - 2, uni = setdiff(ind, c(I1u*ws1, I2u*ws2)))
    }
    else if (homozygote(S)) {                                            			
        Iu <- is.element(allele.names(S), u) 	
        Id <- is.element(allele.names(S), setdiff(alleles, u))
        ws <- which(alleles==allele.names(S)) 	
        res <- k[1]*L(r = 2*x, uni = ind) + 
            k[2]*Iu*L(r = 2*x - 1, uni = setdiff(ind, ws)) +
            k[2]*Id*L(r = 2*x - 1) +
	    k[3]*Id*L(r = 2*x - 2) +
	    k[3]*(Iu - Id)*
	    L(r = 2*x - 2, uni = setdiff(ind, Iu*ws))
    }
    return(res)
}	
    

