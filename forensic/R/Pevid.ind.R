`Pevid.ind` <-
function(alleles, prob, x, u = NULL) {
    require(combinat)
    if (!is.numeric(x) || any(is.na(x)) ||  any(x < 0)) 
        stop("'x' must be a vector of nonnegative integers")
    x0 <- x
    x <- round(x)
    if (any(x != x0)) 
        warning("elements of 'x' have been rounded to integers") 
    if (!is.vector(alleles)) 
        stop("'alleles' must be a vector")
    if (is.vector(prob) == TRUE)
        prob <- as.matrix(prob)
    if (!is.matrix(prob)) 
        stop("'prob' must be a matrix or a vector")
    n_a <- length(alleles)
    if (length(union(alleles, alleles)) < n_a)
        stop("'alleles' must be a vector of distinct elements")
    if (n_a != dim(prob)[1]) 
        stop("the number of rows of 'prob' should equal to the number of distinct alleles in the mixture")
    if (!is.numeric(prob) || any(is.na(prob)) || any(prob <= 0) || any(prob >= 1)) 
        stop("all entries of 'prob' must be numbers between 0 and 1")
    if (any(apply(prob, 2, sum) > 1)) 
        stop("sum of allele proportions in an ethnic group is larger than 1")
    if (length(x) != dim(prob)[2]) 
        stop("error in the length of 'x' or in the dimension of 'prob'")
    u <- union(u, u)
    n_u <- length(u)
    if (n_u > 0) {  
        if (setequal(intersect(u, alleles), u) == FALSE)
            stop("'u' must contain only elements from 'alleles'")
    }
    G <- length(x)
    min.x <- trunc(n_u / 2) + n_u %% 2
    if (sum(x) == 0 && setequal(u, NULL) == FALSE)
        stop(gettextf(ngettext(min.x,
  	    "there must be at least %d unknown person to carry alleles in 'u'",
	    "there must be at least %d unknown people to carry alleles in 'u'"),
	     min.x))
    if (sum(x) > 0 && sum(x) < min.x)
        stop("not enough unknown contributors")
    if (n_u > 0) {
        ind <- rep(NA, n_u)
        for (k in 1:n_u) 
            ind[k] <- which(alleles == u[k])
    }
    T <- rep(0, n_u + 1)
    T[1] <- prod(apply(prob, 2, sum)^(2 * x))
    if (n_u > 1) {
        for (skip in 1:(n_u - 1)) {
            K <- combn(ind, skip)
            Kcol <- ncol(K)
            t <- rep(0, Kcol)
            for (j in 1:Kcol)
                t[j] <- prod(apply(matrix(prob[ -K[,j], ], ncol = G), 2, sum)^(2 * x))
            T[skip + 1] <- (-1)^(skip) * sum(t)
        }
    }
    if (n_u > 0)
        T[n_u + 1] <- (-1)^(n_u) * prod(apply(matrix(prob[ -ind, ], ncol = G), 2, sum)^(2 * x))
    return(sum(T))
}

