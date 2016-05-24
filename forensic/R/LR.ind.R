`LR.ind` <-
function(alleles, prob, x1, x2, u1 = NULL, u2 = NULL) {
    if (!is.numeric(x1) || any(is.na(x1)) ||  any(x1 < 0))
        stop("'x1' must be a vector of nonnegative integers")
    if (!is.numeric(x2) || any(is.na(x2)) || any( x2 < 0))
        stop("'x2' must be a vector of nonnegative integers")
    x10 <- x1
    x1 <- round(x1)
    if (any(x1 != x10)) 
        warning("elements of 'x1' have been rounded to integers")
    x20 <- x2
    x2 <- round(x2)
    if (any(x2 != x20)) 
        warning("elements of 'x2' have been rounded to integers")
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
        stop("the number of rows of 'prob' should equal to the number of alleles in the mixture")
    if (!is.numeric(prob) || any(is.na(prob)) || any(prob <= 0) || any(prob >= 1)) 
        stop("all entries of 'prob' must be numbers between 0 and 1")
    if (any(apply(prob, 2, sum) > 1)) 
        stop("sum of allele proportions in an ethnic group is larger than 1")
    if (length(x1) != dim(prob)[2]) 
        stop("error in the length of 'x1' or in the dimension of 'prob'")
    if (length(x2) != dim(prob)[2]) 
        stop("error in the length of 'x2' or in the dimension of 'prob'")
    u1 <- union(u1, u1)
    n_u1 <- length(u1)
    u2 <- union(u2, u2)
    n_u2 <- length(u2)
    if (n_u1 > 0) {
        if (setequal(intersect(u1, alleles), u1) == FALSE)
            stop("'u1' must contain only elements from 'alleles'")
    }
    if (n_u2 > 0) {
        if (setequal(intersect(u2, alleles), u2) == FALSE)
            stop("'u2' must contain only elements from 'alleles'")
    }
    if (setequal(alleles, u1) == TRUE)
        stop("under the null hypothesis there must be at least one known contributor")
    G <- length(x1)
    min.x1 <- trunc(n_u1 / 2) + n_u1 %% 2
    if (sum(x1) == 0 && setequal(u1, NULL) == FALSE)
        stop(gettextf(ngettext(min.x1,
            "there must be at least %d unknown person to carry alleles in 'u1'",
            "there must be at least %d unknown people to carry alleles in 'u1'"),
             min.x1))
    if (sum(x1) > 0 && sum(x1) < min.x1)
        stop("not enough unknown contributors under the null hypothesis")
    if (sum(x2) < sum(x1) + 1)
        stop("there should be more unknown contributors under the alternative hypothesis than under the null hypothesis")
    min.x2 <- trunc(n_u2 / 2) + n_u2 %% 2
    if (sum(x2) < min.x2)
        stop("not enough unknown contributors under the alternative hypothesis")
    if (setequal(intersect(u1, u2), u1) == FALSE)
        stop("'u1' must be a subset of 'u2'")
    ind1 <- NULL
    ind2 <- NULL
    if (n_u1 > 0) {
        ind1 <- rep(NA, n_u1)
        for (k in 1:n_u1) {
            ind1[k] <- which(alleles == u1[k])
        }
    }
    if (n_u2 > 0) {
        ind2 <- rep(NA, n_u2)
        for (k in 1:n_u2) {
            ind2[k] <- which(alleles == u2[k])
        }
    }
    p_evid <- function(p, x, ind, g = G) {
        require(combinat)
        n_u <- length(ind)
        T <- rep(0, n_u + 1)
        T[1] <- prod(apply(p, 2, sum)^(2 * x))
        if (n_u > 1) {
            for (skip in 1:(n_u - 1)) {
                K <- combn(ind, skip)
                Kcol <- ncol(K)
                t <- rep(0, Kcol)
                for (j in 1:Kcol)
                    t[j] <- prod(apply(matrix(p[ -K[,j], ], ncol = g), 2, sum)^(2 * x))
                T[skip + 1] <- (-1)^(skip) * sum(t)
            }
        }
        if (n_u > 0)
            T[n_u + 1] <- (-1)^(n_u) * prod(apply(matrix(p[ -ind, ], ncol = g), 2, sum)^(2 * x))
        return(sum(T))
    }
    p1 <- p_evid(p = prob, x = x1, ind = ind1)
    p2 <- p_evid(p = prob, x = x2, ind = ind2)
    LR <- p1 / p2
    return(LR)
}

