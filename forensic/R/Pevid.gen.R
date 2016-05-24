`Pevid.gen` <-
function(alleles, prob, x, T = NULL, V = NULL, theta = 0 ) {
    require(genetics)
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
    c <- length(alleles)
    if(c != length(union(alleles, alleles)))
        stop("allele duplicates in 'alleles'")
    if (c != length(prob))
        stop("'alleles' and 'prob' must have the same length")
    if (!is.numeric(prob) || any(is.na(prob)) || any(prob <= 0) || any(prob >= 1))
        stop("all entries of 'prob' must be numbers between 0 and 1")
    if (sum(prob) > 1) 
        stop("sum of allele proportions is larger than 1")
    if (is.null(T) == FALSE) {
        if (is.genotype(T) == FALSE)
            T <- as.genotype(T)
        if (any(is.na(T)))
            stop("entries of T are not correct")
    }      
    if (is.null(V) == FALSE) {
        if (is.genotype(V) == FALSE)
            V <- as.genotype(V)
        if (any(is.na(V)))
            stop("entries of V are not correct")
    }
    if (sum(prob) == 1 && setequal(union(alleles, allele.names(V)), alleles) == FALSE)
        stop("additional alleles in V (mixture contains all alleles of a locus)")
    if (setequal(intersect(allele.names(T), alleles), allele.names(T)) == FALSE)
        stop("unknown alleles in 'T'")  
    if (theta >= 1 || theta < 0) {
        stop("'theta' must be a number between 0 and 1, recommended 0.01 - 0.03") 
    }
    # the known number of declared contributors to the mixture
    n_T <- length(T)
    # the known number of people declared not to be contributors to the mixture 
    # (people who carry at least one allele from 'alleles')
    if (is.null(V) == FALSE) {
        f <- 0
        for (i in 1:length(V)) {
	    if (sum(carrier(V[i], alleles) == TRUE) > 0)
	        f <- f + 1
	}
	n_V <- f
    }
    else
        n_V <- 0
    # the known number of heterozygous declared contributors (in T)
    if (n_T > 0)
        h_T <- sum(heterozygote(T) == TRUE)  
    else 
        h_T <- 0  
    # the known number of heterozygous declared noncontributors (in V)
    if (n_V > 0)
        h_V <- sum(heterozygote(V) == TRUE)  
    else 
        h_V <- 0
    # the known number of copies of allele 'alleles[i]' in T
    t_i <- rep(0, c)
    if (n_T > 0) {
        for (i in 1:c) {
            t_i[i] <- sum(allele.count(T, alleles[i]))
        }
    }
    # the known number of copies of allele 'alleles[i]' in V
    v_i <- rep(0, c)
    if (n_V > 0) {
        for (i in 1:c) {
            v_i[i] <- sum(allele.count(V, alleles[i]))
        } 
    }
    # the known number of distinct alleles carried by n_T declared contributors
    t <- 0
    if (n_T > 0)
        t <- length(allele.names(T))
    # the known number of distinct alleles which must be in U
    u <- c - t
    # the known number of unknown contributors 
    n_U <- x
    # the known number of contributors
    n_C <- n_T + n_U
    min.n_U <- trunc(u / 2) + u %% 2
    min.n_C <- trunc(c / 2) + c %% 2 
    if (n_U < min.n_U || n_C < min.n_C)
        stop("not enough contributors")
    # the set of distinct alleles which must be in U:  U_0 = C \ T
    if (n_T > 0)
        U_0 <- setdiff(alleles, allele.names(T))
    else
        U_0 <- alleles
    # the known number of alleles in U that can be of arbitrary type from C
    r <- 2 * n_U - u
    perm <- function(r, c) {
        M <- NULL
        a <- rep(0, c)
        recurs <- function(r, c, S, i) {
            if (i == c - 1) {
                a[c] <<- r - S
                M <<- rbind(M, a, deparse.level = 0)
            }
            else {
                for (k in 0:r) {
                    if (S + k <= r) {
                        a[i + 1] <<- k
                        Recall(r, c, S + k, i + 1)
                    }
                    else break
                }
            }
        }
        recurs(r, c, 0, 0)
        M
    }
    # r_i[,i] ... the unknown number of copies of allele 'alleles[i]' among 
    # the r unconstrained alleles in U
    if (r>0) { 
        # all possibilities of the vector (r_1, r_2, ... , r_c), where  
        # sum(r_i)=r and r_i>=0, r_i<=r
        r_i <- perm(r, c)
    }
    else
        r_i <- matrix(rep(0, c), nrow = 1, byrow = TRUE)
    # u_i[,i] ... the unknown number of copies of allele 'alleles[i]' in U, 
    # sum(u_i) = 2*n_U 
    u_i <- matrix(0, nrow = nrow(r_i), ncol = c)
    if (n_T > 0) {
        for (i in 1:c) {
            if (sum(allele.count(T, alleles[i])) == 0)  
                u_i[, i] <- r_i[,i] + 1
            else    
                u_i[, i] <- r_i[, i]
        }
    }	   
    else
        u_i <- r_i + 1
    if (n_U == 0)
        const <- factorial(2 * n_U)
    else  
        const <- factorial(2 * n_U) / prod((1 - theta) + 
            ((2 * n_T + 2 * n_V):(2 * n_T + 2 * n_V + 2 * n_U - 1)) * theta) 
    results <- rep(0, nrow(u_i))
    for (d in 1:nrow(u_i)) { 
        prod_p <- rep(0, c)
        for (i in 1:c) {
            if (u_i[d, i] == 0)
	        prod_p[i] <- 1
	    else 
                prod_p[i] <- prod((1 - theta) * prob[i] + 
                    ((t_i[i] + v_i[i]):(t_i[i] + v_i[i] + u_i[d, i] -1)) * theta)  
        }
        results[d] <- prod(prod_p) / prod(factorial(u_i[d, ]))
    }
    return(const * sum(results))
}

