## Randomely generate data conforms to multivariate response distribution 
## Author: Hua Zhou and Yiwen Zhang
##============================================================## 

##============================================================## 
## Generate multinomial data 
##============================================================##

rmn <- function(size, alpha, n) {
    ## alpha can be a vector or a matrix n is optional size can be a scalar;or a
    ## vector; seed is optional
    
    if (is.vector(alpha) && missing(n)) {
        stop("Error: when alpha is a vector, must give n.")
    }
    
    if (is.vector(alpha) && length(alpha) > 1) {
        ## devided by the last column of alpha, so that we can get rid of the
        ## identifiability issue when fitting.  This is equivalent to subtract the last
        ## column of the parameter matrix from the parameter matrix.
        alpha <- alpha/alpha[length(alpha)]
        alpha <- t(matrix(alpha, length(alpha), n))
        if (length(size) == 1) {
            size <- rep(size, n)
        } else {
            stop("The length of size variable doesn't match with alpha or beta")
        }
        
    } else if (missing(n)) {
        n <- nrow(alpha)
    } else if (length(n) != 1) {
        stop("n should be a scalar.")
    } else if (n != dim(alpha)[1]) {
        stop("The sizes of the input don't match. \n\t\t\talpha must be a vector or a matrix with n rows")
    }
    
    k <- ncol(alpha)
    if (k <= 1) {
        stop("The multivariate response data need to have more than one category.")
    }
    
    if (!is.vector(size)) {
        stop("Size must be a scalar, or a column vector matches with \n\t\t\tthe number of rows of alpha")
    } else if (length(size) != n) 
        stop("The length of size should match with n")
    
    ## ----------------------------------------## Generate data
    prob <- alpha/rowSums(alpha)
    
    rdm <- t(sapply(1:n, function(i, size, prob) return(rmultinom(n = 1, size = size[i], 
        prob = prob[i, ])), size, prob))
    
    return(rdm)
}



##============================================================## 
## Generate dirichlet multinomial data
##============================================================##

rdirm <- function(size, alpha, n) {
    ## alpha can be a vector or a matrix n is optional size can be a scalar;or a
    ## vector; seed is optional
    
    if (is.vector(alpha) && missing(n)) {
        stop("When alpha is a vector, must give n.")
    }
    
    if (is.vector(alpha) && length(alpha) > 1) {
        alpha <- matrix(alpha, n, length(alpha), byrow = TRUE)
        if (length(size) == 1) {
            size <- rep(size, n)
        } else if (length(size) != n) {
            stop("The length of size variable doesn't match with sample size")
        }
        
    } else if (missing(n)) {
        n <- nrow(alpha)
    } else if (length(n) != 1) {
        stop("n should be a scalar.")
    } else if (n != nrow(alpha)) {
        stop("The sizes of the input don't match. \n\t\talpha must be a vector or a matrix with n rows")
    }
    
    k <- ncol(alpha)
    if (k <= 1) {
        stop("The multivariate response data need to have more than one category.")
    }
    
    if (!is.vector(size)) {
        stop("Size must be a scalar, or a column vector matches with \n\t\tthe number of rows of alpha")
    } else if (length(size) != n) 
        stop("The length of size should match with n")
    
    ##----------------------------------------## 
    ## Generate data
    ##----------------------------------------## 
    
    G <- sapply(c(alpha), "rgamma", n = 1)
    G <- matrix(G, nrow = n)
    prob <- G/rowSums(G)
    ridx <- rowSums(G) == 0
    if (any(ridx)) {
        if (sum(ridx) > 1) {
            prob[ridx, ] <- t(apply((alpha[ridx, ]/rowSums(alpha[ridx, ])), 1, function(x) return(rmultinom(n = 1, 
                size = 1, prob = x))))
        } else {
            prob[ridx, ] <- rmultinom(n = 1, size = 1, alpha[ridx, ]/sum(alpha[ridx, 
                ]))
        }
    }
    rdm <- t(sapply(1:n, function(i, size, prob) return(rmultinom(n = 1, size = size[i], 
        prob = prob[i, ])), size, prob))
    
    return(rdm)
}


##============================================================## 
## Generate GDM data 
##============================================================##

rgdirm <- function(size, alpha, beta, n) {
    
    if (length(alpha) != length(beta)) 
        stop("The size of alpha and beta should match.")
    
    if (is.vector(alpha) && missing(n)) 
        stop("When alpha and beta are vectors, must give n.")
    
    if (is.vector(alpha)) {
        alpha <- t(matrix(alpha, length(alpha), n))
        beta <- t(matrix(beta, length(beta), n))
        if (length(size) == 1) {
            size <- rep(size, n)
        } else {
            stop("The length of size variable doesn't match with alpha or beta")
        }
    } else {
        if (missing(n)) {
            n <- dim(alpha)[1]
        } else if (length(n) != 1) {
            stop("n should be a scalar.")
        } else if (n != dim(alpha)[1]) 
            stop("The sizes of the input alpha don't match with n")
    }
    
    k <- dim(alpha)[2]
    if (k < 1) {
        stop("The multivariate response data need to have more than one category.")
    }
    
    if (!is.vector(size)) {
        stop("n must be a scalar, or a column vector matches with \n\t\t\tthe number of rows of alpha")
    } else if (length(size) != n) 
        stop("The length of size should match with n")
    
    ##----------------------------------------## 
    ## Generate data
    ##----------------------------------------## 
    
    rdm <- matrix(0, n, (k + 1))
    rdm[, 1] = size
    
    for (i in 1:k) {
        rdm[, c(i, i + 1)] <- rdirm(size = rdm[, i], alpha = cbind(alpha[, i], beta[, 
            i]), n)
    }
    
    return(rdm)
}


##============================================================## 
## Generate NegMN data 
##============================================================##

rnegmn <- function(prob, beta, n) {
    
    if (is.vector(prob) && missing(n)) {
        stop("When prob is a vector, must give n.")
    }
    
    if (is.vector(prob) && length(prob) > 1) {
        prob <- matrix(prob, n, length(prob), byrow = TRUE)
        if (length(beta) == 1) {
            beta <- rep(beta, n)
        } else {
            stop("The length of beta doesn't match with the size of prob")
        }
    } else if (missing(n)) {
        n <- nrow(prob)
    } else if (length(n) != 1) {
        stop("n should be a scalar")
    } else if (n != nrow(prob)) {
        stop("The sizes of the input don't match. \n\t\t\tprob must be a vector or a matrix with n rows")
    }
    
    k <- ncol(prob)
    if (k <= 1) 
        stop("The multivariate response data need to have more than one category.")
    
    if (nrow(prob) != length(beta)) 
        stop("The dimension of prob and the length of beta doesn't match")
    
    if (length(beta) != n) 
        stop("The length of beta should match with n")
    
    ##----------------------------------------## 
    ## Generate data
    ##----------------------------------------## 
    probbeta <- 1 - rowSums(prob)
    prob <- cbind(prob, probbeta)
    scale <- 1/probbeta - 1
    G <- sapply(1:n, function(i, B, P) return(rgamma(n = 1, shape = B[i], rate = 1/P[i])), 
        B = beta, P = scale)
    lambda <- prob[, 1:k] * G/(1 - probbeta)
    rdm <- matrix(sapply(lambda, "rpois", n = 1), n, k)
    
    return(rdm)
} 
