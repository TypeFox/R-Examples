# Functions to calculate log of pdf of each of the multivariate models.  
# Author: Yiwen Zhang
##============================================================## 

##============================================================## 
## Multinomial model 
##============================================================##
dmultn <- function(X, Y, B1, weight) {
    N <- nrow(X)
    d <- ncol(Y)
    batch_sizes <- rowSums(Y)
    B <- cbind(B1, 0)
    alpha <- exp(X %*% B)
    return(sum(Y * (X %*% B) * weight) - sum(Y * log(rowSums(alpha) + 1) * weight) + 
        sum(lgamma(batch_sizes + 1) * weight) - sum(lgamma(Y + 1) * weight))
}


dmn <- function(Y, prob) {
    
    if (is.vector(prob) && length(prob) > 1) {
        
        if (is.vector(Y) && length(Y) > 1) {
            
            if (length(Y) != length(prob)) {
                stop("size of Y and prob doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        
        prob <- t(matrix(prob, length(prob), dim(Y)[1]))
        
    }
    
    if (any(dim(prob) != dim(Y))) 
        stop("dimensions of prob and Y do not match")
    
    m <- rowSums(Y)
    logl <- lgamma(m + 1) - rowSums(lgamma(Y + 1)) + rowSums(Y * log(prob))
    
    return(logl)
    
}


## ============================================================## 
## Dirichlet multinomial model
## ============================================================## 
## this alpha could be X%*%alpha_hat in a regression estimate

ddirm <- function(Y, alpha) {
    
    
    if (is.vector(alpha) && length(alpha) > 1) {
        
        if (is.vector(Y) && length(Y) > 1) {
            if (length(Y) != length(alpha)) {
                stop("size of Y and alpha doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        alpha <- t(matrix(alpha, length(alpha), dim(Y)[1]))
    }
    if (any(dim(alpha) != dim(Y))) {
        stop("dimensions of alpha and Y do not match")
    }
    
    canc <- rowSums(alpha) > 1e+08
    
    ## ----------------------------------------## 
    ## Calculate
    ## ----------------------------------------## 
    alpha_rowsums <- rowSums(alpha)
    m <- rowSums(Y)
    
    logl <- lgamma(m + 1) + rowSums(lgamma(Y + alpha)) + lgamma(alpha_rowsums) - 
        rowSums(lgamma(Y + 1)) - rowSums(lgamma(alpha)) - lgamma(alpha_rowsums + 
        m)
    
    if (sum(canc) > 0) {
        # logMN <- rowSums(Y*log(alpha)) - rowSums(Y*log(alpha_rowsums)) + lgamma(m+1) -
        # rowSums(lgamma(Y+1)) logl[canc] <- logMN[canc]
        logl[canc] <- -Inf
    }
    
    return(logl)
}


## ============================================================## 
## Generalized dirchlet multinomial model
## ============================================================##

dgdirm <- function(Y, alpha, beta) {
    
    if (is.vector(alpha) && is.vector(beta)) {
        
        if (length(alpha) != length(beta)) 
            stop("the sizes of alpha and beta do not match")
        
        if (is.vector(Y) && length(Y) > 1) {
            
            if (length(Y) != (length(alpha) + 1)) {
                stop("size of Y and alpha doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        
        alpha <- t(matrix(alpha, length(alpha), dim(Y)[1]))
        beta <- t(matrix(beta, length(beta), dim(Y)[1]))
    }
    
    if (nrow(alpha) != nrow(Y) | ncol(alpha) != ncol(Y) - 1) 
        stop("dimensions of the parameters and Y do not match")
    
    ## ----------------------------------------## 
    ## Calculate
    ## ----------------------------------------## 
    
    m <- rowSums(Y)
    d <- ncol(Y)
    
    z <- t(apply(apply(apply(Y, 1, rev), 2, cumsum), 2, rev))
    
    logl <- lgamma(m + 1) - rowSums(lgamma(Y + 1)) + rowSums(lgamma(Y[, -d] + alpha)) - 
        rowSums(lgamma(alpha)) + rowSums(lgamma(z[, -1] + beta)) - rowSums(lgamma(beta)) - 
        rowSums(lgamma(alpha + beta + z[, -d])) + rowSums(lgamma(alpha + beta))
    
    return(logl)
}


## ============================================================## 
## Negative Multinomial model pdf calculation
## ============================================================##

dnegmn <- function(Y, prob, beta) {
    if (is.vector(prob) && length(prob) <= 1) 
        stop("The length of alpha should be larger than 2.")
    
    if (is.vector(prob) && length(prob) > 1) {
        if (length(beta) != 1) 
            stop("The sizes of alpha and beta do not match")
        
        if (is.vector(Y) && length(Y) > 1) {
            
            if (length(Y) != length(prob)) {
                stop("size of Y and alpha doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        
        prob <- matrix(prob, nrow(Y), length(prob), byrow = TRUE)
        beta <- matrix(beta, nrow(Y), 1)
    }
    beta <- matrix(beta, , 1)
    if (nrow(beta) != nrow(Y) | nrow(prob) != nrow(Y)) 
        stop("dimensions of the parameters and Y do not match")
    if (any(rowSums(prob) >= 1)) 
        stop("sum of probability should be smaller than 1")
    ## ----------------------------------------## 
    ## Calculate the log ln
    ## ----------------------------------------## 
    
    m <- rowSums(Y)
    d <- ncol(Y)
    logl <- lgamma(beta + rowSums(Y)) - lgamma(beta) - rowSums(lgamma(Y + 1)) + rowSums(Y * 
        log(prob)) + beta * log(1 - rowSums(prob))
    return(logl)
}


## ============================================================## 
## Negative Multinomial model
## ============================================================##

dneg <- function(Y, alpha, beta) {
    
    if (is.vector(alpha) && length(alpha) <= 1) 
        stop("The length of alpha should be larger than 2.")
    
    if (is.vector(alpha) && length(alpha) > 1) {
        if (length(beta) != 1) 
            stop("The sizes of alpha and beta do not match")
        
        if (is.vector(Y) && length(Y) > 1) {
            
            if (length(Y) != length(alpha)) {
                stop("size of Y and alpha doesn't match.")
            } else {
                Y <- matrix(Y, 1, length(Y))
            }
            
        } else if (is.vector(Y) && length(Y) <= 1) {
            stop("Y can not be a scalar")
        }
        
        alpha <- matrix(alpha, nrow(Y), length(alpha), byrow = TRUE)
        beta <- matrix(beta, nrow(Y), 1)
    }
    beta <- matrix(beta, , 1)
    if (nrow(beta) != nrow(Y) | nrow(alpha) != nrow(Y)) 
        stop("dimensions of the parameters and Y do not match")
    
    ## ----------------------------------------## 
    ## Calculate the log ln
    ## ----------------------------------------## 
    m <- rowSums(Y)
    d <- ncol(Y)
    logl <- lgamma(beta + rowSums(Y)) - lgamma(beta) - rowSums(lgamma(Y + 1)) + rowSums(Y * 
        log(alpha)) - (beta + m) * log(rowSums(alpha) + 1)
    
    return(logl)
} 
