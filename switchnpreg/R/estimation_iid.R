## Estimate the p_ij matrix, p_ij = p(y_i|z_i = j, Theta); that is,
## the probabilities of y_i belonging to component j.
pij_iid <- function(y, current) {
    N <- length(y)
    
    alpha <- current$alpha
    f <- current$f
    sigma2 <- current$sigma2
    
    ps <- rep(alpha, each = N) *
        dnorm(y, mean = f, sd = rep(sqrt(sigma2), each = N))
    ps / rowSums(ps)
}


## Estimate alpha
alpha_iid <- function(y, pij, current) {
    N <- nrow(pij)
    colSums(pij) / N
}


criteria_alpha_iid <- function(current, update) {
    abs(update - current)
}


stderr_iid <- function(y, current) {
    pij <- current$pij
    alpha <- current$alpha

    N <- nrow(pij)
    J <- ncol(pij)
    
    if (J==2) {
        pi1 <- pij[,1]
        p1  <- alpha[1]
        Iy_hat <- sum((pi1/p1 - (1-pi1)/(1-p1))^2)
        return(sqrt(1/(Iy_hat)))
    } else {
        ## expected value of the negatives of the 2nd derivatives
        EB <- matrix((N/(1-sum(alpha[-J]))), J-1, J-1)
        diag(EB) <- sapply(alpha[-J], function(x) {
            N * (1/x + 1/(1-sum(alpha[-J])))
        })
        
        ## expected value of L2'L2'^T, gradient x transposed gradient
        ES2 <- matrix(0, nrow = J-1, ncol = J-1)
        for (j in 1:(J-1)) {
            for (l in 1:(J-1)) {
                if (j == l) {
                    ES2[j, l] <- N * ((1/alpha[j]) + (1/(1-sum(alpha[-J])))) -
                        sum(((pij[,j]/alpha[j]) -
                             ((1-rowSums(pij[,-J])) / (1-sum(alpha[-J]))))^2)
                } else {
                    ES2[j, l] <- N / (1-sum(alpha[-J])) -
                        sum(((pij[,j]/alpha[j]) -
                             ((1-rowSums(pij[,-J])) / (1-sum(alpha[-J])))) *
                            ((pij[,l]/alpha[l]) -
                             ((1-rowSums(pij[,-J])) / (1-sum(alpha[-J])))))
                }
            }
        }
        I_hat = EB - ES2
        I_inv <- solve(I_hat)

        return(c(sqrt(diag(I_inv)),
                 sqrt(sum(diag(I_inv)) + 2*sum(I_inv[lower.tri(I_inv)]))))
    }
}
