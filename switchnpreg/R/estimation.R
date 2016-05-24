#################################################
## Functions to estimate the parameters -      ##
## commom to both penalized and Bayesian cases ##
#################################################

## Apply f_j update function for each j=1..J
updateF <- function(x, y, pij, current, lambda, fFn, ...){
    N <- nrow(pij)
    J <- ncol(pij)

    sigma2 <- current$sigma2
    
    f_hat <- matrix(nrow = N, ncol = J)
    trace <- numeric(J)
    for (j in seq(J)) {
        fj <- fFn(x = x,
                        y = y,
                        pij = pij[,j],
                        sigma2 = sigma2[j],
                        lambda = lambda[j],
                        ...)
        f_hat[,j] <- fj$f_hat
        trace[j] <- fj$trace
    }
    list(f_hat = f_hat,
         traces = trace)
}


### Estimate sigma assuming that each component has equal variance
equal.sigma2.update <- function(y, pij, f, traces) {
    J <- ncol(pij)
    sigma2 <- sum(colSums((((y-f)^2)*pij))) / (sum(pij) - sum(traces))
    rep(sigma2, times = J)
}


### Estimate sigma without assuming that each component has equal
### variance
diff.sigma2.update <- function(y, pij, f, traces) {
    colSums(((y-f)^2)*pij) / (colSums(pij) - traces)
}


## Calculate the convergence criteria in the M step
compare <- function(current, updates, criteria.alpha.function){
    c(c.f = colMeans(abs(updates$f - current$f)),
      c.sigma2 = abs(updates$sigma2-current$sigma2),
      c.alpha = criteria.alpha.function(current$alpha,
                                      updates$alpha))
}


### New values of smoothing parameters after the M step
updateLambda <- function(fjFn, x, y, current, ...) {
    J <- ncol(current$pij)
    best_lambda <- numeric(J)
    for (j in seq(J)) {
        optimum <- optimize(GCV,
                            maximum=FALSE,
                            x = x,
                            y=y,
                            fjFn = fjFn,
                            sigma2=current$sigma2[j],
                            pij=current$pij[, j], ...)
        best_lambda[j] <- optimum$minimum
    }
    best_lambda
}



### Cross-validation criterion for smoothing parameter lambda
###
### GCV(lambda_j) is defined as:
###     1/n sum_i{1:n} p_ij (y_i - f_j(x_i))^2 / (1-Hj_ii)^2
###
### fjFn - fj-update function
### x, y - data
### lambda, sigma2, pij - curent parameter values
### ... - optional arguments to `fjFn`
GCV <- function(fjFn, lambda, x, y, sigma2, pij, ...) {
    fj <- fjFn(x, y, pij, sigma2, lambda, ...)
    fhat <- fj$f_hat
    H <- fj$H

    gcv_i <- ((y-fhat) / (1-diag(H)))^2
    
    N <- length(y)
    (1/N) * sum(gcv_i * pij)
}
