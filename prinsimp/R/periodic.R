### Returns the lambda matrix, $\Lambda$, to calculate the periodic
### simplicity measure of a vector $v \in R^d$ via $v^T \Lambda\ v$.
###
### Raises an error if environment levels across periods don't line up
### (i.e., don't advance by the same step in all periods)
###
### Arguments:
###   x - vector containing ordered values of the environment levels
###   period - length of one period 
lambda_periodic <- function(x, period) {
    if (!check_period(x, period)) {
        stop('Measurement steps are not consistent across periods.')
    }
    
    d <- length(x)
    I <- diag(d)
    M <- M_periodic(d, period)
    Q <- Q_periodic(d, period)
    
    diag(d) - crossprod(I - t(M) %*% Q %*% M)
}


M_periodic <- function(d, period) {
    k <- d %% period                         # number of complete periods
    l <- d %/% period                        # remaining periods

    MM <- diag(period)
    cbind(matrix(rep(MM, l), nrow = period),
          MM[ , seq_len(k)])
}


Q_periodic <- function(d, period) {
    k <- d %% period                         # number of complete periods
    l <- d %/% period                        # remaining periods

    diag(1 / rep(c(l+1, l), c(k, period-k)))
}

check_period <- function(x, period) {
    all(diff(diff(x, period)) == 0)
}
