#' Returns the posterior density \eqn{p(\rho|\Theta)}
#' 
#' @param x an array of dimension \code{time x variable x case} of terms.
#' @param y a matrix of dimension \code{time x case} of responses.
#' @param rho vector of quantiles.
#' @param log.p if \code{TRUE}, probabilities are given as \code{log(p)}.
#' @export
p_rho <- function(x, y, rho, log.p = FALSE) {
    
    T <- nrow(y) - 1                 # initial values are not modelled
    N <- ncol(y)
    K <- dim(x)[2]

    ## measure X_i from agent-specific means for each column
    x <- center_x(x)
    
    ## likelihood is calculated for (Y_it - Y_i0)
    y <- center_y(y)

    density <- numeric(length(rho))
    for (i in seq_along(rho)) {
        density[i] <- N * b(rho[i], T) - ((N * (T-1) - K) / 2) *
            log(Q_star(x, w(y, rho[i])))
    
    }

    if (log.p) {
        density
    } else {
        exp(density)
    }
}


## Centers X_i from agent-specific means for each column (dropping the first wave)
center_x <- function(x) {
    T <- nrow(x) - 1                 # initial values are not modelled
    K <- ncol(x)
    N <- dim(x)[3]
    
    ## initial values are not modelled
    x_0 <- matrix(x[1,,], K, N)
    x <- x[-1, , , drop = FALSE]
    
    ## measure X_i from agent-specific means for each column
    sweep(x, 2:3, apply(x, 3, colMeans))
}


## Centers Y_i relative to the first wave (dropping it)
center_y <- function(y) {
    ## initial values are not modelled
    y_0 <- y[1,]
    y <- y[-1,]
    
    ## likelihood is calculated for (Y_it - Y_i0)
    sweep(y, 2, y_0)
}


## Returns $w_i = y_i - \rho y_{i-}
w <- function(y, rho) {
    T <- nrow(y)
    N <- ncol(y)

    w <- matrix(0, T, N)
    for (i in seq_len(N)) {
        w[,i] <- y[,i] - rho * y_(y, i)
    }
    w
}


## Returns $y_{i-} = ( y_{i0}, y_{i1}, \dots, y_{i, T-1} ) $
y_ <- function(y, i) {
    T <- nrow(y)
    
    c(0, y[-T, i])
}

## Returns $\frac{1}{T} \sum_{t=1}^{T-1} \frac{T-t}{t} \rho^t$
b <- function(rho, T) {
    t <- seq_len(T-1)
    sum((T - t) * rho^t / t) / T
}


## Returns $\sum_i w'_i H w_i - \sum_i w'_i H x_i \left(
##    \sum_i x'_i H x_i \right)^{-1} \sum_i x'_i H w_i$
Q_star <- function(x, w) {
    wHw(w) - wHx(x, w) %*% solve(Hstar(x)) %*% xHw(x, w)
}


## Returns $H = I_T - (1/T)_{J J'}$
h <- function(n) {
    diag(n) - matrix(1/n, n, n)
}


## Returns $\sum_i w_i' H w_i$
wHw <- function(w) {
    T <- nrow(w)
    N <- ncol(w)

    H <- h(T)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(w[,i]) %*% H %*% w[,i]
    }
    sum
}


## Returns $\sum_i w_i' H x_i$
wHx <- function(x, w) {
    T <- nrow(w)
    N <- ncol(w)

    H <- h(T)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(w[,i]) %*% H %*% x[,,i]
    }
    sum
}


## Returns $\sum_i x_i' H w_i$
xHw <- function(x, w) {
    T <- nrow(w)
    N <- ncol(w)

    H <- h(T)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(x[,,i]) %*% H %*% w[,i]
    }
    sum
}


## Returns $(\sum_i x_i' H x_i)^{-1}$
Hstar <- function(x) {
    T <- dim(x)[1]
    N <- dim(x)[3]

    H <- h(T)

    sum <- 0
    for (i in seq_len(N)) {
        sum <- sum + t(x[,,i]) %*% H %*% x[,,i]
    }
    sum
}
