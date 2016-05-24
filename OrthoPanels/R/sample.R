## Draw `n` samples of rho, v, and beta parameters from the
## distribution defined by the data `x` and `y`, with possible values
## of rho defined by the set `pts`
##' @importFrom MASS mvrnorm
##' @importFrom stats rgamma
sample_all <- function(x, y, n, pts = seq(-.995, .995, .001)) {
    T <- nrow(y) - 1                 # initial values are not modelled
    N <- ncol(y)
    K <- dim(x)[2]

    ## measure X_i from agent-specific means for each column
    x <- center_x(x)
    
    ## likelihood is calculated for (Y_it - Y_i0)
    y <- center_y(y)

    density <- numeric(length(pts))
    Q_star <- numeric(length(pts))
    H_ast <- Hstar(x)
    beta_hat <- list(length(pts))
    for (i in seq_along(pts)) {
        W <- w(y, pts[i])
        WHW <- wHw(W)
        WHX <- wHx(x, W)
        XHW <- xHw(x, W)

        beta_hat[[i]] <- solve(H_ast, XHW)
        Q_star[i] <- WHW - WHX %*% beta_hat[[i]]
        density[i] <- N * b(pts[i], T) - ((N * (T-1) - K) / 2) *
            log(Q_star[i])
    }

    density <- exp(density - max(density))
    is <- sample(seq_along(pts), n, replace = TRUE, prob = density / sum(density))
    rhos <- pts[is]
    vs <- numeric(n)
    betas <- matrix(, K, n)
    if (K > 1) {
        rownames(betas) <- colnames(x)
    }
    H_ast_inv <- solve(H_ast)
    for (i in seq(n)) {
        j <- is[i]

        vs[i] <- rgamma(1, shape = (N*(T-1) - K)/2,
                        rate = Q_star[j]/2)
        betas[, i] <- mvrnorm(1, beta_hat[[j]], H_ast_inv/vs[i])
    }
    
    list(rho = rhos,
         sig2 = 1/vs,
         beta = t(betas))
}


## Draw `n` samples from `rho`s, proportional to the posterior pdf(rho | x, y)
sample_rho <- function(n, x, y, rho) {
    ps <- p_rho(x, y, rho, log.p = TRUE)
    ps <- exp(ps - max(ps))
    sample(rho, n, replace = TRUE, prob = ps / sum(ps))
}


## Sample variance from the gamma distribution with parameters defined
## by the data `x` and `y` and given a sample of `rho`s
##' @importFrom stats rgamma
sample_sig <- function(x, y, rho) {
    T <- nrow(y) - 1
    N <- ncol(y)
    K <- dim(x)[2]

    x <- center_x(x)
    y <- center_y(y)
    
    sig <- numeric(length(rho))
    rho_rle <- rle(sort(rho))
    j <- 1
    for (i in seq_along(rho_rle$values)) {
        n <- rho_rle$lengths[i]         # number of draws
        sig[j:(j+n-1)] <- rgamma(n, shape = (N*(T-1) - K)/2,
                                 rate = Q_star(x, w(y, rho_rle$values[i]))/2)
        j <- j + n
    }
    sig[order(order(rho))]
}


## Sample beta from the normal distribution with parameters defined by
## the data `x` and `y` and given a sample of `rho`s and variance `v`s
##' @importFrom MASS mvrnorm
sample_beta <- function(x, y, rho, v) {
    stopifnot(length(rho) == length(v))
    
    x <- center_x(x)
    y <- center_y(y)

    K <- dim(x)[2]
    
    h_ast_inv <- solve(Hstar(x))
    beta <- matrix(, K, length(rho))
    
    rho_rle <- rle(sort(rho))
    irho <- order(rho)
    j <- 1
    
    for (i in seq_along(rho_rle$values)) {
        n <- rho_rle$lengths[i]         # number of draws
        beta_hat <- h_ast_inv %*% xHw(x, w(y, rho_rle$values[i]))
        
        for (k in 0:(n-1)) {
            vv <- v[irho[j+k]]
            beta[ , j+k] <- mvrnorm(1, beta_hat, h_ast_inv/vv)
        }
        j <- j + n
    }
    
    t(beta[ , order(irho), drop = FALSE])
}
