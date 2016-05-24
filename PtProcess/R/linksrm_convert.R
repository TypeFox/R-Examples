linksrm_convert <- function (params, abc = TRUE) 
{
    N <- sqrt(length(params) + 1) - 1
    if (abc == FALSE) {
        alpha <- params[1:N]
        nu <- params[(N + 1):(2 * N)]
        rho <- params[(2 * N + 1):(3 * N)]
        theta <- diag(N)
        theta[diag(N) == 0] <- params[(3 * N + 1):(N^2 + 2 * 
            N)]
        theta <- t(theta)
        A <- alpha
        B <- nu * rho
        if (length(rho) == 1) 
            cc <- theta/rho
        else cc <- diag(1/rho) %*% theta
    }
    if (abc == TRUE) {
        A <- params[1:N]
        B <- params[(N + 1):(2 * N)]
        cc <- matrix(params[(2 * N + 1):(N^2 + 2 * N)], byrow = TRUE, 
            ncol = N)
        alpha <- A
        rho <- 1/diag(cc)
        nu <- B/rho
        if (length(rho) == 1) 
            theta <- rho * cc
        else theta <- diag(rho) %*% cc
    }
    return(list(params = params, a = A, b = B, cc = cc, alpha = alpha, 
        nu = nu, rho = rho, theta = theta))
}

