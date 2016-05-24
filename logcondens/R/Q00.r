Q00 <- function(x, a, u, v, gamma, QFhat = FALSE){

    ## x, gamma in R
    ## a, u, v in R^m
    q <- exp(a * (x - u) + a ^ 2 * gamma ^ 2 / 2) *
        (pnorm(v - x - a * gamma ^ 2, mean = 0, sd = gamma) - 
         pnorm(u - x - a * gamma ^ 2, mean = 0, sd = gamma))
         
    # Robustification of q:
    tmp.m <- (0.5 * (u + v) - x) / gamma - a * gamma
    tmp.d <- 0.5 * (v - u) / gamma
    tmp.q1 <- exp(a * (x - u) + a ^ 2 * gamma ^ 2 / 2 - tmp.m ^ 2 / 2) *
        (pnorm(tmp.d) - pnorm(-tmp.d))
    tmp.q2 <- exp(a * (x - u) + a ^ 2 * gamma ^ 2 / 2 - tmp.m^2 / 2 +
            abs(tmp.m * tmp.d) +
            log(1 + exp(-2 * abs(tmp.m * tmp.d))) - log(2)) *
        (pnorm(tmp.d) - pnorm(-tmp.d))
    q <- pmin(pmax(q, tmp.q1), tmp.q2)

    Q <- NA
    if (QFhat == TRUE){
        m <- length(a)
        Q <- rep(NA, m)
        
        II <- (1:m)[abs(a) > 1e-6]
        Q[II] <- q[II] / a[II] + (exp(a[II] * (v[II] - u[II])) * pnorm(x - v[II], mean = 0, sd = gamma) - 
            pnorm(x - u[II], mean = 0, sd = gamma)) / a[II]
    
        II <- (1:m)[abs(a) <= 1e-6]
        Q[II] <- (x - u[II]) * pnorm(x - u[II], mean = 0, sd = gamma) - (x - v[II]) * pnorm(x - v[II], mean = 0, 
            sd = gamma) + gamma ^ 2 * (dnorm(x - u[II], mean = 0, sd = gamma) - dnorm(x - v[II], mean = 0, sd = gamma))
    }

    res <- list("q" = q, "Q" = Q)
    return(res)
}



## old version, without bounds for pathological cases
# Q00 <- function(x, a, u, v, gamma, QFhat = FALSE){
# 
#     ## x, gamma in R
#     ## a, u, v in R^m
#     q <- exp(a * (x - u) + a ^ 2 * gamma ^ 2 / 2) * (pnorm(v - x - a * gamma ^ 2, mean = 0, sd = gamma) - 
#         pnorm(u - x - a * gamma ^ 2, mean = 0, sd = gamma))
# 
#     Q <- NA
#     if (QFhat == TRUE){
#         m <- length(a)
#         Q <- rep(NA, m)
#         
#         II <- (1:m)[abs(a) > 1e-6]
#         Q[II] <- q[II] / a[II] + (exp(a[II] * (v[II] - u[II])) * pnorm(x - v[II], mean = 0, sd = gamma) - 
#             pnorm(x - u[II], mean = 0, sd = gamma)) / a[II]
#     
#         II <- (1:m)[abs(a) <= 1e-6]
#         Q[II] <- (x - u[II]) * pnorm(x - u[II], mean = 0, sd = gamma) - (x - v[II]) * pnorm(x - v[II], mean = 0, sd = gamma) + gamma ^ 2 * (dnorm(x - u[II], mean = 0, sd = gamma) - dnorm(x - v[II], mean = 0, sd = gamma))
#     }
# 
#     res <- list("q" = q, "Q" = Q)
#     return(res)
# }
