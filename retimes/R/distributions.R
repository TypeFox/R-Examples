# ---------------------------------------------
# ex-Gaussian distribution
# ---------------------------------------------

# Density:
dexgauss <- function(q, mu=0, sigma=1, tau=1) {
    D <- (1/tau) * exp(mu/tau+(sigma^2)/(2*tau^2)-q/tau) * pnorm(q,mu+(1/tau)*sigma^2,sigma)
    return(D)
}

# Random number generator:
rexgauss <- function(n, mu=0, sigma=1, tau=1, positive=TRUE) {
    if(positive) {
        while(1) {
            x <- rnorm(n,mu,sigma)+rexp(n,1/tau)
            if(sum(x>0)==n) break
        }
    } else
        x <- rnorm(n,mu,sigma)+rexp(n,1/tau)
    return(x)
}

# -------------------------------------------
# Beta-prime distribution
# -------------------------------------------

# Density:
# dbeta1 <- function(x,shape1,shape2) {
#     D <- x^(shape1-1) * (1+x)^(-shape1-shape2) / beta(shape1,shape2)
#     return(D)
# }

# Quantile:
# qbeta1 <- function(p,shape1,shape2,lower=1e-7,upper=1) {
#     Rsq.Density <- function(x,Density,shape1,shape2) {
#         Density <- get("Density",pos=qbeta1Env)
#         shape1 <- get("shape1",pos=qbeta1Env)
#         shape2 <- get("shape2",pos=qbeta1Env)
#         Residual <- ((x^(shape1-1)*(1+x)^(-shape1-shape2))/beta(shape1,shape2))-Density
#         return(abs(Residual))
#     }
#     qbeta1Env <- new.env()
#     x <- NULL
#     for(i in 1:length(p)) {
#         Density <- p[i]
#         x <- c(x,optimize(f=Rsq.Density,lower=lower,upper=upper)$minimum)
#     }
#     return(x)
# }

# # Random number generator:
# rbeta1 <- function(n,shape1,shape2) {
#     x <- rbeta(n,shape1,shape2)
#     x <- x/(1-x)
#     return(x)
# }
