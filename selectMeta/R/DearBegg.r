DearBegg <- function(y, u, lam = 2, tolerance = 10^-10, maxiter = 1000, trace = TRUE){

# lam: first weight (either 1 or 2, see Dear & Begg, p. 239. They recommend lam = 2.)

# general parameters
n <- length(y)
k <- 1 + floor(n / 2)
teststat <- abs(y) / u
p <- 2 * pnorm(-teststat)
p0 <- p

## data preparation: sort all vectors in decreasing order of p-values
ind <- order(p) 
ind <- rev(ind)
p <- p[ind]
y <- y[ind]
u <- u[ind]
teststat <- teststat[ind]

# start values
w <- 1:k / (k + 1)
theta <- 1
sigma <- 1
hij <- Hij(theta, sigma, y, u, teststat)$Hij
LL0 <- DearBeggLoglik(w, theta, sigma, y, u, hij, lam)$LL
        
LLs <- rep(NA, (k + 2) * maxiter)
LLs[1] <- LL0
iter <- 1
row <- 1
delta <- tolerance + 1

while ((iter <= maxiter) & (delta > tolerance)){

    w_old <- c(w, theta, sigma)
    LL_old <- LLs[is.na(LLs) == FALSE]
    LL_old <- LL_old[length(LL_old)]
    
    ## constrained optimization over each w_i, for fixed theta and sigma
    for (j in 1:k){
    f0 <- function(x, w, theta, sigma, y, u, lam, j){
        w[j] <- x
        hij <- Hij(theta, sigma, y, u, teststat)$Hij
        res <- DearBeggLoglik(w, theta, sigma, y, u, hij, lam)$LL
        return(res)
        }

    d0 <- optimize(f = f0, interval = c(0, 1), w, theta, sigma, y, u, lam, j, maximum = TRUE, tol = 10^-5)
    w[j] <- d0$maximum
    LLs[row] <- d0$objective
    row <- row + 1
    }

    ## constrained optimization over theta, for fixed w and sigma
    f1 <- function(x, w, sigma, y, u, lam){
        hij <- Hij(x, sigma, y, u, teststat)$Hij
        res <- DearBeggLoglik(w, x, sigma, y, u, hij, lam)$LL
        return(res)
        }
    d1 <- optimize(f = f1, interval = c(-10, 10), w, sigma, y, u, lam, maximum = TRUE, tol = 10^-5)

    theta <- d1$maximum
    LLs[row] <- d1$objective
    row <- row + 1

    ## constrained optimization over sigma, for fixed c(w, theta)
    f2 <- function(x, w, theta, y, u, lam){
        hij <- Hij(theta, x, y, u, teststat)$Hij
        res <- DearBeggLoglik(w, theta, x, y, u, hij, lam)$LL
        return(res)
        }
    d2 <- optimize(f = f2, interval = c(0, 100), w, theta, y, u, lam, maximum = TRUE, tol = 10^-5)

    sigma <- d2$maximum
    LLs[row] <- d2$objective
    row <- row + 1

    delta <- max(sqrt(sum((w_old - c(w, theta, sigma)) ^ 2)), abs(LL_old - d2$objective))
    
    if (trace == TRUE){print(paste("run: ", iter, " / LL = ", format(d2$objective, digits = 6, nsmall = 6, scientific = FALSE), " / delta = ", format(delta, digits = 1, scientific = FALSE), sep = ""))}
    iter <- iter + 1
}

res <- list("w" = w, "theta" = theta, "sigma" = sigma, "p" = p, "y" = y, "u" = u, "loglik" = max(LLs, na.rm = TRUE))
return(res)
}
