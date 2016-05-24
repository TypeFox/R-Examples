"iclocpoly" <-
function (x, y = NULL, y.IC, degree=0, h, niter = 10, kernel="normal", gridsize=401) 
{
    if (length(x) != (length(y) + dim(y.IC)[1])) 
        print("Error: Number of responses does not match number of predictors.")
    m <- length(y)
    n <- length(x)
    mu0 <- apply(y.IC, 1, mean)
    muhat <- c(y, mu0)
    y.locpoly <- locpoly(x, muhat, degree=degree, bandwidth=h, kernel=kernel, gridsize=gridsize)
    y.fitted <- approx(y.locpoly$x, y.locpoly$y, x)$y
    sigma <- sd(y.fitted-muhat)
    for (k in 1:niter) {
        y.standard <- (y.IC - y.fitted[(m+1):n])/sigma
        R <- y.standard[, 2]
        L <- y.standard[, 1]
        phi.Rp <- dnorm(R)
        phi.Lp <- dnorm(L)
        Phi.Rp <- pnorm(R)
        Phi.Lp <- pnorm(L)
        den <- Phi.Rp - Phi.Lp
        den[is.na(den)] <- 0
        term1 <- (phi.Rp - phi.Lp)/den
        muhat <- y.fitted - sigma * term1
        mp <- (R + L)/2
        muhat.alt <- mp
        muhat[den < (10^(-8))] <- muhat.alt[den < (10^(-8))]
        term2 <- (phi.Rp * R - phi.Lp * L)/den
        sigmahat2 <- sigma^2 * (1 - term1^2 - term2)
        muhat <- c(y, muhat)
        sigma2 <- var(term1)*sigma^2 + mean(sigmahat2)
        if (m > 0) sigma2 <- sigma2  + var((y-y.fitted[1:m]))
        sigma <- sqrt(sigma2)
        y.locpoly <- locpoly(x, muhat, degree=degree, bandwidth=h, kernel=kernel, gridsize=gridsize)
        y.fitted <- approx(y.locpoly$x, y.locpoly$y, x[(m+1):n])$y
        }
        list(x = x, y = muhat, sigma = sigma)
}

