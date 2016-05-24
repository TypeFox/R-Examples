pureCMAES <- function(par, fun, lower = NULL, upper = NULL, sigma = 0.5,
                      stopfitness = -Inf, stopeval = 1000*length(par)^2, ...)
{
    stopifnot(is.numeric(par))
    if (sigma < 0.1 || sigma > 0.9) sigma <- max(min(sigma, 0.9), 0.1)
    N <- length(par)            # number of objective variables/problem dimension    
    if (length(lower) == 1) lower <- rep(lower, N)
    if (length(upper) == 1) upper <- rep(upper, N)
    fct <- match.fun(fun)
    fun <- function(x) fct(x, ...)

    xmean <- par                # objective variables initial point
                                # coordinate wise standard deviation (step size)
    sigma <- sigma * (upper - lower)
    # stopfitness <- 1e-10      # stop if fitness < stopfitness (minimization)
    # stopeval <- 1e3*N^2       # stop after stopeval number of function evaluations
    
    # Strategy parameter setting: Selection  
    lambda <- 4+floor(3*log(N))             # population size, offspring number
    mu <- lambda/2                          # number of parents/points for recombination
    weights <- log(mu+1/2)-log(1:mu)        # muXone array for weighted recombination
    mu <- floor(mu)         
    weights <- weights/sum(weights)         # normalize recombination weights array
    mueff <- sum(weights)^2/sum(weights^2)  # variance-effectiveness of sum w_i x_i
    
    # Strategy parameter setting: Adaptation
    cc <- (4+mueff/N) / (N+4 + 2*mueff/N)  # time constant for cumulation for C
    cs <- (mueff+2) / (N+mueff+5)          # t-const for cumulation for sigma control
    c1 <- 2 / ((N+1.3)^2+mueff)            # learning rate for rank-one update of C
    cmu <- min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff))   # and for rank-mu update
    damps <- 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs  # damping for sigma 
                                                         # usually close to 1
    # Initialize dynamic (internal) strategy parameters and constants
    pc <- ps <- numeric(N)               # evolution paths for C and sigma
    B <- diag(N)                        # B defines the coordinate system
    D <- rep(1, N)                       # diagonal D defines the scaling
    C <- as.matrix(B) %*% diag(D^2) %*% B             #'covariance matrix C
    invsqrtC <- as.matrix(B) %*% diag(D^-1) %*% B     #C^-1/2 
    eigeneval <- 0                       # track update of B and D
    chiN <- N^0.5 * (1-1/(4*N) + 1/(21*N^2)) # expectation of 
                                         #   ||N(0,I)|| == norm(randn(N,1))

    ml.triu <- function (M, k = 0) {
        if (k == 0) M[lower.tri(M, diag = FALSE)] <- 0
        else        M[col(M) <= row(M) + k - 1] <- 0
        return(M)
    }

    counteval <- 0
    while (counteval < stopeval)
    {
        # Generate and evaluate lambda offspring
        arx <- matrix(0, nrow = N, ncol = lambda)
        arfitness <- numeric(lambda)
        for (k in 1:lambda) {
            # arx[, k] <- xmean + sigma * B %*% (D * rnorm(N)) # m + sig * Normal(0,C)
        arxk <- xmean + sigma * B %*% (D * rnorm(N))
        arxk <- ifelse(arxk > lower, ifelse(arxk < upper, arxk, upper), lower)
        arx[, k] <- arxk
            arfitness[k] <- fun(arx[, k]) # objective function call
            counteval <- counteval + 1
        }

        # Sort by fitness and compute weighted mean
        arindex <- order(arfitness, decreasing = FALSE)
        arfitness <- arfitness[arindex]
        xold <- xmean
        xmean <- arx[, arindex[1:mu]] %*% weights

        # Cumulation: Update evolution paths
        ps   <- (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * invsqrtC %*% (xmean-xold) / sigma
        hsig <- norm(ps, "F")/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1)
        pc   <- (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma

        # Adapt covariance matrix C
        artmp <- (1/sigma) * (arx[, arindex[1:mu]] - matrix(1, 1, mu) %x% xold)
        C <- (1-c1-cmu) * C +                                 # regard old matrix  
                 c1 * (pc %*% t(pc) +                         # plus rank one update
                 (1-hsig) * cc*(2-cc) * C) +                  # minor correction if hsig==0
                 cmu * artmp %*% diag(weights) %*% t(artmp)   # plus rank mu update

        # Adapt step size sigma
        sigma <- sigma * exp((cs/damps)*(norm(ps, "F")/chiN - 1))

        # Decomposition of C into B*diag(D.^2)*B' (diagonalization)
        if (counteval - eigeneval > lambda/(c1+cmu)/N/10) { # to achieve O(N^2)
            eigeneval <- counteval
            C <- ml.triu(C) + t(ml.triu(C,1))       # enforce symmetry
            EigVal <- eigen(C, symmetric = TRUE)    # eigen decomposition
            B <- EigVal$vectors
            D <- sqrt(EigVal$values)    # D is a vector of standard deviations now
            invsqrtC <- B %*% diag(D^-1) %*% t(B)
        }

        # Break, if fitness is good enough or 
        if (arfitness[1] <= stopfitness || max(D) > 1e7 * min(D))
            break
    }

    xmin <- arx[, arindex[1]]   # Return best point of last iteration
    return(list(xmin = xmin, fmin = fun(xmin)))
}
