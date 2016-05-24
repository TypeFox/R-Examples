ps.normal <- function(x, y, lambda, df = 5, method = "df",
                        weights, eval.points, ngrid, weights.penalty = FALSE,
                        nseg, bdeg = 3, pord = 2, display = "lines",
                        increasing = FALSE, decreasing = FALSE, kappa = lambda * 100,
                        fixed, negative = TRUE) {

# Function ps.normal: smooths scatterplot data with P-splines.
# Input: x = abcissae of data.
# Input: y = response (counts).
# Input: nseg = number of intervals for B-splines.
# Input: bdeg = degree of B-splines.
# Input: pord = order of difference penalty.
# Input: lambda = smoothness parameter.
#
# Based on code written by Paul Eilers and Brian Marx, 2007

#     The ordering of matrices, axes etc. needs to be corrected.

    if (is.vector(x)) x <- matrix(x, ncol = 1)
    ndim <- ncol(x)
    n    <- nrow(x)
    if (missing(ngrid))  ngrid  <- switch(ndim, 100, 20, 20)
    lambda.select <- missing(lambda)
    
    if (missing(nseg))   nseg   <- switch(ndim, 100, 17, 10)

#    if (ndim == 1 & (!missing(fixed)) & (increasing | decreasing))
#       stop("monotonic increasing estimation is not available with fixed points.")
    if (ndim == 1 & increasing & decreasing)
       stop("only one of increasing and decreasing can be set to TRUE.")
    if (missing(weights)) weights <- rep(1, length(y))
    weights <- diag(weights)
    
    if (missing(eval.points)) {
       eval.points <- list(length = ndim)
       for (i in 1:ndim)
          eval.points[[i]] <- seq(min(x[,i]) - 0.05 * diff(range(x[,i])), 
                                  max(x[,i]) + 0.05 * diff(range(x[,i])),
                                  length = ngrid)
       evp <- eval.points
       eval.points <- as.matrix(expand.grid(eval.points))
    }
    else if (ndim == 1)
       eval.points <- matrix(eval.points, ncol = 1)
    
    # Compute B-spline basis
    tim <- proc.time()
    b <- list(length = ndim)
    m <- vector(length = ndim)
    for (i in 1:ndim) {
       b[[i]] <- bbase(x[,i], xl = min(x[,i]) - 0.05 * diff(range(x[,i])),
                              xr = max(x[,i]) + 0.05 * diff(range(x[,i])), 
                              nseg = nseg, deg = bdeg)
       m[i]   <- dim(b[[i]])[2]
    }
    B <- b[[1]]
    if (ndim > 1)
       B <- t(apply(cbind(b[[1]], b[[2]]), 1, 
                            function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
    if (ndim == 3)
       B <- t(apply(cbind(B,  b[[3]]), 1, 
                function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
    # cat("bases:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()

    # Construct penalty matrices
    P <- list(length = ndim)
    for (i in 1:ndim) {
       P[[i]] <- diff(diag(m[1]), diff = pord)
       if (weights.penalty)
          P[[i]] <- t(P[[i]]) %*% diag(exp(seq(1:nrow(P[[i]])) * log(100) / (nrow(P[[i]]) - 1))) %*% P[[i]]
       else
          P[[i]] <- t(P[[i]]) %*% P[[i]]
    }
    if (ndim >= 2) {
       P[[1]] <- P[[1]] %x% diag(m[2])
       P[[2]] <- diag(m[2]) %x% P[[2]]
    }
    if (ndim == 3) {
       P[[1]] <- P[[1]] %x% diag(m[3])
       P[[2]] <- P[[2]] %x% diag(m[3])
       P[[3]] <- diag(m[1]) %x% diag(m[2]) %x% P[[3]]
    }
    # cat("penalties:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    
    btb  <- t(B) %*% weights %*% B
    bty  <- t(B) %*% weights %*% y

    # Identify lambda (ndim = 1 only at the moment)
    
    if (lambda.select) {
       if (method == "df") {
          lambda.df <- function(lambda, btb, bty, P) {
             mat  <- 0
             for (i in 1:length(P))
                mat <- mat + lambda[i] * P[[i]]
             B1   <- solve(btb + mat) 
             beta <- as.vector(B1 %*% bty)
             sum(diag(btb %*% B1))
          }
          lambda <- 1
          # print(lambda.df(lambda, btb, bty, P))
          # print(adf)
          while (lambda.df(lambda, btb, bty, P) <= df) lambda <- lambda / 10
          lower <- lambda
          lambda <- 1
          while (lambda.df(lambda, btb, bty, P) >= df) lambda <- lambda * 10
          upper <- lambda
          lambda.crit <- function(lambda, btb, bty, P, df)
             lambda.df(lambda, btb, bty, P) - df
          result <- uniroot(lambda.crit, interval = c(lower, upper), btb, bty, P, df)
          lambda <- result$root
       }
    }
    
    # Fit
    
    btb   <- t(B) %*% weights %*% B
    bty   <- t(B) %*% weights %*% y
    # cat("matrices:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    mat <- 0
    # print(lambda)
    for (i in 1:length(P))
       mat <- mat + lambda[i] * P[[i]]
    B1  <- solve(btb + mat) 
    # cat("inversion:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    beta  <- as.vector(B1 %*% bty)
    # cat("solution:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    # f     <- lsfit(rbind(B, P1, P2, P3), c(y, nix), intercept = FALSE)
    # h     <- hat(f$qr)[1:m]
    # beta  <- f$coef

    mu  <- c(B %*% beta)
    edf <- sum(diag(btb %*% B1))
    # cat("fitting:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    
    # Adjust for monotonicity and for fixed points - one covariate only
    
    if (ndim == 1 & !missing(fixed)) {
       if (any(fixed[,1] <= min(x) - 0.05 * diff(range(x[,i]))) | 
           any(fixed[,1] >= max(x) + 0.05 * diff(range(x[,i]))))
          stop("fixed points must be inside the range of the data.")
       fixed <- matrix(c(fixed), ncol = ndim + 1)
       A     <- bbase(fixed[ , 1:ndim], xl = min(x[,i]) - 0.05 * diff(range(x[,i])),
                              xr = max(x[,i]) + 0.05 * diff(range(x[,i])), 
                              nseg = nseg, deg = bdeg)
       beta  <- beta + 
                 B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% (fixed[ , ndim + 1] - A %*% beta)
       edf   <- NA
    }
    
    if (ndim == 1 & (increasing | decreasing | !negative)) {
       if (!missing(fixed))
          stop("fixed values cannot be used with increasing/decreasing/non-negative constraints")
       D1    <- diff(diag(m[1]), diff = 1)
       delta <- 1
       while (delta > 1e-5) {
          mat1 <- mat
       	  if (increasing | decreasing) {
             if (increasing) v <- as.numeric(diff(beta) <= 0)
             if (decreasing) v <- as.numeric(diff(beta) >= 0)
             mat1 <- mat1 + kappa * t(D1) %*% diag(v) %*% D1
       	  }
          if (!negative) {
             # A    <- bbase(seq(0, 1, length = 20), xl = 0, xr = 1, nseg = nseg, deg = bdeg)
             # v    <- as.numeric(c(A %*% beta) < 0)
             # mat1 <- mat1 + kappa * (t(A) %*% diag(v) %*% A)
             v    <- as.numeric(beta < 0)
             mat1 <- mat1 + kappa * diag(v)
          }
          B1       <- solve(t(B) %*% weights %*% B + mat1)
          beta.old <- beta
          beta     <- as.vector(B1 %*% t(B) %*% weights %*% y)
          delta    <- sum((beta - beta.old)^2) / sum(beta.old^2)
       }
    }
   
    # Cross-validation and dispersion
    # r     <- (y - mu ) / (1 - h)
    # cv    <- sqrt(sum(r ^2))
    # sigma <- sqrt(sum((y - mu) ^2) / (n - sum(h)))

    # Evaluate the estimate at eval.points
    for (i in 1:ndim) {
       b[[i]] <- bbase(eval.points[,i], xl = min(x[,i]) - 0.05 * diff(range(x[,i])),
                                        xr = max(x[,i]) + 0.05 * diff(range(x[,i])), 
                                        nseg = nseg, deg = bdeg)
       m[i]   <- dim(b[[i]])[2]
    }
    B <- b[[1]]
    if (ndim > 1)
       B <- t(apply(cbind(b[[1]], b[[2]]), 1, 
                            function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
    if (ndim == 3)
       B <- t(apply(cbind(B,  b[[3]]), 1, 
                function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
    est <- c(B %*% beta)
    if (ndim > 1) est <- array(est, dim = rep(ngrid, ndim))
    # cat("estimate:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()

    # Plot data and fit
    if (display != "none") {
       if (ndim == 1) {
          plot(x, y, main = '', xlab = '', ylab = '')
          lines(eval.points[,1], est, col = 'blue')
#          if (se > 0 ) {
#             Covb = solve(btb + P[[1]])
#             Covz = sigma^2 * B %*% Covb %*% t(B)
#             seb = se * sqrt(diag(Covz))
#             lines(u, est + seb, lty = 2, col = 'red')
#             lines(u, est - seb, lty = 2, col = 'red')
#          }
       }
       else if (ndim == 2) {
          persp(evp[[1]], evp[[2]], est, 
                  ticktype = "detailed", col = "green", d = 10, theta = 30)
       }
    }

    # Return list
    pp <- list(x = x, y = y, muhat = mu, nseg = nseg,
               bdeg = bdeg, pord = pord, beta = beta, B = B, b = b,
               lambda = lambda, eval.points = eval.points, estimate = est,
               btb = btb, P = P, bty = bty, df = edf, eval.points = eval.points,
               # cv = cv, effdim = sum(h), ed.resid = m - sum(h), sigma = sigma,
               family = "gaussian", link = "identity")
    class(pp) = "pspfit"
    return(invisible(pp))
}

#     Paul Eilers material

tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)


bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3){
# Construct B-spline basis
    dx <- (xr - xl) / nseg
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B }

gauss <- function(x, mu, sig) {
# Gaussian-shaped function
    u <- (x - mu) / sig
    y <- exp(- u * u / 2)
    y }

gbase <- function(x, mus) {
# Construct Gaussian basis
    sig <- (mus[2] - mus[1]) / 2
    G <- outer(x, mus, gauss, sig)
    G }

pbase <- function(x, n) {
# Construct polynomial basis
    u <- (x - min(x)) / (max(x) - min(x))
    u <- 2 * (u - 0.5);
    P <- outer(u, seq(0, n, by = 1), "^")
    P }
