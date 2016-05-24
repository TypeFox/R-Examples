# Gaussian2.R -- version 2011-03-18
p   <- 5L    # number of assets
N   <- 500L  # number of obs
rho <- 0.5   # correlation between two assets

## create uncorrelated observations
X <- rnorm(N * p); dim(X) <- c(N, p)

## check (see ?pairs)
panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1L:2L], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y,  ...)
}
pairs(X, xlim = c(-5, 5), ylim = c(-5,5), labels = NA, 
      diag.panel = panel.hist, col = grey(0.4))
cor(X)

## set correlation matrix
M <- array(rho, dim = c(p, p)); diag(M) <- 1

## induce correlation, check
C <- chol(M); Xc <- X %*% C
pairs(Xc, xlim = c(-5,5), ylim = c(-5,5), labels = NA, 
      diag.panel = panel.hist, col = grey(0.4))
cor(Xc)
