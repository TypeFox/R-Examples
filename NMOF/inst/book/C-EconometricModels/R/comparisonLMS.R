# comparisonLMS.R -- version 2011-05-17
# compare PS and DE for robust regression (LQS)
require(MASS); require(NMOF)

OF <- function(param, data) { 
    X <- data$X
    y <- data$y
    # as.vector(y) for recycling; param is a matrix
    aux <- as.vector(y) - X %*% param 
    aux <- aux * aux
    aux <- apply(aux, 2, sort, partial = data$h)
    colSums(aux[1L:data$h, ])  # LTS
}
OF <- function(param, data) { 
    X <- data$X
    y <- data$y
    # as.vector(y) for recycling; param is a matrix
    aux <- as.vector(y) - X %*% param
    aux <- aux * aux
    aux <- apply(aux, 2, sort, partial = data$h)
    aux[data$h, ]  # LQS
}

createData <- function(n, p, 
    constant = TRUE, sigma = 2, oFrac = 0.1) {
    X <- array(rnorm(n * p), dim = c(n, p))
    if (constant) X[ ,1L] <- 1
    b <- rnorm(p)
    y <- X %*% b + rnorm(n) * 0.5
    nO <- ceiling(oFrac * n)
    when <- sample.int(n, nO)
    X[when, -1L] <- X[when, -1L] + rnorm(nO, sd = sigma)
    list(X = X, y = y)
}

n <- 100  # number of observations
p <- 10   # number of regressors
constant <- TRUE; sigma <- 5; oFrac  <- 0.15
h <- 70L   # ...or use something like floor((n+1)/2)

aux <- createData(n, p, constant, sigma, oFrac) 
X <- aux$X; y <- aux$y
data <- list(y = y, X = X, h = h)

popsize <- 100L; generations <- 400L
ps <- list(
    min = rep(-10, p), 
    max = rep(10, p), 
    c1 = 1.0, 
    c2 = 2.0, 
    iner = 0.8,
    initV = 0.1, 
    maxV = 3, 
    nP = popsize, 
    nG = generations, 
    loopOF = FALSE)
de <- list(
    min = rep(-10, p), 
    max = rep(10, p),
    nP = popsize, 
    nG = generations, 
    F = 0.2, 
    CR = 0.5,
    loopOF = FALSE)

system.time(solPS <- PSopt(OF = OF,algo = ps, data = data))
system.time(solDE <- DEopt(OF = OF,algo = de, data = data))
 
# test: lqs. use X[ ,-1] because lqs includes constant.
system.time(test1 <- lqs(y ~ X[ ,-1L], adjust = TRUE,    
        nsamp = 100000L, method = "lqs", quantile = h))

res1 <- sort((y - X %*% as.matrix(coef(test1)))^2)[h]
res2 <- sort((y - X %*% as.matrix(solPS$xbest))^2)[h]
res3 <- sort((y - X %*% as.matrix(solDE$xbest))^2)[h]
cat("lqs:   ", res1, "\n",
    "PSopt: ", res2, "\n",
    "DEopt: ", res3, "\n", sep = "")