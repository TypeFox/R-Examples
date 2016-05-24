# exampleLoop.R -- version 2010-12-30
# set up X matrix with n rows and p columns
n <- 100L; p <- 5L
X <- array(rnorm(n * p), dim = c(n, p))
# set up population
nP <- 100L 
mP <- array(rnorm(p * nP), dim = c(p, nP))

Y <- array(0, dim = c(n, nP)); Z <- Y

system.time({
        for (r in 1L:1000L) 
            for (i in 1L:nP) 
                Y[ ,i] <- X %*% mP[ ,i]
    })
system.time({
        for (r in 1L:1000L) 
            Z <- X %*% mP
    })
all.equal(Y, Z)	# ... should be TRUE