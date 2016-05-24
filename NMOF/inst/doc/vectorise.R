### R code from vignette source 'vectorise.Rnw'

###################################################
### code chunk number 1: vectorise.Rnw:22-23
###################################################
options(continue = " ", digits = 5)


###################################################
### code chunk number 2: vectorise.Rnw:85-86
###################################################
require("NMOF")


###################################################
### code chunk number 3: vectorise.Rnw:96-97
###################################################
tfRosenbrock


###################################################
### code chunk number 4: vectorise.Rnw:100-104
###################################################
OF <- tfRosenbrock     ## see ?testFunctions
size <- 5L             ## set dimension
x <- rep.int(1, size)  ## the known solution ...
OF(x)                  ## ... should give zero


###################################################
### code chunk number 5: vectorise.Rnw:108-115
###################################################
algo <- list(printBar = FALSE,
                   nP = 50L,
                   nG = 500L,
                    F = 0.6,
                   CR = 0.9,
                  min = rep(-100, size),
                  max = rep( 100, size))


###################################################
### code chunk number 6: vectorise.Rnw:118-124
###################################################
## a vectorised OF: works only with *matrix* x
OF2 <- function(x) {
    n <- dim(x)[1L]
    xi <- x[1L:(n - 1L), ]
    colSums(100 * (x[2L:n, ] - xi * xi)^2 + (1 - xi)^2)
}


###################################################
### code chunk number 7: vectorise.Rnw:127-131
###################################################
x <- matrix(rnorm(size * algo$nP), size, algo$nP)
c(OF(x[ ,1L]), OF(x[ ,2L]), OF(x[ ,3L]))
OF2(x)[1L:3L]  ## should give the same result
all.equal(OF2(x)[1L:3L], c(OF(x[ ,1L]), OF(x[ ,2L]), OF(x[ ,3L])))


###################################################
### code chunk number 8: vectorise.Rnw:137-143
###################################################
set.seed(1223445)
(t1 <- system.time(sol <- DEopt(OF = OF, algo = algo)))

algo$loopOF <- FALSE
set.seed(1223445)
(t2 <- system.time(sol2 <- DEopt(OF = OF2, algo = algo)))


###################################################
### code chunk number 9: vectorise.Rnw:146-149
###################################################
sol$OFvalue    ## both should be zero (with luck)
sol2$OFvalue
t1[[3L]]/t2[[3L]]  ## speedup


###################################################
### code chunk number 10: vectorise.Rnw:173-184
###################################################
na <- 100L  ## number of assets
np <- 100L  ## size of population
trials <- seq_len(100L)  ## for speed test

## a covariance matrix
Sigma <- array(0.7, dim = c(na, na)); diag(Sigma) <- 1

## set up population
W <- array(runif(na * np), dim = c(na, np))
## budget constraint
scaleFun <- function(x) x/sum(x); W <- apply(W, 2L, scaleFun)


###################################################
### code chunk number 11: vectorise.Rnw:187-207
###################################################
## variant 1
t1 <- system.time({
for (i in trials) {
    res1 <- numeric(np)
    for (j in seq_len(np)) {
        w <- W[ ,j]
        res1[j] <- w %*% Sigma %*% w
    }
}
})

## variant 2
t2 <- system.time({
    for (i in trials) res2 <- diag(t(W) %*% Sigma %*% W)
})

## variant 3
t3 <- system.time({
    for (i in trials) res3 <- colSums(Sigma %*% W * W)
})


###################################################
### code chunk number 12: vectorise.Rnw:210-212
###################################################
all.equal(res1,res2)
all.equal(res2,res3)


###################################################
### code chunk number 13: vectorise.Rnw:215-222
###################################################
## time required
#  ... variant 1
t1
## ... variant 2
t2
## ... variant 3
t3


###################################################
### code chunk number 14: vectorise.Rnw:243-257
###################################################
n  <- 100L  # number of observation
p  <- 5L    # number of regressors
np <- 100L  # population size
trials <- seq_len(1000L)

## random data
X <- array(rnorm(n * p), dim = c(n, p))
y <- rnorm(n)

## random population
Theta <- array(rnorm(p * np), dim = c(p, np))

## empty residuals matrix
R1 <- array(NA, dim = c(n, np))


###################################################
### code chunk number 15: vectorise.Rnw:260-269
###################################################
system.time({
    for (i in trials)
        for (j in seq_len(np))
            R1[ ,j] <- y - X %*% Theta[ ,j]
})
system.time({
    for (i in trials)
        R2 <- y - X %*% Theta
})


###################################################
### code chunk number 16: vectorise.Rnw:274-275
###################################################
all.equal(R1, R2)  ## ... should be TRUE


