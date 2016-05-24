### R code from vignette source 'PSlms.Rnw'

###################################################
### code chunk number 1: PSlms.Rnw:21-22
###################################################
options(continue = " ", digits = 5)


###################################################
### code chunk number 2: PSlms.Rnw:45-48
###################################################
require("NMOF")
require("MASS")
set.seed(11223344)


###################################################
### code chunk number 3: PSlms.Rnw:52-64
###################################################
createData <- function(n, p, constant = TRUE,
                       sigma = 2, oFrac = 0.1) {
    X <- array(rnorm(n * p), dim = c(n, p))
    if (constant) 
        X[ ,1L] <- 1L
    b <- rnorm(p)
    y <- X %*% b + rnorm(n)*0.5
    nO <- ceiling(oFrac*n)
    when <- sample.int(n, nO)
    X[when, -1L] <- X[when, -1L] + rnorm(nO, sd = sigma)
    list(X = X, y = y, outliers = when)
}


###################################################
### code chunk number 4: PSlms.Rnw:78-86
###################################################
n <- 100L   ## number of observations
p <- 10L    ## number of regressors
constant <- TRUE; sigma <- 5; oFrac  <- 0.1
h <- 75L    ## ... or use something like floor((n+1)/2)

aux <- createData(n, p, constant, sigma, oFrac)
X <- aux$X; y <- aux$y
Data <- list(y = as.vector(y), X = X, h = h)


###################################################
### code chunk number 5: PSlms.Rnw:89-93
###################################################
par(bty = "n", las = 1, tck = 0.01, mar = c(4,4,1,1))
plot(X[ ,2L], type = "h", ylab = "X values", xlab = "observation")
lines(aux$outliers, X[aux$outliers ,2L], type = "p", pch = 21, 
      col = "blue", bg = "blue")


###################################################
### code chunk number 6: PSlms.Rnw:99-113
###################################################
OF <- function(param, Data) {
    X <- Data$X; y <- Data$y
    aux <- y - X %*% param
    aux <- aux * aux
    aux <- apply(aux, 2L, sort, partial = Data$h)
    colSums(aux[1:Data$h, ])  ## LTS
}
OF <- function(param, Data) {
    X <- Data$X; y <- Data$y
    aux <- y - X %*% param
    aux <- aux * aux
    aux <- apply(aux, 2L, sort, partial = Data$h)
    aux[Data$h, ]  ## LQS
}


###################################################
### code chunk number 7: PSlms.Rnw:122-160
###################################################
popsize <- 100L; generations <- 500L
ps <- list(min = rep(-10,p),
           max = rep( 10,p),
           c1 = 0.9,
           c2 = 0.9,
           iner = 0.9,
           initV = 1,
           nP = popsize,
           nG = generations,
           maxV = 5,
           loopOF = FALSE,
           printBar = FALSE,
           printDetail = FALSE)
de <- list(min = rep(-10,p),
           max = rep( 10,p),
           nP = popsize,
           nG = generations,
           F = 0.7,
           CR = 0.9,
           loopOF = FALSE,
           printBar = FALSE,
           printDetail = FALSE)

system.time(solPS <- PSopt(OF = OF, algo = ps, Data = Data))
system.time(solDE <- DEopt(OF = OF, algo = de, Data = Data))

if (require("MASS", quietly = TRUE)) {
    system.time(test1 <- lqs(y ~ X[ ,-1L], adjust = TRUE,
                             nsamp = 100000L, method = "lqs",
                             quantile = h))
    res1 <- sort((y - X %*% as.matrix(coef(test1)))^2)[h]
} else 
    res1 <- NA
res2 <- sort((y - X %*% as.matrix(solPS$xbest))^2)[h]
res3 <- sort((y - X %*% as.matrix(solDE$xbest))^2)[h]
cat("lqs:   ", res1, "\n",
    "PSopt: ", res2, "\n",
    "DEopt: ", res3, "\n", sep = "")


###################################################
### code chunk number 8: PSlms.Rnw:168-175
###################################################
popsize <- 100L; generations <- 20L
de$nP <- popsize; de$nG <- generations
ps$nP <- popsize; ps$nG <- generations

de$loopOF <- TRUE; ps$loopOF <- TRUE
t1ps <- system.time(solPS <- PSopt(OF = OF, algo = ps, Data = Data))
t1de <- system.time(solDE <- DEopt(OF = OF, algo = de, Data = Data))


###################################################
### code chunk number 9: PSlms.Rnw:179-182
###################################################
de$loopOF <- FALSE; ps$loopOF <- FALSE
t2ps <- system.time(solPS <- PSopt(OF = OF, algo = ps, Data = Data))
t2de <- system.time(solDE <- DEopt(OF = OF, algo = de, Data = Data))


###################################################
### code chunk number 10: PSlms.Rnw:185-187
###################################################
t1ps[[3L]]/t2ps[[3L]]  ## PS
t1de[[3L]]/t2de[[3L]]  ## DE


