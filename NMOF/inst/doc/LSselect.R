### R code from vignette source 'LSselect.Rnw'

###################################################
### code chunk number 1: LSselect.Rnw:24-25
###################################################
options(continue = " ", digits = 5)


###################################################
### code chunk number 2: LSselect.Rnw:43-44 (eval = FALSE)
###################################################
## showExample("exampleLS", chapter = "Portfolio")


###################################################
### code chunk number 3: LSselect.Rnw:47-49
###################################################
require("NMOF")
set.seed(112233)


###################################################
### code chunk number 4: LSselect.Rnw:71-77
###################################################
na <- 500L                          ## number of assets
C <- array(0.6, dim = c(na, na))    ## correlation matrix
diag(C) <- 1
minVol <- 0.20; maxVol <- 0.40      ## covariance matrix
Vols <- (maxVol - minVol) * runif(na) + minVol
Sigma <- outer(Vols, Vols) * C   


###################################################
### code chunk number 5: LSselect.Rnw:80-85
###################################################
OF <- function(x, Data) {
    w <- x/sum(x)
    res <-  crossprod(w[x], Data$Sigma[x, x])
    tcrossprod(w[x], res)
}


###################################################
### code chunk number 6: LSselect.Rnw:88-92
###################################################
OF2 <- function(x, Data) {
    w <- 1/sum(x)
    sum(w * w * Data$Sigma[x, x])
}


###################################################
### code chunk number 7: LSselect.Rnw:95-106
###################################################
neighbour <- function(xc, Data) {
    xn <- xc
    p <- sample.int(Data$na, Data$nc, replace = FALSE)
    xn[p] <- !xn[p]

    ## reject infeasible solution
    if (sum(xn) > Data$Kmax || sum(xn) < Data$Kmin)
        xc 
    else 
        xn
}


###################################################
### code chunk number 8: LSselect.Rnw:115-120
###################################################
Data <- list(Sigma = Sigma,
              Kmin = 30L,
              Kmax = 60L,
              na = na,
              nc = 1L)


###################################################
### code chunk number 9: LSselect.Rnw:126-130
###################################################
card0 <- sample(Data$Kmin:Data$Kmax, 1L, replace = FALSE)
assets <- sample.int(na, card0, replace = FALSE)
x0 <- logical(na)
x0[assets] <- TRUE


###################################################
### code chunk number 10: LSselect.Rnw:137-142
###################################################
algo <- list(x0 = x0,
              neighbour = neighbour,
              nS = 5000L,
              printDetail = FALSE,
              printBar = FALSE)


###################################################
### code chunk number 11: LSselect.Rnw:145-150
###################################################
system.time(sol1 <- LSopt(OF, algo, Data))
sqrt(sol1$OFvalue)
par(ylog = TRUE, bty = "n", las = 1, tck = 0.01, mar = c(4,4,1,1))
plot(sqrt(sol1$Fmat[ ,2L]), main = "",
     type = "l", ylab = "portfolio volatility", xlab = "iterations")


###################################################
### code chunk number 12: LSselect.Rnw:156-164
###################################################
nRuns <- 3L
allRes <- restartOpt(LSopt, n = nRuns, OF, algo = algo, Data = Data)
allResOF <- numeric(nRuns)
for (i in seq_len(nRuns))
    allResOF[i] <- sqrt(allRes[[i]]$OFvalue)
par(bty = "n", las = 1, tck = 0.01, mar = c(4,4,1,1))
plot(ecdf(allResOF), xlab = "x: Portfolio volatility", pch = 21,
     main = "")


