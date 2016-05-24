### R code from vignette source 'NMOFman.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: NMOFman.Rnw:145-147
###################################################
version <- as.Date("2014-01-27")
options(continue = " ", digits = 3, width = 60)


###################################################
### code chunk number 2: NMOFman.Rnw:229-231 (eval = FALSE)
###################################################
## code <- system.file("NMOFex/NMOFman.R", package = "NMOF")
## file.show(code, title = "NMOF manual")


###################################################
### code chunk number 3: NMOFman.Rnw:291-293
###################################################
require("NMOF")
set.seed(123321)


###################################################
### code chunk number 4: NMOFman.Rnw:319-327
###################################################
n <- 100L         ## set number of rows
rho <- 0.7        ## set a correlation between a and b

C <- array(rho,   ## create correlation matrix
            dim = c(2L, 2L))
diag(C) <- 1
ab <- array(rnorm(n * 2), dim = c(n, 2L)) %*% chol(C)
colnames(ab) <- c("a", "b")


###################################################
### code chunk number 5: NMOFman.Rnw:330-340
###################################################
par(mfrow = c(1L, 3L), bty = "n", mar = c(4, 4, 0, 0),
    tck = 0.01, las = 1, ps = 9, mgp = c(1.8, 0.5, 0),
    cex.axis = 1.2, cex.lab = 1.2, lab = c(3,3,7))
prx <- pretty(ab)
lims <- c(min(prx), max(prx))
plot(ab, pch = 19, cex = 0.5, asp = 1,xlim = lims, ylim = lims,
     xaxt = "n", yaxt = "n", xlab = "a", ylab = "")
axis(1, col = grey(0.5))
axis(2, col = grey(0.5))
mtext("b", side = 2, line = 1.8, cex = 0.8)


###################################################
### code chunk number 6: NMOFman.Rnw:355-396
###################################################
plotSubsets <- function(ab, subset1, subset2, linesAtZero = FALSE) {
    if (missing(subset2))
        subset2 <- !subset1
    if (cor(ab[subset1, ])[1L,2L] < cor(ab[subset2, ])[1L,2L]) {
        tmp <- subset1
        subset1 <- subset2
        subset2 <- tmp
    }

    par(mfrow = c(1L, 3L), bty = "n", mar = c(4, 4, 0, 0),
        tck = 0.01, las = 1, ps = 9, mgp = c(1.8,0.5,0), 
        cex.axis = 1.2, cex.lab = 1.2, lab = c(3,3,7))

    prx <- pretty(ab)
    lims <- c(min(prx), max(prx))

    plot(ab[subset1, ], xlim = lims, ylim = lims, xlab = "a", ylab = "",
         col = grey(0.2), pch = 19, cex = 0.5, xaxt = "n", yaxt = "n", asp = 1)
    axis(1, col = grey(0.6))
    lines(ab[subset2, ], col = grey(0.75), pch = 19, type = "p", cex = 0.5, asp = 1)
    axis(2, col = grey(0.6))
    mtext("b", side = 2, line = 1.8, cex = 0.8)
    if (linesAtZero)
        abline(v = 0, h = 0, col=grey(0.6))

    par(mar = c(4,2,0,2))
    plot(ab[subset1, ], xlim=lims, ylim = lims, xlab = "a", ylab = "",
         col = grey(0.2), pch = 19, cex = 0.5, xaxt = "n", yaxt = "n", asp = 1)
    axis(1, col = grey(0.6))
    if (linesAtZero)
        abline(v=0,h=0, col=grey(0.5))

    par(mar = c(4,0,0,4))
    plot(ab[subset2, ], xlim=lims, ylim = lims, xlab = "a", ylab = "",
         col = grey(0.75), pch = 19, cex = 0.5, xaxt = "n", yaxt = "n", asp = 1)
    axis(1, col = grey(0.6))
    if (linesAtZero)
        abline(v=0,h=0, col=grey(0.5))
    invisible(NULL)
}
plotSubsets(ab, 1:50, 51:100)


###################################################
### code chunk number 7: NMOFman.Rnw:402-404
###################################################
cor(ab[ 1: 50, ])
cor(ab[51:100, ])


###################################################
### code chunk number 8: NMOFman.Rnw:424-425
###################################################
x0 <- rep(c(TRUE, FALSE), each = 50L)


###################################################
### code chunk number 9: NMOFman.Rnw:429-430
###################################################
x0 <- runif(n) > 0.5


###################################################
### code chunk number 10: NMOFman.Rnw:433-434 (eval = FALSE)
###################################################
## ab[x, ]


###################################################
### code chunk number 11: NMOFman.Rnw:437-438 (eval = FALSE)
###################################################
## ab[!x, ]


###################################################
### code chunk number 12: NMOFman.Rnw:457-459
###################################################
OF <- function(x, ab)
    -abs(cor(ab[x, ])[1L, 2L] - cor(ab[!x, ])[1L, 2L])


###################################################
### code chunk number 13: NMOFman.Rnw:463-466
###################################################
x0 <- rep(c(TRUE, FALSE), each = 50L)
OF( x0, ab)
OF(!x0, ab) ## should give the same result


###################################################
### code chunk number 14: NMOFman.Rnw:486-498
###################################################
trials <- 1e5
OFvalues <- numeric(trials)
solutions <- vector("list", trials)
for (i in seq_len(trials)) {

    c1 <- sample(21:80, 1L)    ## cardinality of subset 1
    x0 <- logical(n)
    x0[sample.int(n, c1)] <- TRUE

    OFvalues[i] <- OF(x0, ab)  ## store results
    solutions[[i]] <- x0
}


###################################################
### code chunk number 15: NMOFman.Rnw:502-503
###################################################
summary(OFvalues)


###################################################
### code chunk number 16: NMOFman.Rnw:506-508
###################################################
xbest <- which.min(OFvalues)
OFvalues[xbest]


###################################################
### code chunk number 17: NMOFman.Rnw:515-516
###################################################
xRandom <- solutions[[xbest]]


###################################################
### code chunk number 18: NMOFman.Rnw:534-536
###################################################
subset1 <- ab[ ,1L] * ab[ ,2L] >  0
subset2 <- ab[ ,1L] * ab[ ,2L] <= 0


###################################################
### code chunk number 19: NMOFman.Rnw:541-542
###################################################
plotSubsets(ab, subset1, linesAtZero=TRUE)


###################################################
### code chunk number 20: NMOFman.Rnw:555-556
###################################################
OF(subset1, ab)


###################################################
### code chunk number 21: NMOFman.Rnw:561-563
###################################################
sum(subset1)
sum(subset2)


###################################################
### code chunk number 22: NMOFman.Rnw:570-577
###################################################
cr <- order(ab[ ,1L] * ab[ ,2L])
OFvalues <- numeric(n)
for (i in 20:80) {
    x0 <- logical(n)
    x0[cr[seq_len(i)]] <- TRUE
    OFvalues[i] <- OF(x0, ab)
}


###################################################
### code chunk number 23: NMOFman.Rnw:580-585
###################################################
cutoff <- which.min(OFvalues)
subset1 <- logical(n)
subset1[cr[seq_len(n) <= cutoff]] <- TRUE
subset2 <- !subset1
OF(subset1, ab)


###################################################
### code chunk number 24: NMOFman.Rnw:590-591
###################################################
xConstr <- subset1


###################################################
### code chunk number 25: NMOFman.Rnw:628-663
###################################################
greedy <- function(fun, x0, ab, n, nmin, maxit = 1000L) {
    done <- FALSE
    xbest <- xc <- x0
    xbestF <- xcF <- fun(xbest, ab)
    ic <- 0

    while (!done) {
        if (ic > maxit)
            break
        else
            ic <- ic + 1L

        done <- TRUE
        xc <- xbest
        for (i in seq_len(n)) {

            xn <- xc            ## create a new solution
            xn[i] <- !xn[i]

            sxn <- sum(xn)      ## check constraints 
            enough <- sxn >= nmin
            notTooMany <- sxn <= n - nmin

            if (enough && notTooMany) {
                xnF <- fun(xn, ab)
                if (xnF < xbestF) {
                    xbest <- xn
                    xbestF <- xnF
                    done <- FALSE
                }
            }
        }
    }
    list(xbest = xbest, OFvalue = xbestF, ic = ic)
}


###################################################
### code chunk number 26: NMOFman.Rnw:670-677
###################################################
x0 <- rep(c(TRUE, FALSE), each = 50L)
result <- greedy(fun = OF, x0 = x0, ab = ab,
                 n = 100, nmin = 20, maxit = 1000L)
xGreedy <- result$xbest
OF(x0, ab)
OF(xGreedy, ab)
result$OFvalue


###################################################
### code chunk number 27: NMOFman.Rnw:680-681
###################################################
result$ic


###################################################
### code chunk number 28: NMOFman.Rnw:687-702
###################################################
trials <- 1000
OFvalues <- numeric(trials)
solutions <- vector("list", trials)
moves <- numeric(trials)
for (i in seq_len(trials)) {

    ## create a random solution and improve it
    x0 <- runif(n) > 0.5
    result <- greedy(fun = OF, x0 = x0, ab = ab,
                     n = 100, nmin = 20, maxit = 1000L)

    OFvalues[i] <- result$OFvalue
    solutions[[i]] <- result$xbest
    moves[i] <- result$ic
}


###################################################
### code chunk number 29: NMOFman.Rnw:706-707
###################################################
summary(OFvalues)


###################################################
### code chunk number 30: NMOFman.Rnw:714-715
###################################################
summary(moves)


###################################################
### code chunk number 31: NMOFman.Rnw:723-726
###################################################
xbest <- which.min(OFvalues)
OFvalues[xbest]
xGreedy <- solutions[[xbest]]


###################################################
### code chunk number 32: NMOFman.Rnw:750-751
###################################################
Data <- list(ab = ab, n = n, nmin = 20L)


###################################################
### code chunk number 33: NMOFman.Rnw:756-757
###################################################
x0 <- runif(n) > 0.5


###################################################
### code chunk number 34: NMOFman.Rnw:764-779
###################################################
neighbour <- function(xc, Data) {
    xn <- xc
    
    p <- sample.int(Data$n, size = 1L) 
    xn[p] <- !xn[p]

    sxn <- sum(xn)
    enough <- sxn >= Data$nmin
    notTooMany <- sxn <= (Data$n - Data$nmin)

    if (enough && notTooMany)
        xn
    else
        xc
}


###################################################
### code chunk number 35: NMOFman.Rnw:788-789
###################################################
table(x0 == neighbour(x0, Data))


###################################################
### code chunk number 36: NMOFman.Rnw:792-794
###################################################
OF <- function(x, Data)
    -abs(cor(Data$ab[x, ])[1L, 2L] - cor(Data$ab[!x, ])[1L, 2L])


###################################################
### code chunk number 37: NMOFman.Rnw:797-799
###################################################
OF(x0, Data)
OF(neighbour(x0, Data), Data)


###################################################
### code chunk number 38: NMOFman.Rnw:803-809
###################################################
algo <- list(nS = 3000L,             ## number of steps to make
              neighbour = neighbour,  ## neighbourhood function
              x0 = x0,                ## initial solution
              printBar = FALSE)       
sol1 <- LSopt(OF, algo = algo, Data = Data)
(sol1$OFvalue)


###################################################
### code chunk number 39: NMOFman.Rnw:814-818
###################################################
subset1 <-  sol1$xbest
subset2 <- !sol1$xbest
c1 <- cor(ab[subset1, ])[1L, 2L]
c2 <- cor(ab[subset2, ])[1L, 2L]


###################################################
### code chunk number 40: NMOFman.Rnw:822-823
###################################################
plotSubsets(ab, subset1)


###################################################
### code chunk number 41: ranfun
###################################################
makex0 <- function() {
    x <- logical(100)
    while (sum(x) > 80 || sum(x) < 20)
        x <- runif(100) > 0.5
    x
}


###################################################
### code chunk number 42: NMOFman.Rnw:848-859
###################################################
algo <- list(nS = 4000L,            ## number of steps to make
             neighbour = neighbour,  ## neighbourhood function
             x0 = makex0,            ## initial solution
             printBar = FALSE,
             printDetail = FALSE)
restarts1 <- restartOpt(LSopt, 100, OF = OF, algo = algo, Data)
restarts1OFvalues <- sapply(restarts1, `[[`, "OFvalue")

algo$nS <- 8000L
restarts2 <- restartOpt(LSopt, 100, OF = OF, algo = algo, Data)
restarts2OFvalues <- sapply(restarts2, `[[`, "OFvalue")


###################################################
### code chunk number 43: NMOFman.Rnw:864-870
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0), ps = 8, tck = 0.001)
plot( ecdf(restarts1OFvalues), main = "", ylab = "", xlab = "",
     cex = 0.4, pch = 19, col = grey(.2), xlim = c(-2,-1))
lines(ecdf(restarts2OFvalues),
     cex = 0.4, pch = 19, col = grey(.6))
abline(v = OF(xConstr, Data))


###################################################
### code chunk number 44: NMOFman.Rnw:890-895
###################################################
algo$q <- 0.9
algo$nT <- 10
algo$nS <- 400
sol2 <- TAopt(OF, algo = algo, Data = Data)
sol2$OFvalue


###################################################
### code chunk number 45: NMOFman.Rnw:899-904
###################################################
subset1 <-  sol2$xbest
subset2 <- !sol2$xbest
c1 <- cor(ab[subset1, ])[1L, 2L]
c2 <- cor(ab[subset2, ])[1L, 2L]
OF(sol2$xbest, Data)


###################################################
### code chunk number 46: NMOFman.Rnw:908-909
###################################################
plotSubsets(ab, subset1)


###################################################
### code chunk number 47: NMOFman.Rnw:915-923
###################################################
algo$printBar <- FALSE
algo$printDetail <- FALSE
restarts3 <- restartOpt(TAopt, 100, OF = OF, algo = algo, Data)
restarts3OFvalues <- sapply(restarts3, `[[`, "OFvalue")

algo$nS <- 1500
restarts4 <- restartOpt(TAopt, 100, OF = OF, algo = algo, Data)
restarts4OFvalues <- sapply(restarts4, `[[`, "OFvalue")


###################################################
### code chunk number 48: NMOFman.Rnw:928-938
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0), ps = 8, tck = 0.001)
plot( ecdf(restarts1OFvalues), main = "", ylab = "", xlab = "",
     cex = 0.4, pch = 19, col = grey(.2), xlim = c(-2,-1))
lines(ecdf(restarts2OFvalues),
     cex = 0.4, pch = 19, col = grey(.6))
lines(ecdf(restarts3OFvalues),
     cex = 0.4, pch = 19, col = grey(.2), lty=2)
lines(ecdf(restarts4OFvalues),
     cex = 0.4, pch = 19, col = grey(.6), lty=2)
abline(v = OF(xConstr, Data))


###################################################
### code chunk number 49: NMOFman.Rnw:956-971
###################################################
neighbour <- function(xc, Data) {
    xn <- xc
    
    p <- sample.int(Data$n, size = Data$size) 
    xn[p] <- !xn[p]

    sxn <- sum(xn)
    enough <- sxn >= Data$nmin
    notTooMany <- sxn <= (Data$n - Data$nmin)

    if (enough && notTooMany)
        xn
    else
        xc
}


###################################################
### code chunk number 50: NMOFman.Rnw:978-998
###################################################
neighbour.maker <- function(n, nmin, size) {
    force(n)
    force(nmin)
    
    function(xc) {
        xn <- xc
    
        p <- sample.int(n, size = size) 
        xn[p] <- !xn[p]
        
        sxn <- sum(xn)
        enough <- sxn >= nmin
        notTooMany <- sxn <= (n - nmin)
        
        if (enough && notTooMany)
            xn
        else
            xc
    }
}


###################################################
### code chunk number 51: NMOFman.Rnw:1002-1004
###################################################
N <- neighbour.maker(n = 10, nmin = 2, size = 3)
N


###################################################
### code chunk number 52: NMOFman.Rnw:1006-1015
###################################################
compareLogicals <- function(x, y, ...) {
    argsL <- list(...)
    if (!("sep" %in% names(argsL))) 
        argsL$sep <- ""
    do.call("cat",
            c(list("\n",as.integer(x), "\n", as.integer(y), "\n",
                   ifelse(x == y, " ", "^"), "\n"), argsL))
    message("The vectors differ in ", sum(x != y), " place(s).")
}


###################################################
### code chunk number 53: NMOFman.Rnw:1019-1022
###################################################
x0 <- rep(c(TRUE, FALSE), each = 5L)
Data <- list(n = 10, nmin = 2, size = 1)
compareLogicals(x0, neighbour(x0, Data))


###################################################
### code chunk number 54: NMOFman.Rnw:1026-1030
###################################################
N <- neighbour.maker(n = 10, nmin = 2, size = 1)
compareLogicals(x0, N(x0))
compareLogicals(x0, N(x0))
compareLogicals(x0, N(x0))


###################################################
### code chunk number 55: NMOFman.Rnw:1034-1038
###################################################
N <- neighbour.maker(n = 10, nmin = 2, size = 3)
compareLogicals(x0, N(x0))
compareLogicals(x0, N(x0))
compareLogicals(x0, N(x0))


###################################################
### code chunk number 56: NMOFman.Rnw:1051-1054
###################################################
require("NMOF")
require("rbenchmark")
set.seed(46457)


###################################################
### code chunk number 57: NMOFman.Rnw:1080-1089
###################################################
randomData <- function(p, n, rscale = 0.5) {
    X <- array(rnorm(n * p), dim = c(n, p))
    k <- sample.int(p, 1L)    ## the number of regressors
    K <- sample.int(p, k)     ## the set of regressors
    betatrue <- numeric(p)
    betatrue[K] <- rnorm(k)   ## the true coefficients
    y <- X %*% betatrue + rnorm(n)*rscale
    list(X = X, y = y, betatrue = betatrue, K = K, n = n, p = p)
}


###################################################
### code chunk number 58: NMOFman.Rnw:1097-1100
###################################################
n <- 60L
p <- 5L
rD <- randomData(p, n)


###################################################
### code chunk number 59: NMOFman.Rnw:1112-1113
###################################################
b0 <- rnorm(p)


###################################################
### code chunk number 60: NMOFman.Rnw:1124-1128
###################################################
Data <- list(X = rD$X,
             y = rD$y,
             p = rD$p,
             n = rD$n)


###################################################
### code chunk number 61: NMOFman.Rnw:1134-1138
###################################################
OFls <- function(b, Data) {
    tmp <- Data$y - Data$X %*% b
    sum(tmp * tmp)
}


###################################################
### code chunk number 62: NMOFman.Rnw:1144-1149
###################################################
tmp <- rnorm(1e4)
benchmark(ignore1 <- tmp * tmp * tmp,
          ignore2 <- tmp^3,
          columns = c("test", "elapsed", "relative"),
          replications = 1e3, order = "relative")


###################################################
### code chunk number 63: NMOFman.Rnw:1154-1155
###################################################
all.equal(ignore1, ignore2)


###################################################
### code chunk number 64: NMOFman.Rnw:1160-1166
###################################################
algo <- list(nG = 200,  ## number of generations
             nP = 50,   ## population size
             min = rep(-20, p),
             max = rep( 20, p),
             printBar = FALSE)
resDE <- DEopt(OFls, algo = algo, Data = Data)


###################################################
### code chunk number 65: NMOFman.Rnw:1172-1174
###################################################
data.frame(QR = qr.solve(Data$X, Data$y),
           DE = resDE$xbest)


###################################################
### code chunk number 66: NMOFman.Rnw:1189-1192
###################################################
algo <- list(nP = 50,
             min = rep(-20, p), max = rep( 20, p),
             printBar = FALSE, printDetail = FALSE)


###################################################
### code chunk number 67: NMOFman.Rnw:1198-1206
###################################################
algo$nG <- 25
results1 <- restartOpt(DEopt, n = 50, OF= OFls, algo = algo, Data = Data)
algo$nG <- 50
results2 <- restartOpt(DEopt, n = 50, OF= OFls, algo = algo, Data = Data)
algo$nG <- 100
results3 <- restartOpt(DEopt, n = 50, OF= OFls, algo = algo, Data = Data)
algo$nG <- 200
results4 <- restartOpt(DEopt, n = 50, OF= OFls, algo = algo, Data = Data)


###################################################
### code chunk number 68: NMOFman.Rnw:1212-1217
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0), ps = 8, tck = 0.001)
plot( ecdf(sapply(results1, `[[`, "OFvalue")), main = "", ylab = "", xlab = "",
     cex = 0.4, pch = 19, col = grey(.6), xlim = c(0,150), ylim = c(0,1))
lines(ecdf(sapply(results2, `[[`, "OFvalue")),
     cex = 0.4, pch = 19, col = grey(.4))


###################################################
### code chunk number 69: NMOFman.Rnw:1222-1227
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0), ps = 8, tck = 0.001)
plot( ecdf(sapply(results2, `[[`, "OFvalue")), main = "", ylab = "", xlab = "",
     cex = 0.4, pch = 19, col = grey(.6), ylim = c(0,1))
lines(ecdf(sapply(results3, `[[`, "OFvalue")),
     cex = 0.4, pch = 19, col = grey(.4))


###################################################
### code chunk number 70: NMOFman.Rnw:1232-1237
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0), ps = 8, tck = 0.001)
plot( ecdf(sapply(results3, `[[`, "OFvalue")), main = "", ylab = "", xlab = "",
     cex = 0.4, pch = 19, col = grey(.4), ylim = c(0,1))
lines(ecdf(sapply(results4, `[[`, "OFvalue")),
     cex = 0.4, pch = 19, col = grey(.2))


###################################################
### code chunk number 71: NMOFman.Rnw:1249-1254
###################################################
OFlts <- function(b, Data) {
    tmp <- Data$y - Data$X %*% b
    tmp <- sort(tmp * tmp, partial = Data$h)
    sum(tmp[seq_len(Data$h)])
}


###################################################
### code chunk number 72: NMOFman.Rnw:1262-1266
###################################################
require("robustbase")
alpha <- 0.9
Data <- list(X = rD$X, y = rD$y, p = rD$p, n = rD$n,
             h = h.alpha.n(alpha, n = n, p = p))


###################################################
### code chunk number 73: NMOFman.Rnw:1273-1278
###################################################
resDE  <- DEopt(OFlts, algo = algo, Data = Data)
resLTS <- ltsReg(rD$y ~ -1 + rD$X, alpha = alpha,
                use.correction = FALSE)
data.frame(fastLTS = resLTS$raw.coefficients,
           DE = resDE$xbest)


###################################################
### code chunk number 74: NMOFman.Rnw:1282-1288
###################################################
cLTS <- resLTS$raw.coefficients
cat("LTS")
sum(sort((Data$X %*%cLTS - Data$y)^2)[1:Data$h])
cDE <- resDE$xbest
cat("DEopt")
sum(sort((Data$X %*%cDE -  Data$y)^2)[1:Data$h])


###################################################
### code chunk number 75: NMOFman.Rnw:1296-1298
###################################################
any(b0 < 0)
sum(b0)


###################################################
### code chunk number 76: NMOFman.Rnw:1307-1311
###################################################
repair <- function(b, Data) {
    b <- abs(b)
    b/sum(b)
}


###################################################
### code chunk number 77: NMOFman.Rnw:1318-1321
###################################################
b1 <- repair(b0, Data)
all(b1 >= 0)  ## should be TRUE
sum(b1)       ## should be 1


###################################################
### code chunk number 78: NMOFman.Rnw:1331-1334
###################################################
b0
b0 - abs(b0)
sum(b0)


###################################################
### code chunk number 79: NMOFman.Rnw:1341-1347
###################################################
b <- rnorm(1000L)
benchmark(ignore1 <- sum(b - abs(b))/2,
          ignore2 <- sum(b[b < 0]),
          columns = c("test", "elapsed", "relative"),
          replications = 1e4, order = "relative")
all.equal(ignore1, ignore2)


###################################################
### code chunk number 80: NMOFman.Rnw:1351-1352
###################################################
abs(sum(b0) - 1)


###################################################
### code chunk number 81: NMOFman.Rnw:1357-1363
###################################################
Data$pw1 <- 100
Data$pw2 <- 100
penalty <- function(b, Data)
    Data$pw1 * -sum(b - abs(b)) + Data$pw2 * abs(sum(b) - 1)
penalty(b0, Data)
penalty(b1, Data)


###################################################
### code chunk number 82: NMOFman.Rnw:1368-1372
###################################################
algo <- list(nG = 300, nP = 50,
             min = rep(-20, p), max = rep( 20, p),
             printBar = FALSE)
resDE <- DEopt(OFls, algo = algo, Data = Data)


###################################################
### code chunk number 83: NMOFman.Rnw:1377-1381
###################################################
round(resDE$xbest, 5)
resDE$OFvalue
sum(resDE$xbest)
all(resDE$xbest >= 0)


###################################################
### code chunk number 84: NMOFman.Rnw:1385-1391
###################################################
algo$repair <- repair
resDE <- DEopt(OFls, algo = algo, Data = Data)
round(resDE$xbest,5)
resDE$OFvalue
sum(resDE$xbest)
all(resDE$xbest >= 0)


###################################################
### code chunk number 85: NMOFman.Rnw:1394-1401
###################################################
algo$repair <- NULL
algo$pen <- penalty
resDE <- DEopt(OFls, algo = algo, Data = Data)
round(resDE$xbest,20)
resDE$OFvalue
sum(resDE$xbest)
all(resDE$xbest >= 0)


###################################################
### code chunk number 86: NMOFman.Rnw:1415-1416
###################################################
OFlts


###################################################
### code chunk number 87: NMOFman.Rnw:1419-1422
###################################################
b0 <- rnorm(p)
b1 <- rnorm(p)
P <- cbind(b0 = b0, b1 = b1)


###################################################
### code chunk number 88: NMOFman.Rnw:1426-1429
###################################################
head(Data$y - Data$X %*% b0)
head(Data$y - Data$X %*% b1)
head(drop(Data$y) - Data$X %*% P)


###################################################
### code chunk number 89: NMOFman.Rnw:1438-1441
###################################################
head(Data$y - Data$X %*% b0)^2
head(Data$y - Data$X %*% b1)^2
head((drop(Data$y) - Data$X %*% P)*(drop(Data$y) - Data$X %*% P))


###################################################
### code chunk number 90: NMOFman.Rnw:1444-1450
###################################################
OFlts2 <- function(b, Data) {
    tmp <- drop(Data$y) - Data$X %*% b
    tmp <- tmp * tmp
    tmp <- apply(tmp, 2L, sort, partial = Data$h)
    .colSums(tmp[seq_len(Data$h), ,drop = FALSE], Data$h, ncol(b))
}


###################################################
### code chunk number 91: NMOFman.Rnw:1454-1457
###################################################
nP <- 100
P <- array(rnorm(p * nP), dim = c(p, nP))
sol0 <- OFlts2(P, Data)


###################################################
### code chunk number 92: NMOFman.Rnw:1460-1467
###################################################
sol1 <- numeric(nP)
benchmark(for (i in seq_len(nP))
          sol1[i] <- OFlts(P[ , i, drop = FALSE], Data),
          sol2 <- OFlts2(P, Data),
          columns = c("test", "elapsed", "relative"),
          replications = 100, order = "relative")
all.equal(ignore1, ignore2)


###################################################
### code chunk number 93: NMOFman.Rnw:1600-1607
###################################################
OF <- tfTrefethen
n <- 100L
surf <- matrix(NA, n, n)
x1 <- seq(from = -10, to = 10, length.out = n)
for (i in seq_len(n))
    for (j in seq_len(n))
        surf[i, j] <- tfTrefethen(c(x1[i], x1[j]))


###################################################
### code chunk number 94: NMOFman.Rnw:1615-1621
###################################################
par(bty = "n", las = 1, mar = c(3,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))

## the actual minimum
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))


###################################################
### code chunk number 95: NMOFman.Rnw:1627-1638
###################################################
algo <- list(nP = 50L,
             nG = 300L,
             F = 0.6,
             CR = 0.9,
             min = c(-10,-10),
             max = c(10,10),
             printDetail = FALSE,
             printBar = FALSE,
             storeF = TRUE,
             storeSolutions = TRUE)
sol <- DEopt(OF = OF, algo = algo)


###################################################
### code chunk number 96: NMOFman.Rnw:1643-1649
###################################################
names(sol)
sd(sol$popF)
ts.plot(sol$Fmat, xlab = "generations", ylab = "OF")

length(sol$xlist)
xlist <- sol$xlist[[1L]]


###################################################
### code chunk number 97: NMOFman.Rnw:1658-1673
###################################################
## show solution 1 (column 1) in population over time
xlist[[  1L]][ ,1L]  ## at the end of generation 1
## ...
xlist[[ 10L]][ ,1L]  ## at the end of generation 10
## ...
xlist[[300L]][ ,1L]  ## at the end of generation 300

res  <- sapply(xlist, `[`, 1:2, 1)  ## get row 1 and 2 from column 1
res2 <- sapply(xlist, `[`, TRUE, 1) ## simpler
all.equal(res, res2)

dim(res)
res[ ,1L]
res[ ,2L]
res[ ,300L]


###################################################
### code chunk number 98: NMOFman.Rnw:1679-1693
###################################################
## show parameter 2 (row 2) in population over time
xlist[[  1L]][2L, ]  ## at the end of generation 1
## ...
xlist[[ 10L]][2L, ]  ## at the end of generation 10
## ...
xlist[[300L]][2L, ]  ## at the end of generation 300

res <- sapply(xlist, `[`, 2, 1:50)
res <- sapply(xlist, `[`, 2, TRUE)  ## simpler

dim(res)
res[ ,1L]
res[ ,2L]
res[ ,300L]


###################################################
### code chunk number 99: NMOFman.Rnw:1699-1708 (eval = FALSE)
###################################################
## ## transposing xlist[[i]] gives a two-column matrix -- see ?points
## ## initial solutions
## points(t(xlist[[1L]]), pch = 21, bg=grey(0.9), col = grey(.2))
## 
## ## solutions at the end of generation 100
## points(t(xlist[[100L]]), pch = 21, bg=grey(0.9), col = grey(.2))
## 
## ## solutions at the end of generation 100
## points(t(xlist[[300L]]), pch = 21, bg=grey(0.9), col = grey(.2))


###################################################
### code chunk number 100: NMOFman.Rnw:1712-1741
###################################################
setEPS()
postscript(file = "figures/c1.eps", width = 2, height = 2)
par(bty = "n", las = 1,mar = c(2,2,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))
points(t(xlist[[1L]]), pch = 21, bg=grey(0.9), col = grey(.2))
invisible(dev.off())

setEPS()
postscript(file = "figures/c2.eps", width = 2, height = 2)
par(bty = "n", las = 1,mar = c(2,2,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))
points(t(xlist[[100L]]), pch = 21, bg=grey(0.9), col = grey(.2))
invisible(dev.off())

setEPS()
postscript(file = "figures/c3.eps", width = 2, height = 2)
par(bty="n", las = 1,mar = c(2,2,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))
points(t(xlist[[300L]]), pch = 21, bg=grey(0.9), col = grey(.2))
invisible(dev.off())
cat("\\includegraphics{figures/c1.eps}",
    "\\includegraphics{figures/c2.eps}",
    "\\includegraphics{figures/c3.eps}", sep = "")


###################################################
### code chunk number 101: NMOFman.Rnw:1755-1770
###################################################
OF <- function(par, Data) {
    ## compute model yields
    y <- Data$model(par, Data$tm)

    ## all rates finite?
    validRates <- !any(is.na(y))

    if (validRates) {
        ## any rates negative? if yes, add penalty
        pen1 <- sum(abs(y - abs(y))) * Data$ww

        F <- max(abs(y - Data$yM)) + pen1
    } else F <- 1e8
    F
}


###################################################
### code chunk number 102: NMOFman.Rnw:1777-1794
###################################################
algo <- list(nP = 200L, nG = 100L,
             F = 0.50, CR = 0.99,
             min = c( 0,-10,-10,  0),
             max = c( 1, 10, 10, 10),
             storeSolutions = TRUE, printBar = FALSE)


## set up yield curve and put information in Data
tm <- 1:20                ## times to maturity
parTRUE <- c(5, 3, 2, 1)  ## true parameters
yM <- NS(parTRUE, tm)     ## true market yields
Data <- list(yM = yM, tm = tm, model = NS, ww = 0.1, maxb1 = 4)

## solve with DEopt
sol <- DEopt(OF = OF, algo = algo, Data = Data)
P <- sol$xlist[[1L]] ## all population matrices
p1 <- sapply(P, `[`, 1L, TRUE)


###################################################
### code chunk number 103: NMOFman.Rnw:1801-1810
###################################################
par(bty = "n", las = 1, mar = c(4,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
plot(jitter(rep(seq_len(algo$nG), each = algo$nP), factor = 5),
     p1,
     pch = 21, cex = 0.01, ylim = c(-5,10),
     xlab = "", ylab = "")

mtext("generation", 1, line = 2)
mtext("parameter\nvalue", 2, line = 1)


###################################################
### code chunk number 104: NMOFman.Rnw:1818-1853
###################################################
OF2 <- function(par, Data) {
    ## compute model yields
    y <- Data$model(par, Data$tm)

    ## all rates finite?
    validRates <- !any(is.na(y))

    if (validRates) {
        ## any rates negative? if yes, add penalty
        pen1 <- sum(abs(y - abs(y))) * Data$ww

        ## is b1 greater than Data$maxb1? if yes, add penalty
        pen2 <- par[1L] - Data$maxb1
        pen2 <- pen2 + abs(pen2)
        pen2 <- pen2

        F <- max(abs(y - Data$yM)) + pen1 + pen2
    } else F <- 1e8
    F
}

## solve with DEopt
sol <- DEopt(OF = OF2, algo = algo, Data = Data)
P <- sol$xlist[[1L]] ### all population matrices
p1 <- sapply(P, `[`, 1, TRUE)

par(bty = "n", las = 1, mar = c(4,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
plot(jitter(rep(seq_len(algo$nG), each = algo$nP), factor = 5),
     p1,
     pch = 21, cex = 0.01, ylim = c(-5,10),
     xlab = "", ylab = "")
abline(h = 4, col=grey(0.5))
mtext("generation", 1, line = 2)
mtext("parameter\nvalue", 2, line = 1)


###################################################
### code chunk number 105: NMOFman.Rnw:1944-1945
###################################################
tfRosenbrock


###################################################
### code chunk number 106: NMOFman.Rnw:1952-1956
###################################################
OF <- tfRosenbrock     ## see ?testFunctions
size <- 5L             ## set dimension
x <- rep.int(1, size)  ## the known solution ...
OF(x)                  ## ... should give zero


###################################################
### code chunk number 107: NMOFman.Rnw:1964-1971
###################################################
algo <- list(printBar = FALSE,
             nP = 50L,
             nG = 500L,
             F = 0.6,
             CR = 0.9,
             min = rep(-100, size),
             max = rep( 100, size))


###################################################
### code chunk number 108: NMOFman.Rnw:1978-1984
###################################################
## a vectorised OF: works only with *matrix* x
OF2 <- function(x) {
    n <- dim(x)[1L]
    xi <- x[1L:(n - 1L), ]
    colSums(100 * (x[2L:n, ] - xi * xi)^2 + (1 - xi)^2)
}


###################################################
### code chunk number 109: NMOFman.Rnw:1989-1993
###################################################
x <- matrix(rnorm(size * algo$nP), size, algo$nP)
c(OF(x[ ,1L]), OF(x[ ,2L]), OF(x[ ,3L]))
OF2(x)[1L:3L]  ## should give the same result
all.equal(OF2(x)[1L:3L], c(OF(x[ ,1L]), OF(x[ ,2L]), OF(x[ ,3L])))


###################################################
### code chunk number 110: NMOFman.Rnw:2004-2010
###################################################
set.seed(1223445)
(t1 <- system.time(sol <- DEopt(OF = OF, algo = algo)))

algo$loopOF <- FALSE
set.seed(1223445)
(t2 <- system.time(sol2 <- DEopt(OF = OF2, algo = algo)))


###################################################
### code chunk number 111: NMOFman.Rnw:2015-2018
###################################################
sol$OFvalue    ## both should be zero (with luck)
sol2$OFvalue
t1[[3L]]/t2[[3L]]  ## speedup


###################################################
### code chunk number 112: NMOFman.Rnw:2050-2061
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
### code chunk number 113: NMOFman.Rnw:2065-2085
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
### code chunk number 114: NMOFman.Rnw:2090-2092
###################################################
all.equal(res1,res2)
all.equal(res2,res3)


###################################################
### code chunk number 115: NMOFman.Rnw:2097-2100
###################################################
t1  ##  speedup for variant 1
t2  ##  speedup for variant 2
t3  ##  speedup for variant 3


###################################################
### code chunk number 116: NMOFman.Rnw:2125-2139
###################################################
n <- 100L   ## number of observation
p <- 5L     ## number of regressors
np <- 100L  ## population size
trials <- seq_len(1000L)

## random data
X <- array(rnorm(n * p), dim = c(n, p))
y <- rnorm(n)

## random population
Theta <- array(rnorm(p * np), dim = c(p, np))

## empty residuals matrix
R1 <- array(NA, dim = c(n, np))


###################################################
### code chunk number 117: NMOFman.Rnw:2143-2152
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
### code chunk number 118: NMOFman.Rnw:2159-2160
###################################################
all.equal(R1, R2)  ## ... should be TRUE


###################################################
### code chunk number 119: NMOFman.Rnw:2195-2204
###################################################
testFun <- function(x) {
    Sys.sleep(0.1)  ## wasting time...
    cos(1/x^2)
}
system.time(sol1 <- bracketing(testFun, interval = c(0.3, 0.9),
                               n = 100L))
system.time(sol2 <- bracketing(testFun, interval = c(0.3, 0.9),
                               n = 100L, method = "snow", cl = 2))
all.equal(sol1, sol2)


###################################################
### code chunk number 120: NMOFman.Rnw:2228-2232
###################################################
testFun  <- function(x) {
    Sys.sleep(0.1)  ## wasting time...
    x[1L] + x[2L]^2
}


###################################################
### code chunk number 121: NMOFman.Rnw:2236-2240
###################################################
lower <- c(1, 3); upper <- 5; n <- 5L
system.time(sol1 <- gridSearch(fun = testFun,
                               lower = lower, upper = upper,
                               n = n, printDetail = TRUE))


###################################################
### code chunk number 122: NMOFman.Rnw:2244-2246
###################################################
seq(from = 1, to = 5, length.out= n)  ## x_1
seq(from = 3, to = 5, length.out= n)  ## x_2


###################################################
### code chunk number 123: NMOFman.Rnw:2250-2252
###################################################
sol1$minfun
sol1$minlevels


###################################################
### code chunk number 124: NMOFman.Rnw:2256-2263
###################################################
system.time(sol2 <- gridSearch(fun = testFun,
                               lower = lower,
                               upper = upper,
                               n = n, printDetail = FALSE,
                               method = "snow",  ## use 'snow' ...
                               cl = 2L))         ## ... with 2 cores
all.equal(sol1, sol2)


###################################################
### code chunk number 125: NMOFman.Rnw:2303-2304
###################################################
cfBSM


###################################################
### code chunk number 126: NMOFman.Rnw:2308-2327
###################################################
S <- 100    ## spot
X <- 100    ## strike
tau <- 1    ## time-to-maturity
r <- 0.02   ## interest rate
q <- 0.08   ## dividend rate
v <- 0.2    ## volatility

## the closed-form solution
callBSM <- function(S,X,tau,r,q,v) {
    d1 <- (log(S/X) + (r - q + v^2 / 2)*tau) / (v*sqrt(tau))
    d2 <- d1 - v*sqrt(tau)
    S * exp(-q * tau) * pnorm(d1) -  X * exp(-r * tau) * pnorm(d2)
}
callBSM(S,X,tau,r,q,v)

## with the characteristic function
callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
       v = v^2,  ## variance, not vol
       implVol = TRUE)


###################################################
### code chunk number 127: NMOFman.Rnw:2343-2348
###################################################
cf <- c(5, 5, 5, 5, 5, 105) ## cashflows
times <- 1:6                ## times to payment
y <- 0.047                  ## the "true" yield
b0 <- sum(cf/(1 + y)^times)
b0


###################################################
### code chunk number 128: NMOFman.Rnw:2357-2364
###################################################
vanillaBond <- function(cf, times, df, yields) {
    if (missing(times))
        times <- seq_len(length(cf))
    if (missing(df))
        df <- 1/(1+yields)^times      
    drop(cf %*% df)
}


###################################################
### code chunk number 129: NMOFman.Rnw:2369-2372
###################################################
cf <- c(rep(5, 9), 105)
vanillaBond(cf, yields = 0.05)
vanillaBond(cf, yields = 0.03)


###################################################
### code chunk number 130: NMOFman.Rnw:2381-2382
###################################################
2^(1:5)


###################################################
### code chunk number 131: NMOFman.Rnw:2393-2394
###################################################
vanillaBond(cf, 1:10, yield = NS(c(0.03,0,0,2), 1:10))


###################################################
### code chunk number 132: NMOFman.Rnw:2403-2409
###################################################
  cf <- c(5, 5, 5, 5, 5, 105) ## cashflows
  times <- 1:6                ## times to payment
  y <- 0.047                  ## the "true" yield
  b0 <- sum(cf/(1 + y)^times)
  cf <- c(-b0, cf); times <- c(0, times)  
  data.frame(times=times, cashflows=cf)


###################################################
### code chunk number 133: NMOFman.Rnw:2415-2451
###################################################
ytm <- function(cf, times, y0 = 0.05,
                tol = 1e-05, h = 1e-05, maxit = 1000L) {        
    dr <- 1
    for (i in seq_len(maxit)) {
        y1 <- 1 + y0
        g <- cf / y1 ^ times
        g <- sum(g)
        t1 <- times - 1
        dg <- times * cf * 1/y1 ^ t1
        dg <- sum(dg)
        dr <- g/dg
        y0 <- y0 + dr
        if (abs(dr) < tol)
            break
    }
    y0
}
ytm2 <- function(cf, times, y0 = 0.05,
                 tol = 1e-04, h = 1e-08, maxit = 1000L) {        
    dr <- 1
    for (i in seq_len(maxit)) {
        y1 <- 1 + y0
        g <- sum(cf/y1^times)
        y1 <- y1 + h
        dg <- (sum(cf/y1^times) - g)/h
        dr <- g/dg
        y0 <- y0 - dr
        if (abs(dr) < tol)
            break
    }
    y0
}
system.time(for (i in 1:2000) ytm(cf, times, y0=0.06))
system.time(for (i in 1:2000) ytm2(cf, times, y0=0.06))
ytm(cf, times, y0=0.062, maxit = 5000)
ytm2(cf, times, y0=0.062, maxit = 5000)


###################################################
### code chunk number 134: NMOFman.Rnw:2462-2463
###################################################
(initial.value <- 5/b0)


###################################################
### code chunk number 135: NMOFman.Rnw:2466-2470
###################################################
ytm(cf,  times, y0 = 0.7, maxit = 5000)
ytm(cf,  times, y0 = initial.value)
ytm2(cf, times, y0 = 0.7, maxit = 5000)
ytm2(cf, times, y0 = initial.value)


###################################################
### code chunk number 136: NMOFman.Rnw:2516-2532
###################################################
mu <- 1
sigma <- 2
a <- -1
b <-  4

at <- (a - mu)/sigma
bt <- (b - mu)/sigma

## "throw away" strategy
x0 <- rnorm(10000L, mean = mu, sd = sigma)
x0 <- x0[x0>=a & x0<=b]

## inverse strategy
u <- runif(length(x0))
z <- qnorm(pnorm(at) + u*(pnorm(bt) - pnorm(at)))
x1 <- z * sigma + mu


###################################################
### code chunk number 137: NMOFman.Rnw:2537-2542
###################################################
par(mfrow = c(1, 2), mar = c(3, 3, 1, 1),
    bty = "n", las = 1, ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))

hist(x0, xlab = "")
hist(x1, xlab = "")


###################################################
### code chunk number 138: NMOFman.Rnw:2546-2550
###################################################
x1 <- x1[1:750]
x2 <- rnorm(200)
x3 <- runif(500)
x4 <- rbinom(100, size = 50, prob = 0.4)


###################################################
### code chunk number 139: NMOFman.Rnw:2554-2556
###################################################
cormat <- array(0.5, dim = c(4, 4))
diag(cormat) <- 1


###################################################
### code chunk number 140: NMOFman.Rnw:2560-2563
###################################################
results <- resampleC(x1 = x1, x2 = x2, x3 = x3, x4 = x4,
                     size = 50, cormat = cormat)
cor(results, method = "spearman")


###################################################
### code chunk number 141: NMOFman.Rnw:2566-2580
###################################################
## this function is taken from ?pairs
panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1L], y, col = grey(.5))
}
par(mar = c(3, 3, 1, 1),
    bty = "n", las = 1, ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
pairs(results,
      diag.panel = panel.hist,
      gap = 0, pch = 19, cex = 0.5)


###################################################
### code chunk number 142: NMOFman.Rnw:2585-2595
###################################################
par(mfrow = c(2, 4), mar = c(3, 5, 1, 1),
    bty = "n", las = 1, ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
hist(x1, xlab = "", ylab = "original")
hist(x2, xlab = "", ylab = "")
hist(x3, xlab = "", ylab = "")
hist(x4, xlab = "", ylab = "")
hist(results[ ,"x1"], xlab = "", ylab = "resampled")
hist(results[ ,"x2"], xlab = "", ylab = "")
hist(results[ ,"x3"], xlab = "", ylab = "")
hist(results[ ,"x4"], xlab = "", ylab = "")


###################################################
### code chunk number 143: NMOFman.Rnw:2659-2667
###################################################
na <- 500L                      ## number of assets

C <- array(0.6, dim = c(na,na)) ## correlation matrix
diag(C) <- 1

minVol <- 0.20; maxVol <- 0.40  ## covariance matrix
Vols <- (maxVol - minVol) * runif(na) + minVol
Sigma <- outer(Vols, Vols) * C


###################################################
### code chunk number 144: NMOFman.Rnw:2678-2684
###################################################
OF <- function(x, Data) {
    sx <- sum(x)
    w <- rep.int(1/sx, sx)
    res <- crossprod(w, Data$Sigma[x, x])
    tcrossprod(w, res)
}


###################################################
### code chunk number 145: NMOFman.Rnw:2687-2697
###################################################
neighbour <- function(xc, Data) {
    xn <- xc
    p <- sample.int(Data$na, Data$nn, replace = FALSE)
    xn[p] <- !xn[p]

    ## reject infeasible solution
    sumx <- sum(xn)
    if ( (sumx > Data$Ksup) || (sumx < Data$Kinf) )
        xc else xn
}


###################################################
### code chunk number 146: NMOFman.Rnw:2703-2708
###################################################
Data <- list(Sigma = Sigma,   ## cov-matrix
             Kinf  = 30L,     ## min cardinality
             Ksup  = 60L,     ## max cardinality
             na    = na,      ## number of assets
             nn    = 1L)      ## how many assets to change per iteration


###################################################
### code chunk number 147: NMOFman.Rnw:2711-2715
###################################################
card0 <- sample(Data$Kinf:Data$Ksup, 1L, replace = FALSE)
assets <- sample.int(na, card0, replace = FALSE)
x0 <- logical(na)
x0[assets] <- TRUE


###################################################
### code chunk number 148: NMOFman.Rnw:2720-2728
###################################################
## Local Search
algo <- list(x0 = x0, neighbour = neighbour, nS = 5000L,
             printDetail = FALSE, printBar = FALSE)
system.time(solLS <- LSopt(OF, algo = algo, Data = Data))

## Threshold Accepting
algo$nT <- 10L; algo$nS <- trunc(algo$nS/algo$nT); algo$q <- 0.2
system.time(solTA <- TAopt(OF, algo = algo, Data = Data))


###################################################
### code chunk number 149: NMOFman.Rnw:2812-2819
###################################################
OF2 <- function(x, Data) {
    res <- colSums(Data$Sigma %*% x * x)
    n <- colSums(x); res <- res / n^2
    ## penalise
    p <- pmax(Data$Kinf - n, 0) + pmax(n - Data$Ksup, 0)
    res + p
}


###################################################
### code chunk number 150: NMOFman.Rnw:2826-2829
###################################################
algo <- list(nB = na, nP = 100L, nG = 500L, prob = 0.002,
             printBar = FALSE, loopOF = FALSE)
system.time(solGA <- GAopt(OF = OF2, algo = algo, Data = Data))


###################################################
### code chunk number 151: NMOFman.Rnw:2834-2838
###################################################
cat("Local Search        ", format(sqrt(solLS$OFvalue), digits = 4), "\n",
    "Threshold Accepting ", format(sqrt(solTA$OFvalue), digits = 4), "\n",
    "Genetic Algorithm   ", format(sqrt(solGA$OFvalue), digits = 4), "\n",
    sep = "")


###################################################
### code chunk number 152: NMOFman.Rnw:2867-2875
###################################################
na <- 100L                                 ## number of assets
ns <- 200L                                 ## number of scenarios
vols <- runif(na, min = 0.2, max = 0.4)    ## marginal vols
C <- matrix(0.6, na, na); diag(C) <- 1     ## correlation matrix
R <- rnorm(ns * na)/16                     ## random returns
dim(R) <- c(ns, na)
R <- R %*% chol(C)
R <- R %*% diag(vols)


###################################################
### code chunk number 153: NMOFman.Rnw:2900-2910
###################################################
Data <- list(R = t(R),              ## scenarios
             theta = 0.005,         ## return threshold
             na = na,               ## number of assets
             ns = ns,               ## number of scenarios
             max = rep( 0.05, na),  ## DE: vector of max. weight
             min = rep(-0.05, na),  ## DE: vector of min. weight
             wsup =  0.05,          ## TA: max weight
             winf = -0.05,          ## TA: min weight
             eps = 0.5/100,         ## TA: step size
             w = 1)                 ## penalty weight


###################################################
### code chunk number 154: NMOFman.Rnw:2917-2920
###################################################
x0 <- Data$min + runif(Data$na)*(Data$max - Data$min)
x0[1:5]
sum(x0)


###################################################
### code chunk number 155: NMOFman.Rnw:2926-2930
###################################################
temp <- R %*% x0             ## compute portfolio returns
temp <- temp - Data$theta    ## subtract return threshold
temp <- (temp[temp < 0])^2   ## select elements below threshold
sum(temp)/ns                 ## compute semivariance


###################################################
### code chunk number 156: NMOFman.Rnw:2936-2943
###################################################
OF <- function(x, Data) {
    Rx <- crossprod(Data$R, x)
    Rx <- Rx - Data$theta
    Rx <- Rx - abs(Rx)
    Rx <- Rx * Rx
    colSums(Rx) /(4*Data$ns)
}


###################################################
### code chunk number 157: NMOFman.Rnw:2951-2953
###################################################
OF(x0, Data)
OF(cbind(x0, x0, x0), Data)


###################################################
### code chunk number 158: NMOFman.Rnw:2961-2974
###################################################
repair <- function(x, Data) {
    myFun <- function(x)
        x/sum(x)
    if (is.null(dim(x)[2L]))
        myFun(x) else apply(x, 2L, myFun)
}

repair2 <- function(x, Data) {
    myFun <- function(x)
        x + (1 - sum(x))/Data$na
    if (is.null(dim(x)[2L]))
        myFun(x) else apply(x, 2L, myFun)
}


###################################################
### code chunk number 159: NMOFman.Rnw:2980-2986
###################################################
sum(x0)
sum(repair(x0, Data))
sum(repair2(x0, Data))

colSums(repair( cbind(x0, x0, x0), Data))
colSums(repair2(cbind(x0, x0, x0), Data))


###################################################
### code chunk number 160: NMOFman.Rnw:2992-2994
###################################################
summary(repair (x0, Data)-x0)
summary(repair2(x0, Data)-x0)


###################################################
### code chunk number 161: NMOFman.Rnw:2999-3010
###################################################
penalty <- function(x, Data) {
    up <- Data$max
    lo <- Data$min
    xadjU <- x - up
    xadjU <- xadjU + abs(xadjU)
    xadjL <- lo - x
    xadjL <- xadjL + abs(xadjL)
    if (is.null(dim(x)[2L]))
        Data$w * (sum(xadjU) + sum(xadjL)) else
    Data$w * (colSums(xadjU) + colSums(xadjL))
}


###################################################
### code chunk number 162: NMOFman.Rnw:3018-3024
###################################################
x0[1L] <- 0.30
penalty(x0, Data)
penalty(cbind(x0, x0, x0), Data)
x0[1L] <- 0
penalty(x0, Data)
penalty(cbind(x0, x0, x0), Data)


###################################################
### code chunk number 163: NMOFman.Rnw:3029-3042
###################################################
algo <- list(nP = 100,        ## population size
             nG = 1000,       ## number of generations
             F = 0.25,        ## step size
             CR = 0.9,
             min = Data$min,
             max = Data$max,
             repair = repair,
             pen = penalty,
             printBar = FALSE,
             printDetail = TRUE,
             loopOF = TRUE,      ## do not vectorise
             loopPen = TRUE,     ## do not vectorise
             loopRepair = TRUE)  ## do not vectorise


###################################################
### code chunk number 164: NMOFman.Rnw:3048-3055
###################################################
system.time(sol <- DEopt(OF = OF,algo = algo,Data = Data))
16 * 100 * sqrt(sol$OFvalue)   ## solution quality

## check constraints
all(all.equal(sum(sol$xbest), 1),  ## budget constraint
    sol$xbest <= Data$max,         ## holding size constraints
    sol$xbest >= Data$min)


###################################################
### code chunk number 165: NMOFman.Rnw:3065-3075
###################################################
## looping over the population
algo$loopOF <- TRUE; algo$loopPen <- TRUE; algo$loopRepair <- TRUE
t1 <- system.time(sol <- DEopt(OF = OF,algo = algo, Data = Data))

## evaluating the population in one step
algo$loopOF <- FALSE; algo$loopPen <- FALSE; algo$loopRepair <- FALSE
t2 <- system.time(sol <- DEopt(OF = OF,algo = algo, Data = Data))

## speedup
t1[[3L]]/t2[[3L]]


###################################################
### code chunk number 166: NMOFman.Rnw:3086-3099
###################################################
algo$printDetail <- FALSE
restartsDE <- restartOpt(fun = DEopt,      ## what function
                         n = 20L,          ## how many restarts
                         OF = OF,
                         algo = algo,
                         Data = Data,
                         method = "snow",  ## using package snow
                         cl = 2)           ## 2 cores

## extract best solution
OFvaluesDE <- sapply(restartsDE, `[[`, "OFvalue")
OFvaluesDE <- 16 * 100 * sqrt(OFvaluesDE)
weightsDE  <- sapply(restartsDE, `[[`, "xbest")


###################################################
### code chunk number 167: NMOFman.Rnw:3104-3109
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001)
plot(sort(OFvaluesDE), (seq_len(length(OFvaluesDE))) / length(OFvaluesDE),
     type = "S", ylim = c(0, 1), xlab = "", ylab = "")
mtext("OF value",  side = 1, line = 2)


###################################################
### code chunk number 168: NMOFman.Rnw:3113-3119
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001)
boxplot(t(weightsDE),
        outline = FALSE, boxwex = 0.4, ylim = c(-0.06,0.06))
mtext("assets",  side = 1, line = 2)
mtext("weights", side = 2, line = 1.3, las = 1, padj = -4)


###################################################
### code chunk number 169: NMOFman.Rnw:3132-3147
###################################################
algo$printDetail <- FALSE;  algo$nP <- 200L; restarts <- 20L
nGs <- c(500L, 1500L, 3000L)
lstOFvaluesDE <- list()
for (i in 1:3) {
    algo$nG <- nGs[i]
    restartsDE <- restartOpt(fun = DEopt,
                             n = restarts,
                             OF = OF,  algo = algo, Data = Data,
                             method = "snow", cl = 2)
    ## extract best solution
    OFvaluesDE <- sapply(restartsDE, `[[`, "OFvalue")
    OFvaluesDE <- 16 * 100 * sqrt(OFvaluesDE)
    lstOFvaluesDE[[i]] <- OFvaluesDE
}
res <- simplify2array(lstOFvaluesDE)


###################################################
### code chunk number 170: NMOFman.Rnw:3152-3166
###################################################
algo$repair <- repair2
lstOFvaluesDE <- list()
for (i in 1:3) {
    algo$nG <- nGs[i]
    restartsDE <- restartOpt(fun = DEopt,
                             n = restarts,
                             OF = OF,  algo = algo, Data = Data,
                             method = "snow", cl = 2)
    ## extract best solution
    OFvaluesDE <- sapply(restartsDE, `[[`, "OFvalue")
    OFvaluesDE <- 16 * 100 * sqrt(OFvaluesDE)
    lstOFvaluesDE[[i]] <- OFvaluesDE
}
res2 <- simplify2array(lstOFvaluesDE)


###################################################
### code chunk number 171: NMOFman.Rnw:3171-3181
###################################################
allres <- as.vector(rbind(res,res2))
xlims <- pretty(allres); xlims <- c(min(xlims), max(xlims))
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001)
plot(ecdf(res[ ,3L]), xlim = xlims, cex = 0.4,
     main = "", ylab = "", xlab = "")
for (i in 1:2)
    lines(ecdf(res[  ,i]), cex = 0.4)
for (i in 1:3)
    lines(ecdf(res2[ ,i]), col = "blue", cex = 0.4)


###################################################
### code chunk number 172: NMOFman.Rnw:3191-3198
###################################################
weightsDE <- sapply(restartsDE, `[[`, "xbest")
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
ps = 8, tck = 0.001)
boxplot(t(weightsDE),
outline = FALSE, boxwex = 0.4, ylim = c(-0.06, 0.06))
mtext("assets",  side = 1, line = 2)
mtext("weights", side = 2, line = 1.3, las = 1, padj = -4)


###################################################
### code chunk number 173: NMOFman.Rnw:3212-3229
###################################################
algo <- list(nP = 100L,      ## population size
             nG = 1000L,     ## number of generations
             c1 = 0.5,       ## weight for individually best solution
             c2 = 1.5,       ## weight for overall best solution
             min = Data$min,
             max = Data$max,
             repair = repair, pen = penalty,
             iner = 0.7, initV = 1, maxV = 0.2,
             printBar = FALSE, printDetail = TRUE)

system.time(sol <- PSopt(OF = OF,algo = algo,Data = Data))
16 * 100 * sqrt(sol$OFvalue)      ## solution quality

## check constraints
all(all.equal(sum(sol$xbest),1),  ## budget constraint
sol$xbest <= Data$max,
sol$xbest >= Data$min)


###################################################
### code chunk number 174: NMOFman.Rnw:3237-3244
###################################################
changeV <- function(x, Data) {
    myFun <- function(x) x - (sum(x))/Data$na
    if (is.null(dim(x)[2L]))
        myFun(x) else apply(x, 2L, myFun)
}
sum(changeV(x0, Data))
colSums(changeV(cbind(x0, x0, x0), Data))


###################################################
### code chunk number 175: NMOFman.Rnw:3249-3253
###################################################
initP <- Data$min + diag(Data$max - Data$min) %*%
    array(runif(length(Data$min) * algo$nP),
          dim = c(length(Data$min),  algo$nP))
colSums(initP <- repair(initP,Data))[1:10] ## check


###################################################
### code chunk number 176: NMOFman.Rnw:3259-3265
###################################################
algo$changeV <- changeV        ## function to adjust velocity
algo$initP <- initP            ## initial population
algo$repair <- NULL            ## not needed anymore

system.time(sol <- PSopt(OF = OF,algo = algo, Data = Data))
16 * 100 * sqrt(sol$OFvalue)   ## solution quality


###################################################
### code chunk number 177: NMOFman.Rnw:3270-3273
###################################################
all(all.equal(sum(sol$xbest), 1), ## budget constraint
sol$xbest <= Data$max,
sol$xbest >= Data$min)


###################################################
### code chunk number 178: NMOFman.Rnw:3276-3279
###################################################
algo$loopOF <- FALSE; algo$loopPen <- FALSE
algo$loopRepair <- FALSE; algo$loopChangeV <- FALSE
system.time(sol <- PSopt(OF = OF, algo = algo, Data = Data))


###################################################
### code chunk number 179: NMOFman.Rnw:3284-3300
###################################################
algo$printDetail <- FALSE
restartsPS <- restartOpt(fun = PSopt,
                         n = 20L,
                         OF = OF,
                         algo = algo, Data = Data,
                         method = "snow", cl = 2)

## extract best solution
OFvaluesPS <- sapply(restartsPS, `[[`, "OFvalue")
OFvaluesPS <- 16 * 100 * sqrt(OFvaluesPS)
par(bty = "n", las = 1,mar = c(3,4,0,0),
    ps = 8, tck = 0.001)
plot(sort(OFvaluesPS),
     (seq_len(length(OFvaluesPS))) / length(OFvaluesPS),
     type = "S", ylim = c(0, 1), xlab = "", ylab = "")
mtext("OF value",  side = 1, line = 2)


###################################################
### code chunk number 180: NMOFman.Rnw:3311-3333
###################################################
Data$R <- R  ## not transposed any more

neighbourU <- function(sol, Data){
    resample <- function(x, ...)
        x[sample.int(length(x), ...)]
    wn <- sol$w
    toSell <- wn > Data$winf
    toBuy  <- wn < Data$wsup
    i <- resample(which(toSell), size = 1L)
    j <- resample(which(toBuy), size = 1L)
    eps <- runif(1) * Data$eps
    eps <- min(wn[i] - Data$winf, Data$wsup - wn[j], eps)
    wn[i] <- wn[i] - eps
    wn[j] <- wn[j] + eps
    Rw <- sol$Rw + Data$R[,c(i,j)] %*% c(-eps,eps)
    list(w = wn, Rw = Rw)
}
OF <- function(x, Data) {
    Rw <- x$Rw - Data$theta
    Rw <- Rw - abs(Rw)
    sum(Rw*Rw) / (4*Data$ns)
}


###################################################
### code chunk number 181: NMOFman.Rnw:3337-3349
###################################################
w0 <- runif(Data$na); w0 <- w0/sum(w0)
x0 <- list(w = w0, Rw = R %*% w0)
algo <- list(x0 = x0,
             neighbour = neighbourU,
             nS = 2000L,
             nT = 10L,
             nD = 5000L,
             q = 0.20,
             printBar = FALSE,
             printDetail = FALSE)
system.time(sol2 <- TAopt(OF,algo,Data))
16 * 100 * sqrt(sol2$OFvalue)


###################################################
### code chunk number 182: NMOFman.Rnw:3354-3376
###################################################
restartsTA <- restartOpt(fun = TAopt,
                         n = 20L,
                         OF = OF,
                         algo = algo,
                         Data = Data,
                         method = "snow",
                         cl = 2)

OFvaluesTA <- sapply(restartsTA, `[[`, "OFvalue") ## extract best solution
OFvaluesTA <- 16 * 100 * sqrt(OFvaluesTA)
weightsTA <- sapply(restartsTA, `[[`, "xbest")
par(bty = "n", las = 1,mar = c(3,4,0,0), ps = 8,
    tck = 0.001, mgp = c(3, 0.5, 0))

## blue: DE solution with nP = 200 and nG = 2000
xlims <- pretty(c(res2[,3], OFvaluesTA))
plot(ecdf(res2[,3]), col = "blue", cex = 0.4,
     main = "", ylab = "", xlab = "",
     xlim = c(min(xlims), max(xlims)) )

## black: TA
lines(ecdf(OFvaluesTA), cex = 0.4)


###################################################
### code chunk number 183: NMOFman.Rnw:3400-3406
###################################################
na <- 50
no <- 5000
D <- array(rnorm(na*no)*0.01, dim = c(no,na))
w <- runif(na)
w <- w/sum(w)
R <- D %*% w


###################################################
### code chunk number 184: NMOFman.Rnw:3411-3412
###################################################
require("compiler")


###################################################
### code chunk number 185: NMOFman.Rnw:3419-3425
###################################################
var1 <- function(R) {
    n <- NROW(R)
    m <- sum(R)/n
    crossprod(R)/(n-1) - m^2
}
var(R) - var1(R)


###################################################
### code chunk number 186: NMOFman.Rnw:3431-3444
###################################################
var1 <- function(R) {
    n <- NROW(R)
    m <- sum(R)/n
    crossprod(R)/(n - 1) - m^2
}
var2 <- cmpfun(var1)
runs <- 10000
system.time(for (i in seq_len(runs))
            ignore <- var(R))
system.time(for (i in seq_len(runs))
            ignore <- var1(R))  
system.time(for (i in seq_len(runs))
            ignore <- var2(R))  


###################################################
### code chunk number 187: NMOFman.Rnw:3451-3478
###################################################
pm0 <- function(x, xp = 2, threshold = 0, lower = TRUE) {    
    n <- NROW(x)
    x <- x - threshold
    if (lower)
        x <- x[x < 0] else x <- x[x > 0]
    sum(x^xp)/n    
}
pm1 <- function(x, xp = 2, threshold = 0, lower = TRUE, keep.sign = FALSE) {
    x <- x - threshold
    if (lower)
        x <- x - abs(x)
    else
        x <- x + abs(x)        
    sx <- sign(x)
    x <- abs(x)
    if (xp == 1L)
        sum(x)/2/length(x)
    else if (xp == 2L)
        sum(x*x)/4/length(x)
    else if (xp == 3L)
        sum(x*x*x)/8/length(x)
    else if (xp == 4L)
        sum(x*x*x*x)/16/length(x)
    else 
        sum(x^xp)/2^xp/length(x)
}
pm2 <- cmpfun(pm1)


###################################################
### code chunk number 188: NMOFman.Rnw:3484-3505
###################################################
pm0(R)
pm1(R)
pm2(R)
runs <- 1000
system.time(for (i in seq_len(runs))
            ignore <- pm0(R))
system.time(for (i in seq_len(runs))
            ignore <- pm1(R))  
system.time(for (i in seq_len(runs))
            ignore <- pm2(R))
 
pm0(R,2.2)
pm1(R,2.2)
pm2(R,2.2)
runs <- 1000
system.time(for (i in seq_len(runs))
            ignore <- pm0(R, 2.5))
system.time(for (i in seq_len(runs))
            ignore <- pm1(R, 2.5))  
system.time(for (i in seq_len(runs))
            ignore <- pm2(R, 2.5))  


###################################################
### code chunk number 189: NMOFman.Rnw:3542-3549
###################################################
size <- 20L
x <- logical(size)
x[runif(size) > 0.5] <- TRUE

## store information
Data <- list()
Data$size <- size


###################################################
### code chunk number 190: NMOFman.Rnw:3553-3562
###################################################
compareLogicals <- function(x, y, ...) {
    argsL <- list(...)
    if (!("sep" %in% names(argsL))) 
        argsL$sep <- ""
    do.call("cat",
            c(list("\n",as.integer(x), "\n", as.integer(y), "\n",
                   ifelse(x == y, " ", "^"), "\n"), argsL))
    message("The vectors differ in ", sum(x != y), " place(s).")
}


###################################################
### code chunk number 191: NMOFman.Rnw:3569-3575
###################################################
## there should be no difference
compareLogicals(x, x)

## change the second element
z <- x; z[2L] <- !z[2L]
compareLogicals(x, z)


###################################################
### code chunk number 192: NMOFman.Rnw:3584-3591
###################################################
Data$n <- 5L  ## how many elements to change
neighbour <- function(x, Data) {
    ii <- sample.int(Data$size, Data$n)
    x[ii] <- !x[ii]
    x
}
compareLogicals(x, neighbour(x, Data))


###################################################
### code chunk number 193: NMOFman.Rnw:3601-3612
###################################################
neighbour <- function(x, Data) {
    ## required: x must have at least one TRUE and one FALSE
    Ts <- which(x)
    Fs <- which(!x)
    lenTs <- length(Ts)
    O <- sample.int(lenTs,  1L)
    I <- sample.int(Data$size - lenTs, 1L)
    x[c(Fs[I], Ts[O])] <- c(TRUE, FALSE)
    x
}
compareLogicals(x, neighbour(x, Data))


###################################################
### code chunk number 194: NMOFman.Rnw:3619-3637
###################################################
size <- 5L
x0 <- runif(size)
xTRUE <- runif(size)
Data <- list(xTRUE = xTRUE,
             step = 0.02)
OF <- function(x, Data)
    max(abs(x - Data$xTRUE))
neighbour <- function(x, Data)
    x + runif(length(Data$xTRUE))*Data$step - Data$step/2

algo <- list(q = 0.05, nS = 1000L, nT = 10L,
             neighbour = neighbour, x0 = x0,
             printBar = FALSE,
             printDetail = FALSE,
             storeSolutions = TRUE,
             storeF = TRUE)
res <- TAopt(OF, algo = algo, Data = Data)
res$OFvalue < 0.005


###################################################
### code chunk number 195: NMOFman.Rnw:3658-3929
###################################################
## N1: This neighbour enforces a budget constraint, a minimum
## holding size (which need not be zero) and a maximum holding size

Data <- list(wmin = 0,     ## the minimal weight
             wmax = 0.22,  ## the maximal weight
             eps = 0.2/100,  ## the step size
             ## resample = function(x, ...)
             ##            x[sample.int(length(x), ...)],
             na = dim(fundData)[2L],
             R = fundData)

cat("The portfolio will consist of at least ",
    ceiling(1/Data$wmax), " assets.\n", sep = "")

neighbour1 <- function(w, Data){
    toSell <- which(w > Data$wmin)
    toBuy  <- which(w < Data$wmax)
    i <- toSell[sample.int(length(toSell), size = 1L)]
    j <- toBuy[ sample.int(length(toBuy),  size = 1L)]
    eps <- runif(1) * Data$eps
    eps <- min(w[i] - Data$wmin, Data$wmax - w[j], eps)
    w[i] <- w[i] - eps
    w[j] <- w[j] + eps
    w
}

neighbour1U <- function(x, Data){
    wn <- x$w
    toSell <- which(wn > Data$wmin)
    toBuy  <- which(wn < Data$wmax)
    i <- toSell[sample.int(length(toSell), size = 1L)]
    j <- toBuy[ sample.int(length(toBuy),  size = 1L)]
    eps <- runif(1) * Data$eps
    eps <- min(wn[i] - Data$wmin, Data$wmax - wn[j], eps)
    wn[i] <- wn[i] - eps
    wn[j] <- wn[j] + eps
    Rw <- x$Rw + Data$R[ ,c(i,j)] %*% c(-eps,eps)
    list(w = wn, Rw = Rw)
}

## create a random solution
makex <- function(Data) {
    resample <- function(x, ...) 
        x[sample.int(length(x), ...)]
    w0 <- numeric(Data$na)
    nAssets <- resample(ceiling(1/Data$wmax):Data$na, 1L)
    w0[sample(seq_len(Data$na), nAssets)] <- runif(nAssets)
    w0/sum(w0)
}

isOK <- function(w, Data) {
    tooBig   <- any(w > Data$wmax)
    tooSmall <- any(w < Data$wmin)
    sumToOne <- abs(sum(w)-1) < 1e-12
    if (!tooBig && !tooSmall && sumToOne)
        TRUE
    else
        FALSE
}

## TEST 1
w0 <- makex(Data)
x0 <- list(w = w0, Rw = fundData %*% w0)
isOK(w0, Data)
isOK(x0$w, Data)

set.seed(545)
w0 <- makex(Data)
nTests <- 1e3
for (i in seq(nTests)) {
    w1 <- neighbour1(w0,  Data)
    if (isOK(w1, Data))
        w0 <- w1
    else
        stop("error")
}

set.seed(545)
w0 <- makex(Data)
x0 <- list(w = w0, Rw = fundData %*% w0)
nTests <- 1e3
for (i in seq(nTests)) {
    x1 <- neighbour1U(x0,  Data)
    if (isOK(x1$w, Data))
        x0 <- x1
    else
        stop("error")
}
all.equal(fundData %*% w1, x1$Rw)


## TEST 2: reach a target solution
makeOF <- function(wt)
    function(w0, Data)
        sum(abs(wt - w0))
wt <- makex(Data)
OF <- makeOF(wt)
w0 <- makex(Data)
OF(w0, Data)
TAsettings <- list(neighbour = neighbour1,
                   x0 = w0, nS = 5000, q = 0.1,
                   printBar = FALSE)
res <- TAopt(OF, algo = TAsettings, Data)
round(head(sort(abs(res$xbest-wt), decreasing = TRUE),5),6)




## N2: This long-only neighbour enforces a budget constraint, a
## minimum holding size (which need not be zero), and a maximum holding
## size and a maximum cardinality.

Data <- list(wmax = 0.3, ## the maximal weight
             Kmax = 10,  ## max cardinality
             eps = 1/100, ## the step size
             ## resample = function(x, ...)
             ##                x[sample.int(length(x), ...)],
             na = dim(fundData)[2L],
             R = fundData)
cat("The portfolio will consist of at least ",
    ceiling(1/Data$wmax), " assets.\n", sep = "")

neighbour2 <- function(w, Data){
    tol <- 1e-12
    J <- sum(w > tol)
    if (J == Data$Kmax)
        toBuy <- which(w > tol & w < Data$wmax)
    else
        toBuy <- which(w < Data$wmax)
    toSell <- which(w > tol)
    i <- toSell[sample.int(length(toSell), size = 1L)]
    j <- toBuy[ sample.int(length(toBuy),  size = 1L)]
    eps <- runif(1) * Data$eps
    eps <- min(w[i], Data$wmax - w[j], eps)
    w[i] <- w[i] - eps
    w[j] <- w[j] + eps
    w
}
neighbour2U <- function(x, Data){
    tol <- 1e-12
    w <- x$w
    J <- sum(w > tol)
    if (J == Data$Kmax)
        toBuy <- which(w > tol & w < Data$wmax)
    else
        toBuy <- which(w < Data$wmax)
    toSell <- which(w > tol)
    i <- toSell[sample.int(length(toSell), size = 1L)]
    j <- toBuy[ sample.int(length(toBuy),  size = 1L)]
    eps <- runif(1) * Data$eps
    eps <- min(w[i], Data$wmax - w[j], eps)
    w[i] <- w[i] - eps
    w[j] <- w[j] + eps
    Rw <- x$Rw + Data$R[ ,c(i,j)] %*% c(-eps, eps)
    list(w = w, Rw = Rw)
}


makex <- function(Data) {
    w0 <- numeric(Data$na)
    nAssets <- sample(ceiling(1/Data$wmax):Data$Kmax, 1L)
    w0[sample(seq_len(Data$na), nAssets)] <- runif(nAssets)
    w0/sum(w0)
}

isOK <- function(w, Data) {
    tooBig   <- any(w > Data$wmax)
    tooMany <- sum(w > 1e-12) > Data$Kmax
    sumToOne <- abs(sum(w)-1) < 1e-12
    if (!tooBig && !tooMany && sumToOne)
        TRUE
    else
        FALSE
}

## TEST 1
w0 <- makex(Data)
x0 <- list(w = w0, Rw = fundData %*% w0)
isOK(w0, Data)
isOK(x0$w, Data)

set.seed(545)
w0 <- makex(Data)
nTests <- 1e3
for (i in seq(nTests)) {
    w1 <- neighbour2(w0,  Data)
    if (isOK(w1, Data))
        w0 <- w1
    else
        stop("error")
}

set.seed(545)
w0 <- makex(Data)
x0 <- list(w = w0, Rw = fundData %*% w0)
nTests <- 1e3
for (i in seq(nTests)) {
    x1 <- neighbour2U(x0,  Data)
    if (isOK(x1$w, Data))
        x0 <- x1
    else
        stop("error")
}

all.equal(fundData %*% w1, x1$Rw)

## TEST 2: reach a target solution
makeOF <- function(wt)
    function(w0, Data)
        sum(abs(wt - w0))
wt <- makex(Data)
OF <- makeOF(wt)
w0 <- makex(Data)
OF(w0, Data)
OF(wt, Data)

TAsettings <- list(neighbour = neighbour2,
                   x0 = w0, nS = 5000, q = 0.1,
                   printBar = FALSE)
res <- TAopt(OF, algo = TAsettings, Data)
isOK(res$xbest, Data)

df <- data.frame(target=wt, w0 = w0, wTAopt = res$xbest)
tmpfun <- function(x)
    !all(abs(x) < 1e-14)
df[apply(df,1,tmpfun),]
apply(df, 2, sum)

wt <- numeric(200)
wt[1:4] <- c(0.3,0.3,0.3,0.1)
OF <- makeOF(wt)

TAsettings <- list(neighbour = neighbour2,
                   x0 = w0, nS = 5000, q = 0.1,
                   printBar = FALSE)
res <- TAopt(OF, algo = TAsettings, Data)
isOK(res$xbest, Data)
df <- data.frame(target=wt, w0 = w0, wTAopt = res$xbest)
tmpfun <- function(x)
    !all(abs(x) < 1e-14)
df[apply(df,1,tmpfun),]
apply(df, 2, sum)


w0 <- makex(Data)
x0 <- list(w = w0, Rw = fundData %*% w0)

## the N is slower
system.time(for (i in 1:10000) neighbour2(w0, Data))
system.time(for (i in 1:10000) neighbour2U(x0, Data))


TAsettings2 <- list(neighbour = neighbour2,
                    x0 = w0, nS = 500, q = 0.1,
                    printBar = FALSE, printDetail = FALSE)
TAsettings2U <- list(neighbour = neighbour2U,
                     x0 = x0, nS = 500, q = 0.1,
                     printBar = FALSE, printDetail = FALSE)

ofun <- function(w, Data) {
    Rw <- Data$R %*% w
    crossprod(Rw)
}
ofunU <- function(sol, Data)
    crossprod(sol$Rw)
ign <- TAopt(ofun, TAsettings2, Data)
ign <- TAopt(ofunU, TAsettings2U, Data)

##benchmark(ign <- TAopt(ofun, TAsettings2, Data),
##          ign <- TAopt(ofunU, TAsettings2U, Data),
##          replications = 1, order = "relative")


###################################################
### code chunk number 196: NMOFman.Rnw:3964-3966
###################################################
maxabs <- function(x, lim)
    max(sum(abs(x)) - lim, 0)


###################################################
### code chunk number 197: NMOFman.Rnw:3971-4054
###################################################
require("compiler")

wmin <- 0.01
wmax <- 0.10


w <- numeric(50)
w[1:10] <- 0.1
step <- 0.01

N <- function(w) {

    wo <- w
    ## initial sell
    sel <- which(w >= wmin)
    i <- sel[sample.int(length(sel), size = 1)]

    if (w[i] == wmin) {
        eps <- wmin
        w[i] <- 0
    } else {
        eps <- min(runif(1)*step , w[i] - wmin)
        w[i] <- w[i] - eps
    }

    cash <- eps

    iter <- 0
    while (abs(cash) > 1e-14) {
        iter <- iter + 1
        if (iter > 10) {
            return(wo)
        }
        ##message("cash ", cash)
        if (cash > 0) {  ## buy something
            sel <- which(w < wmax)
            i <- sel[sample.int(length(sel), size = 1)]            
            if (w[i] == 0) {
                w[i] <- eps <- wmin
            } else {
                eps <- min(runif(1)*step , wmax - w[i], cash)
                w[i] <- w[i] + eps
            }
            cash <- cash - eps            
        } else { ## sell something
            sel <- which(w >= wmin)
            i <- sel[sample.int(length(sel), size = 1)]
            if (w[i] == wmin) {
                eps <- wmin
                w[i] <- 0
            } else {
                eps <- min(runif(1)*step , w[i] - wmin)
                w[i] <- w[i] - eps
            }
            cash <- cash + eps
        }
    }
    ##message(iter)
    w
}

##N <- cmpfun(N)
##w
system.time(for (i in 1:10000) w <- N(w))

## goal <- numeric(50)
## goal[41:50] <- 0.1

## OF <- function(x) {
##     tmp <- x - goal
##     sum(tmp * tmp)    
## }
## ans <- LSopt(OF, algo = list(nS = 1000000, neighbour = N, x0 = w))
## sum(ans$xbest > 0)
## ans$xbest

## ans <- numeric(10000)
## for (i in seq_along(ans))
##     ans[i] <- ceiling(runif(1)*10)

## system.time(for (i in 1:10000) ignore <- ceiling(runif(5)*99))
## system.time(for (i in 1:10000) ignore <- sample.int(99, 5))



###################################################
### code chunk number 198: NMOFman.Rnw:4069-4080
###################################################
makeCashFlows <- function(coupon, T) {
    t1 <- T - floor(T)                  ## time to first coupon
    tm <- seq(ifelse(t1 > 1e-5, t1, 1), ## 1e-5 is less than a calendar day
              T,
              by = 1)
    cf <- rep.int(coupon, length(tm))
    cf[length(cf)] <- cf[length(cf)] + 100
    list(cf = cf, tm = tm)
}

makeCashFlows(3, 10.2)


###################################################
### code chunk number 199: NMOFman.Rnw:4084-4104
###################################################
cf1 <- c(rep(5.75,  8), 105.75); tm1 <- 0:8 + 0.5
cf2 <- c(rep(4.25, 17), 104.25); tm2 <- 1:18
cf3 <- c(3.5, 103.5); tm3 <- 0:1 + 0.5
cf4 <- c(rep(3.00, 15), 103.00); tm4 <- 1:16
cf5 <- c(rep(3.25, 11), 103.25); tm5 <- 0:11 + 0.5
cf6 <- c(rep(5.75, 17), 105.75); tm6 <- 0:17 + 0.5
cf7 <- c(rep(3.50, 14), 103.50); tm7 <- 1:15
cf8 <- c(rep(5.00,  8), 105.00); tm8 <- 0:8 + 0.5
cf9 <- 105; tm9 <- 1
cf10 <- c(rep(3.00, 12), 103.00); tm10 <- 0:12 + 0.5
cf11 <- c(rep(2.50,  7), 102.50); tm11 <- 1:8
cf12 <- c(rep(4.00, 10), 104.00); tm12 <- 1:11
cf13 <- c(rep(3.75, 18), 103.75); tm13 <- 0:18 + 0.5
cf14 <- c(rep(4.00, 17), 104.00); tm14 <- 1:18
cf15 <- c(rep(2.25,  8), 102.25); tm15 <- 0:8 + 0.5
cf16 <- c(rep(4.00,  6), 104.00); tm16 <- 1:7
cf17 <- c(rep(2.25, 12), 102.25); tm17 <- 1:13
cf18 <- c(rep(4.50, 19), 104.50); tm18 <- 0:19 + 0.5
cf19 <- c(rep(2.25,  7), 102.25); tm19 <- 1:8
cf20 <- c(rep(3.00, 14), 103.00); tm20 <- 1:15


###################################################
### code chunk number 200: NMOFman.Rnw:4110-4125
###################################################
cfList <- list( cf1, cf2, cf3, cf4, cf5, cf6, cf7, cf8, cf9,cf10,
               cf11,cf12,cf13,cf14,cf15,cf16,cf17,cf18,cf19,cf20)
tmList <- list( tm1, tm2, tm3, tm4, tm5, tm6, tm7, tm8, tm9,tm10,
               tm11,tm12,tm13,tm14,tm15,tm16,tm17,tm18,tm19,tm20)
tm <- unlist(tmList, use.names = FALSE)
tm <- sort(unique(tm))
nR <- length(tm)
nC <- length(cfList)

cfMatrix <- array(0, dim = c(nR, nC))
for(j in seq(nC))
    cfMatrix[tm %in% tmList[[j]], j] <- cfList[[j]]
rownames(cfMatrix) <- tm

cfMatrix[1:10, 1:10]


###################################################
### code chunk number 201: NMOFman.Rnw:4132-4136
###################################################
betaTRUE <- c(5,-2,1,10,1,3)
yM <- NSS(betaTRUE,tm)
diFa <- 1 / ( (1 + yM/100)^tm )
bM <- diFa %*% cfMatrix


###################################################
### code chunk number 202: NMOFman.Rnw:4140-4144
###################################################
Data <- list(bM = bM, tm = tm, cfMatrix = cfMatrix, model = NSS,
             ww = 1,
             min = c( 0,-15,-30,-30,0  ,2.5),
             max = c(15, 30, 30, 30,2.5,5  ))


###################################################
### code chunk number 203: NMOFman.Rnw:4154-4164
###################################################
OF2 <- function(param, Data) {
  tm <- Data$tm
  bM <- Data$bM
  cfMatrix <- Data$cfMatrix
  diFa  <- 1 / ((1 + Data$model(param, tm)/100)^tm)
  b <- diFa %*% cfMatrix
  aux <- b - bM; aux <- max(abs(aux))
  if (is.na(aux)) aux <- 1e10
  aux
}


###################################################
### code chunk number 204: NMOFman.Rnw:4168-4183
###################################################
penalty <- function(mP, Data) {
    minV <- Data$min
    maxV <- Data$max
    ww <- Data$ww
    ## if larger than maxV, element in A is positiv
    A <- mP - as.vector(maxV)
    A <- A + abs(A)
    ## if smaller than minV, element in B is positiv
    B <- as.vector(minV) - mP
    B <- B + abs(B)
    ## beta 1 + beta2 > 0
    C <- ww*((mP[1L, ] + mP[2L, ]) - abs(mP[1L, ] + mP[2L, ]))
    A <- ww * colSums(A + B) - C
    A
}


###################################################
### code chunk number 205: NMOFman.Rnw:4188-4204
###################################################
algo <- list(nP  = 200L,
             nG  = 1000L,
             F   = 0.50,
             CR  = 0.99,
             min = c( 0,-15,-30,-30,0  ,2.5),
             max = c(15, 30, 30, 30,2.5,5  ),
             pen = penalty,
             repair = NULL,
             loopOF = TRUE,
             loopPen = FALSE,
             loopRepair = FALSE,
             printBar = FALSE,
             printDetail = FALSE,
             storeF = FALSE)

sol <- DEopt(OF = OF2, algo = algo, Data = Data)


###################################################
### code chunk number 206: NMOFman.Rnw:4211-4213
###################################################
max( abs(Data$model(sol$xbest, tm) - Data$model(betaTRUE, tm)))
sol$OFvalue


###################################################
### code chunk number 207: NMOFman.Rnw:4218-4235
###################################################
s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
system.time(sol2 <- nlminb(s0,OF2,Data = Data,
                           lower = Data$min,
                           upper = Data$max,
                           control = list(eval.max = 50000,
                           iter.max = 50000)))
max(abs(Data$model(sol2$par,tm) - Data$model(betaTRUE,tm)))
sol2$objective

par(ps = 8, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.5, 0), mar = c(4, 4, 1, 1))
plot(tm, yM, xlab = "maturities in years", ylab = "yields in %")
lines(tm,Data$model(sol$xbest,tm), col = "blue")
lines(tm,Data$model(sol2$par,tm), col = "darkgreen", lty = 2)
legend(x = "bottom", legend = c("true yields", "DE", "nlminb"),
       col = c("black", "blue", "darkgreen"),
       pch = c(1, NA, NA), lty = c(0, 1, 2))


###################################################
### code chunk number 208: NMOFman.Rnw:4240-4243
###################################################
diFa <- 1 / ((1 + NSS(sol$xbest,tm)/100)^tm)
b <- diFa %*% cfMatrix
b - bM


###################################################
### code chunk number 209: NMOFman.Rnw:4248-4252
###################################################
par(ps = 8, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.5, 0), mar = c(4, 4, 1, 1))
plot(tm, NSS(sol$xbest,tm) - NSS(betaTRUE,tm),
xlab = "maturities in years", ylab = "yield error in %")


###################################################
### code chunk number 210: NMOFman.Rnw:4259-4263
###################################################
par(ps = 8, bty = "n", las = 1, tck = 0.01,
mgp = c(3, 0.5, 0), mar = c(4, 4, 1, 1))
plot(as.numeric(unlist(lapply(tmList, max))), as.vector(b - bM),
xlab = "maturities in years", ylab = "price error in %")


###################################################
### code chunk number 211: NMOFman.Rnw:4271-4275
###################################################
beta <- c(5,-2,1,10,1,3)
yM <- NSS(beta,tm)
diFa <- 1 / ( (1 + yM/100)^tm )
b <- diFa %*% cfMatrix


###################################################
### code chunk number 212: NMOFman.Rnw:4280-4286
###################################################
B <- cbind(c(5,-2,1,10,1,3), c(4,-2,1,10,1,3))
Y <- array(0, dim = c(length(tm), ncol(B)))
for (i in 1:ncol(Y))
    Y[ ,i] <- NSS(B[ ,i], tm)
D <- 1/((1+Y/100)^tm)
t(cfMatrix) %*% D - as.vector(b)


###################################################
### code chunk number 213: NMOFman.Rnw:4297-4316
###################################################
fy <- function(ytm, cf, tm)
    sum( cf / ( (1 + ytm)^tm ) )
compYield <- function(cf, tm, guess = NULL) {
    logik <- cf != 0
    cf <- cf[logik]
    tm <- tm[logik]
    if (is.null(guess)) {ytm <- 0.05} else {ytm <- guess}
    h <- 1e-8;	dF <- 1; ci <- 0
    while (abs(dF) > 1e-5) {
        ci <- ci + 1; if (ci > 5) break
        FF  <-  fy(ytm, cf, tm)
        dFF <- (fy(ytm + h, cf, tm) - FF) / h
        dF <- FF / dFF
        ytm <- ytm - dF
    }
    if (ytm < 0)
        ytm <- 0.99
    ytm
}


###################################################
### code chunk number 214: NMOFman.Rnw:4320-4346
###################################################
OF3 <- function(param, Data) {
    tm <- Data$tm
    rM <- Data$rM
    cfMatrix<- Data$cfMatrix
    nB <- dim(cfMatrix)[2L]
    zrates <- Data$model(param,tm); aux <- 1e10
    if ( all(zrates > 0,
             !is.na(zrates))
        ) {
        diFa <- 1 / ((1 + zrates/100)^tm)
        b <- diFa %*% cfMatrix
        r <- numeric(nB)
        if ( all(!is.na(b),
                 diFa < 1,
                 diFa > 0,
                 b > 1)
            ) {
            for (bb in 1:nB) {
                r[bb] <- compYield(c(-b[bb], cfMatrix[ ,bb]), c(0,tm))
            }
            aux <- abs(r - rM)
            aux <- sum(aux)
        }
    }
    aux
}


###################################################
### code chunk number 215: NMOFman.Rnw:4356-4361
###################################################
betaTRUE <- c(5,-2,1,10,1,3)
yM <- NSS(betaTRUE, tm)
diFa <- 1 / ( (1 + yM/100)^tm )
bM <- diFa %*% cfMatrix
rM <- apply(rbind(-bM, cfMatrix), 2, compYield, c(0, tm))


###################################################
### code chunk number 216: NMOFman.Rnw:4366-4386
###################################################
Data <- list(rM = rM, tm = tm,
             cfMatrix = cfMatrix,
             model = NSS,
             min = c( 0,-15,-30,-30,0  ,2.5),
             max = c(15, 30, 30, 30,2.5,5  ),
             ww = 0.1,
             fy = fy)
algo <- list(nP = 100L,
             nG = 1000L,
             F  = 0.50,
             CR = 0.99,
             min = c( 0,-15,-30,-30,0  ,2.5),
             max = c(15, 30, 30, 30,2.5,5  ),
             pen = penalty,
             repair = NULL,
             loopOF = TRUE,
             loopPen = FALSE,
             loopRepair = FALSE,
             printBar = FALSE,
             printDetail = FALSE)


###################################################
### code chunk number 217: NMOFman.Rnw:4389-4392
###################################################
sol <- DEopt(OF = OF3, algo = algo, Data = Data)
max(abs(Data$model(sol$xbest,tm) - Data$model(betaTRUE,tm)))
sol$OFvalue


###################################################
### code chunk number 218: NMOFman.Rnw:4397-4405
###################################################
s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
sol2 <- nlminb(s0, OF3, Data = Data,
               lower = algo$min,
               upper = algo$max,
               control = list(eval.max = 50000L,
               iter.max = 50000L))
max(abs(Data$model(sol2$par,tm) - Data$model(betaTRUE,tm)))
sol2$objective


###################################################
### code chunk number 219: NMOFman.Rnw:4408-4417
###################################################
par(ps = 8, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.5, 0), mar = c(4, 4, 1, 1))
plot(tm, yM, xlab = "maturities in years", ylab = "yields in %")
lines(tm,Data$model(sol$xbest,tm), col = "blue")
lines(tm,Data$model(sol2$par,tm), col = "darkgreen", lty = 2)

legend(x = "bottom", legend = c("true yields","DE","nlminb"),
       col = c("black", "blue", "darkgreen"),
       pch = c(1, NA, NA), lty = c(0,1,2))


###################################################
### code chunk number 220: NMOFman.Rnw:4421-4423
###################################################
betaTRUE
round(sol$xbest,3)


###################################################
### code chunk number 221: NMOFman.Rnw:4436-4437
###################################################
rm(list = ls())


###################################################
### code chunk number 222: NMOFman.Rnw:4441-4443
###################################################
require("NMOF")
set.seed(94679)


###################################################
### code chunk number 223: NMOFman.Rnw:4465-4484
###################################################
randomData <- function(p = 200L,      ## number of available regressors
                        n = 200L,      ## number of observations
                        maxReg = 10L,  ## max. number of included regressors
                        s = 1,         ## standard deviation of residuals
                        constant = TRUE ) {

    X <- array(rnorm(n * p), dim = c(n, p))
    if (constant) 
        X[ ,1L] <- 1

    k <- sample.int(maxReg, 1L)   ## the number of true regressors
    K <- sort(sample.int(p, k))   ## the set of true regressors
    betatrue <- rnorm(k)          ## the true coefficients

    ## the response variable y
    y <- X[ ,K] %*% as.matrix(betatrue) + rnorm(n, sd = s)

    list(X = X, y = y, betatrue = betatrue, K = K, n = n, p = p)
}


###################################################
### code chunk number 224: NMOFman.Rnw:4488-4490
###################################################
rD <- randomData(p = 100L, n = 200L, s = 1,
                 constant = TRUE, maxReg = 10L)


###################################################
### code chunk number 225: NMOFman.Rnw:4493-4499
###################################################
Data <- list(X = rD$X,
             y = rD$y,
             n = rD$n,
             p = rD$p,
             maxk  = 30L,  ## maximum number of regressors included in model
             lognn = log(rD$n)/rD$n)


###################################################
### code chunk number 226: NMOFman.Rnw:4502-4506
###################################################
x0 <- logical(Data$p)
temp <- sample.int(Data$maxk, 1L)
temp <- sample.int(Data$p, temp)
x0[temp] <- TRUE


###################################################
### code chunk number 227: NMOFman.Rnw:4515-4516
###################################################
rD$K


###################################################
### code chunk number 228: NMOFman.Rnw:4519-4520
###################################################
which(x0)


###################################################
### code chunk number 229: NMOFman.Rnw:4533-4537
###################################################
result1 <- lm(Data$y ~ -1 + Data$X[ ,x0])
result2 <- qr.solve(Data$X[ ,x0], Data$y)
## ... coefficients should be the same
all.equal(as.numeric(coef(result1)), as.numeric(result2))


###################################################
### code chunk number 230: NMOFman.Rnw:4540-4546
###################################################
require("rbenchmark")
benchmark(lm(Data$y ~ -1 + Data$X[ ,x0]),
          qr.solve(Data$X[ ,x0], Data$y), 
         columns = c("test", "elapsed", "relative"),
          order = "test",
          replications = 1000L)


###################################################
### code chunk number 231: NMOFman.Rnw:4559-4564
###################################################
OF <- function(x, Data) {
    q <- qr(Data$X[ ,x])
    e <- qr.resid(q, Data$y)
    log(crossprod(e)/Data$n) + sum(x) * Data$lognn
}


###################################################
### code chunk number 232: NMOFman.Rnw:4567-4568
###################################################
OF(x0, Data)


###################################################
### code chunk number 233: NMOFman.Rnw:4574-4582
###################################################
neighbour <- function(xc, Data) {
    xn <- xc
    ex <- sample.int(Data$p, 1L)
    xn[ex] <- !xn[ex]
    sumx <- sum(xn)
    if ( sumx < 1L || (sumx > Data$maxk) )
        xc else xn
}


###################################################
### code chunk number 234: NMOFman.Rnw:4586-4589
###################################################
OF(neighbour(x0, Data), Data)
OF(neighbour(x0, Data), Data)
OF(neighbour(x0, Data), Data)


###################################################
### code chunk number 235: NMOFman.Rnw:4594-4601
###################################################
algo <- list(nT = 10L,    ## number of thresholds
             nS = 200L,   ## number of steps per threshold
             nD = 1000L,  ## number of random steps to compute thresholds
             neighbour = neighbour,
             x0 = x0,
             printBar = FALSE)
system.time(sol1 <- TAopt(OF, algo = algo, Data = Data))


###################################################
### code chunk number 236: NMOFman.Rnw:4607-4610
###################################################
sol1$OFvalue
which(sol1$xbest)  ## the selected regressors
rD$K               ## the true regressors


###################################################
### code chunk number 237: NMOFman.Rnw:4617-4621
###################################################
xtrue <- logical(Data$p)
xtrue[rD$K] <- TRUE
OF(sol1$xbest, Data)
OF(xtrue, Data)


###################################################
### code chunk number 238: NMOFman.Rnw:4628-4637
###################################################
restarts <- 50L
algo$printDetail <- FALSE
res <- restartOpt(TAopt, n = restarts,
                  OF = OF, algo = algo, Data = Data,
                  method = "snow", cl = 2)
par(bty = "n", las = 1,mar = c(3,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.5, 0))
plot(ecdf(sapply(res, `[[`, "OFvalue")),  ## extract solution quality
     cex = 0.4, main = "", ylab = "", xlab = "")


###################################################
### code chunk number 239: NMOFman.Rnw:4641-4649
###################################################
xbestAll <- sapply(res, `[[`, "xbest")    ## extract all solutions
inclReg  <- which(rowSums(xbestAll) > 0L) ## get included regressors
inclReg  <- sort(union(rD$K, inclReg))
data.frame(regressor   = inclReg,
           `included`  = paste(rowSums(xbestAll)[inclReg], "/", 
                                    restarts, sep = ""),
           `true regressor?` = inclReg %in% rD$K, 
           check.names = FALSE)


###################################################
### code chunk number 240: NMOFman.Rnw:4671-4708
###################################################
dim(rD$X)

neighbour2 <- function(xc, Data) {
    if ((sumx <- sum(x0)) >= Data$maxk) 
        ex <- sample(which(x0), 1L)
    else if (sumx == 1L)
        ex <- sample(which(!x0), 1L)
    else             
        ex <- sample.int(Data$p, 1L)
    xc[ex] <- !xc[ex]
    xc
}    

neighbour <- function(xc, Data) {
    xn <- xc
    ex <- sample.int(Data$p, 1L)
    xn[ex] <- !xn[ex]
    sumx <- sum(xn)
    if ( sumx < 1L || (sumx > Data$maxk) )
        xc else xn
}

algo <- list(nT = 20L,    ## number of thresholds
             nS = 200L,   ## number of steps per threshold
             nD = 1000L,  ## number of random steps to compute thresholds
             neighbour = neighbour,
             x0 = x0, q= 0.5,
             printBar = FALSE)
system.time(sol1 <- TAopt(OF, algo = algo, Data = Data))
plot(cummin(sol1$Fmat[ ,2L]), type = "l", log = "y")

##rD
OF <- function(x, Data) {
    q <- qr(Data$X[ ,x])
    e <- qr.resid(q, Data$y)
    crossprod(e)
}


###################################################
### code chunk number 241: NMOFman.Rnw:4730-4737
###################################################
bsm <- function(S, X, tau, r, q, vol, I = 1) {
    d1 <- (log(S/X) + (r - q + vol^2/2) * tau)/(vol * sqrt(tau))
    d2 <- d1 - vol * sqrt(tau)
    list(value = I * (S * exp(-q * tau) * pnorm(I * d1) -
                      X * exp(-r * tau) * pnorm(I * d2)),
         vega  = S * exp(-q*tau) * dnorm(d1 * I) * sqrt(tau))
}


###################################################
### code chunk number 242: NMOFman.Rnw:4743-4751
###################################################
S <- 99   ## sport
X <- 100  ## strike
r <- 0.01
q <- 0.0
tau <- 0.25
vol <- 0.2
I <- 1  ## a call (-1 for a put)
unlist(bsm(S, X, tau, r, q, vol, I))


###################################################
### code chunk number 243: NMOFman.Rnw:4759-4764
###################################################
tmp <- unlist(vanillaOptionEuropean(S = S, X = X, tau = tau,
                                    r = r, q = q, v = vol^2,
                                    type = 
                                    ifelse(I == 1, "call", "put")))
tmp[c("value","vega")]


###################################################
### code chunk number 244: NMOFman.Rnw:4772-4780
###################################################
S <- 99
X <- 100
r <- 0.01
q <- 0.01
tau <- 0.1
I <- 1
vol <- 0.247
(price <- bsm(S, X, tau, r, q, vol, I)$value)


###################################################
### code chunk number 245: NMOFman.Rnw:4788-4800
###################################################
impliedVol <- function(price, S, X, tau, r, q, vol0 = 0.15, I = 1,
                       tol = 1e-4, maxit = 10) {

    for (i in seq_len(maxit)) {
        tmp <- bsm(S, X, tau, r, q, vol0, I)
        step <- (tmp$value - price)/tmp$vega
        vol0 <- vol0 - step
        if (all(abs(step) < tol))
            break
    }
    vol0
}


###################################################
### code chunk number 246: NMOFman.Rnw:4807-4808
###################################################
impliedVol(price, S, X, tau, r, q, vol, I)


###################################################
### code chunk number 247: NMOFman.Rnw:4813-4815
###################################################
vanillaOptionImpliedVol(exercise = "european", price, S, X, tau, r,
                          q, type = "call")


###################################################
### code chunk number 248: NMOFman.Rnw:4822-4829 (eval = FALSE)
###################################################
## benchmark(iV = impliedVol(price, S, X, tau, r, q, runif(1L) + 0.01,I),
##           vanOptIV = vanillaOptionImpliedVol(exercise = "european",
##           price, S, X, tau, r,
##           q, tauD = 0, D = 0, type = "call",
##           M = 101, uniroot.info = FALSE),
##           columns = c("test", "elapsed", "relative"),
##           replications = 1e2, order = "relative")


###################################################
### code chunk number 249: NMOFman.Rnw:4839-4844
###################################################
S <- rep(99, 21)  ## spot
X <- 90:110      ## strike
r <- 0.01; q <- 0.02
tau <- 0.2; vol <- 0.24; I <- 1
data.frame(S = S, X = X, bsm(S, X, tau, r, q, vol, I))


###################################################
### code chunk number 250: NMOFman.Rnw:4850-4853
###################################################
Xvec <- 80:120
tauvec <- c(c(3, 6, 9)/12,  ## 3, 6, 9 months
            1, 2, 3, 4, 5)  ## 1..5 years


###################################################
### code chunk number 251: NMOFman.Rnw:4859-4866
###################################################
loop <- function() {
    callprices <- array(NA, dim = c(length(Xvec), length(tauvec)))
    for (X in Xvec)
        for (tau in tauvec)
            callprices[X == Xvec, tau == tauvec] <- bsm(S,X,tau,r,q,vol)$value
    callprices
}


###################################################
### code chunk number 252: NMOFman.Rnw:4871-4877
###################################################
vect <- function() {
    tmp <- expand.grid(Xvec,tauvec)
    callprices <- bsm(S, tmp[[1L]], tmp[[2L]], r, q, vol, I)$value
    dim(callprices) <- c(length(Xvec), length(tauvec))
    callprices
}


###################################################
### code chunk number 253: NMOFman.Rnw:4882-4893
###################################################
S <- 101
Xvec <- 80:120
tauvec <- c(c(3, 6, 9)/12,  ## 3, 6, 9 months
            1, 2, 3, 4, 5)  ## 1..5 years
  r <- 0.01; q <- 0.01
  tau <- 0.25; vol <- 0.2; I <- 1

  callprices1 <- loop()
  callprices2 <- vect()

  all.equal(callprices1, callprices2)


###################################################
### code chunk number 254: NMOFman.Rnw:4899-4902 (eval = FALSE)
###################################################
## benchmark(loop(), vect(),
##           columns = c("test", "elapsed", "relative"),
## replications = 1e2, order = "relative")


###################################################
### code chunk number 255: NMOFman.Rnw:4909-4920
###################################################
S <- rep(99,21)  ## spot
X <- 90:110
r <- 0.01
q <- 0.02
tau <- runif(21)
vol <- (runif(21)+0.2)/3
ivol <- impliedVol(bsm(S, X, tau, r, q, vol, I)$value,
S, X, tau, r, q, vol = 0.2,
I, tol = 1e-5, maxit = 5)

data.frame(S = S, X = X, vol = vol, ivol = ivol, diff = abs(vol-ivol))


###################################################
### code chunk number 256: NMOFman.Rnw:4927-4931 (eval = FALSE)
###################################################
## benchmark(iv = impliedVol(bsm(S, X, tau, r, q, vol, I)$value,
##           S, X, tau, r, q, tol = 1e-5, maxit = 5),
##           columns = c("test", "elapsed"),
##           replications = 1e3, order = "relative")


###################################################
### code chunk number 257: NMOFman.Rnw:5027-5053 (eval = FALSE)
###################################################
##   valDate <- as.Date("2012-02-10")
## 
##   calls <- optionData$pricesCall
##   puts  <- optionData$pricesPut
##   spot  <- optionData$index
## 
##   ## strikes
##   strikes <- as.numeric(row.names(calls))
##   
##   ## compute interest rates
##   expDates <- as.Date(paste0(colnames(calls),"15"),
##                       format = "%Y%m%d")
##   tau <- as.numeric((expDates - valDate)/365)
##   r <- NSS(optionData$NSSpar, tau)
## 
##   ## compute basis
##   notNA <- !is.na(calls) & !is.na(puts)
##   i <- 2
##   q <- numeric(ncol(calls))
##   for (i in seq_along(q))
##       q[i] <- mean(-log((calls[,i] - puts[,i] +
##                          exp(-r[i]*tau[i]) * strikes)/spot)/tau[i],
##                    na.rm = TRUE)
## 
##   ## check with future
##   ((log(spot)-log(optionData$future[2]))+r[2]*tau[2])/tau[2]


###################################################
### code chunk number 258: NMOFman.Rnw:5057-5071 (eval = FALSE)
###################################################
##   pricevec <- as.numeric(calls)
##   notNA <- !is.na(pricevec)
##   Xvec   <- rep(strikes, ncol(calls))
##   tauvec <- rep(tau, each = nrow(calls))
##   rvec   <- rep(r,   each = nrow(calls))
##   qvec   <- rep(q,   each = nrow(calls))
##   ivol <- impliedVol(pricevec[notNA],
##                      spot, Xvec[notNA], tauvec[notNA],
##                      rvec[notNA], qvec[notNA], vol = 0.25, maxit = 5)
## 
##   vols <- numeric(length(pricevec))
##   vols[notNA] <- ivol
##   data.frame(pricevec[notNA],
##              spot, Xvec[notNA], tauvec[notNA], vols=ivol)


###################################################
### code chunk number 259: NMOFman.Rnw:5080-5099
###################################################
S <- 100    ## spot
X <- 100    ## strike
tau <- 1    ## time-to-maturity
r <- 0.02   ## interest rate
q <- 0.02   ## dividend rate
v <- 0.2    ## volatility

## the closed-form solution
callBSM <- function(S,X,tau,r,q,v) {
    d1 <- (log(S/X) + (r - q + v^2 / 2)*tau) / (v*sqrt(tau))
    d2 <- d1 - v*sqrt(tau)
    S * exp(-q * tau) * pnorm(d1) -  X * exp(-r * tau) * pnorm(d2)
}
callBSM(S,X,tau,r,q,v)

## with the characteristic function
callCF(cf = cfBSM, S = S, X = X, tau = tau, r = r, q = q,
       v = v^2,  ## variance, not vol
       implVol = TRUE)


###################################################
### code chunk number 260: NMOFman.Rnw:5105-5108
###################################################
Xvec <- 80:120
tauvec <- c(c(3, 6, 9)/12,  ## 3, 6, 9 months
            1, 2, 3, 4, 5)  ## 1..5 years


###################################################
### code chunk number 261: NMOFman.Rnw:5112-5119
###################################################
fun1 <- function() {
    callprices <- array(NA, dim = c(length(Xvec), length(tauvec)))
    for (X in Xvec)
        for (tau in tauvec)
            callprices[X == Xvec, tau == tauvec] <- callBSM(S,X,tau,r,q,v)
    callprices
}


###################################################
### code chunk number 262: NMOFman.Rnw:5123-5129
###################################################
fun2 <- function() {
    tmp <- expand.grid(Xvec,tauvec)
    callprices <- callBSM(S, tmp[[1]], tmp[[2]], r, q, v)
    dim(callprices) <- c(length(Xvec), length(tauvec))
    callprices
}


###################################################
### code chunk number 263: NMOFman.Rnw:5132-5136
###################################################
callprices1 <- fun1()
callprices2 <- fun2()

all.equal(callprices1, callprices2)


###################################################
### code chunk number 264: NMOFman.Rnw:5140-5146
###################################################
system.time(
    for (i in 1:100)
        ignore <- fun1() )
system.time(
    for (i in 1:100)
        ignore <- fun2() )


###################################################
### code chunk number 265: NMOFman.Rnw:5151-5199
###################################################
priceMatrix <- function(cf, S, Xvec, tauvec, r, q = 0, ...,
                        nodes = NULL, weights = NULL, n = 200) {

    if (is.null(nodes)) {
        tmp <- xwGauss(n)
        tmp <- changeInterval(tmp$nodes, tmp$weights,
                              oldmin = -1, oldmax = 1,
                              newmin =  0, newmax = 200)
        nodes <- tmp$nodes
        weights <- tmp$weights
    }

    callprices <- array(NA, dim = c(length(Xvec), length(tauvec)))
    tmpmat <- array(NA, dim = c(length(Xvec), length(weights)))

    inodes <- 1i * nodes
    itau <- 0L
    for (tau in tauvec) {
        itau <- itau + 1L
        cfi <- S * exp((r - q) * tau)
        cf1 <- cf(nodes - 1i, S, tau, r, q, ...)/inodes/cfi
        cf2 <- cf(nodes,      S, tau, r, q, ...)/inodes
        iX <- 0L
        for (X in Xvec) {
            iX <- iX + 1L
            if (itau == 1L)
                tmpmat[iX, ] <- exp(-inodes * log(X))
            P1 <- 0.5 + weights %*% Re(tmpmat[iX, ] * cf1)/pi
            P2 <- 0.5 + weights %*% Re(tmpmat[iX, ] * cf2)/pi
            callprices[iX, itau] <-
                exp(-q*tau) * S * P1 - exp(-r*tau) * X * P2
        }
   }
    callprices
}
v <- .2

system.time(
    for (i in 1:100)
        ignore <- priceMatrix(cf = cfBSM, S, Xvec, tauvec, r, q = q,
                   v = 0.2^2, n = 50) )
require("compiler")
priceMatrix2 <- cmpfun(priceMatrix)
system.time(
    for (i in 1:100)
        ignore <- priceMatrix2(cf = cfBSM, S, Xvec, tauvec, r, q = q,
                   v = 0.2^2, n = 50) )



###################################################
### code chunk number 266: NMOFman.Rnw:5202-5232
###################################################
cfp <- priceMatrix(cf = cfBSM, S, Xvec, tauvec, r, q = q,
                   v = 0.2^2, n = 100)

callprices1[1:5, 1:5]
callprices2[1:5, 1:5]
cfp[1:5, 1:5]
all.equal(callprices1, cfp)
system.time(
    for (i in 1:100)
        ignore <- fun1() )
system.time(
    for (i in 1:100)
        ignore <- fun2() )
system.time(
    for (i in 1:100)
        ignore <- priceMatrix(cf = cfBSM, S, Xvec, tauvec, r, q = q,
                   v = 0.2^2, n = 50) )

tmp <- xwGauss(50)
tmp <- changeInterval(tmp$nodes, tmp$weights,
                      oldmin = -1, oldmax = 1,
                      newmin =  0, newmax = 200)
nodes <- tmp$nodes
weights <- tmp$weights
system.time(
    for (i in 1:100)
        ignore <- priceMatrix(cf = cfBSM, S, Xvec, tauvec, r, q = q,
                              v = 0.2^2,
                              nodes = nodes, weights = weights) )



###################################################
### code chunk number 267: NMOFman.Rnw:5322-5324 (eval = FALSE)
###################################################
## install.packages("NMOF") ## CRAN
## install.packages("NMOF", repos = "http://enricoschumann.net/R")


###################################################
### code chunk number 268: NMOFman.Rnw:5332-5334 (eval = FALSE)
###################################################
## require("NMOF")
## showExample("exampleOF.R")


###################################################
### code chunk number 269: NMOFman.Rnw:5354-5355
###################################################
toLatex(sessionInfo())


