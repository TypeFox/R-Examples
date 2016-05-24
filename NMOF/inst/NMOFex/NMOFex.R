### R code from vignette source 'NMOFex.Rnw'

###################################################
### code chunk number 1: NMOFex.Rnw:74-76
###################################################
## version 2012-01-15
options(continue = " ", digits = 3, width = 70)


###################################################
### code chunk number 2: NMOFex.Rnw:94-96 (eval = FALSE)
###################################################
## install.packages("NMOF") ### CRAN
## install.packages("NMOF", repos = "http://R-Forge.R-project.org")


###################################################
### code chunk number 3: NMOFex.Rnw:99-101
###################################################
require("NMOF")
set.seed(1112233344)


###################################################
### code chunk number 4: NMOFex.Rnw:108-110 (eval = FALSE)
###################################################
## whereToLook <- system.file("NMOFex/NMOFex.R", package = "NMOF")
## file.show(whereToLook, title = "NMOF examples")


###################################################
### code chunk number 5: NMOFex.Rnw:141-145
###################################################
testFun  <- function(x) {
    Sys.sleep(0.1) ## just wasting time :-)
    x[1L] + x[2L]^2
}


###################################################
### code chunk number 6: NMOFex.Rnw:149-153
###################################################
lower <- c(1,3); upper <- 5; n <- 5L
system.time(sol1 <- gridSearch(fun = testFun,
                               lower = lower, upper = upper,
                               n = n, printDetail = TRUE))


###################################################
### code chunk number 7: NMOFex.Rnw:157-159
###################################################
seq(from = 1, to = 5, length.out= n)  ## x_1
seq(from = 3, to = 5, length.out= n)  ## x_2


###################################################
### code chunk number 8: NMOFex.Rnw:163-165
###################################################
sol1$minfun
sol1$minlevels


###################################################
### code chunk number 9: NMOFex.Rnw:169-175
###################################################
system.time(sol2 <- gridSearch(fun = testFun,
                               lower = lower, upper = upper,
                               n = n, printDetail = FALSE,
                               method = "snow",  ### use 'snow' ...
                               cl = 2L))         ### ... with 2 cores
all.equal(sol1, sol2)


###################################################
### code chunk number 10: NMOFex.Rnw:191-198
###################################################
## create random data with vols between minVol and maxVol
## and pairwise correlation of 0.6
na <- 500L
C <- array(0.6, dim = c(na,na)); diag(C) <- 1
minVol <- 0.20; maxVol <- 0.40
Vols <- (maxVol - minVol) * runif(na) + minVol
Sigma <- outer(Vols,Vols) * C


###################################################
### code chunk number 11: NMOFex.Rnw:210-228
###################################################
## objective function for LS/TA
OF <- function(x, data) {
    sx <- sum(x)
    w <- rep.int(1/sx, sx)
    res <- crossprod(w, data$Sigma[x, x])
    tcrossprod(w, res)
}

## neighbourhood function for LS/TA
neighbour <- function(xc, data) {
    xn <- xc
    p <- sample.int(data$na, data$nn, replace = FALSE)
    xn[p] <- !xn[p]
    ## reject infeasible solution
    sumx <- sum(xn)
    if ( (sumx > data$Ksup) || (sumx < data$Kinf) )
        xc else xn
}


###################################################
### code chunk number 12: NMOFex.Rnw:235-241
###################################################
## data
data <- list(Sigma = Sigma,   ### cov-matrix
              Kinf = 30L,      ### min cardinality
              Ksup = 60L,      ### max cardinality
              na = na,         ### number of assets
              nn = 1L)         ### how many assets to change per iteration


###################################################
### code chunk number 13: NMOFex.Rnw:245-250
###################################################
## a random solution x0
card0 <- sample(data$Kinf:data$Ksup, 1L, replace = FALSE)
assets <- sample.int(na, card0, replace = FALSE)
x0 <- logical(na)
x0[assets] <- TRUE


###################################################
### code chunk number 14: NMOFex.Rnw:256-264
###################################################
## *Local Search*
algo <- list(x0 = x0, neighbour = neighbour, nS = 5000L,
             printDetail = FALSE, printBar = FALSE)
system.time(solLS <- LSopt(OF, algo = algo, data = data))

## *Threshold Accepting*
algo$nT <- 10L; algo$nS <- trunc(algo$nS/algo$nT); algo$q <- 0.2
system.time(solTA <- TAopt(OF, algo = algo, data = data))


###################################################
### code chunk number 15: NMOFex.Rnw:346-354
###################################################
## *Genetic Algorithm*
OF2 <- function(x, data) {
    res <- colSums(data$Sigma %*% x * x)
    n <- colSums(x); res <- res / n^2
    ## penalise
    p <- pmax(data$Kinf - n, 0) + pmax(n - data$Ksup, 0)
    res + p
}


###################################################
### code chunk number 16: NMOFex.Rnw:360-363
###################################################
algo <- list(nB = na, nP = 100L, nG = 500L, prob = 0.002,
             printBar = FALSE, loopOF = FALSE)
system.time(solGA <- GAopt(OF = OF2, algo = algo, data = data))


###################################################
### code chunk number 17: NMOFex.Rnw:366-371
###################################################
cat(
    "Local Search        ", format(sqrt(solLS$OFvalue), digits = 4), "\n",
    "Threshold Accepting ", format(sqrt(solTA$OFvalue), digits = 4), "\n",
    "Genetic Algorithm   ", format(sqrt(solGA$OFvalue), digits = 4), "\n",
    sep = "")


###################################################
### code chunk number 18: NMOFex.Rnw:395-403
###################################################
na <- 100L                                 ### number of assets
ns <- 200L                                 ### number of scenarios
vols <- runif(na, min = 0.2, max = 0.4)    ### marginal vols
C <- matrix(0.6, na, na); diag(C) <- 1     ### correlation matrix
R <- rnorm(ns * na)/16                     ### random returns
dim(R) <- c(ns, na)
R <- R %*% chol(C)
R <- R %*% diag(vols)


###################################################
### code chunk number 19: NMOFex.Rnw:424-435
###################################################
data <- list(R = t(R),              ### scenarios
              theta = 0.005,         ### return threshold
              na = na,               ### number of assets
              ns = ns,               ### number of scenarios
              max = rep( 0.05, na),  ### DE: vector of max. weight
              min = rep(-0.05, na),  ### DE: vector of min. weight
              wsup =  0.05,          ### TA: max weight
              winf = -0.05,          ### TA: min weight
              eps = 0.5/100,         ### TA: step size
              w = 1                  ### penalty weight
             )


###################################################
### code chunk number 20: NMOFex.Rnw:440-443
###################################################
x0 <- data$min + runif(data$na)*(data$max - data$min)
x0[1:5]
sum(x0)


###################################################
### code chunk number 21: NMOFex.Rnw:448-452
###################################################
temp <- R %*% x0             ### compute portfolio returns
temp <- temp - data$theta
temp <- (temp[temp < 0])^2
sum(temp)/ns                 ### semivariance


###################################################
### code chunk number 22: NMOFex.Rnw:457-464
###################################################
OF <- function(x, data) {
    Rx <- crossprod(data$R, x)
    Rx <- Rx - data$theta
    Rx <- Rx - abs(Rx)
    Rx <- Rx * Rx
    colSums(Rx) /(4*data$ns)
}


###################################################
### code chunk number 23: NMOFex.Rnw:473-475
###################################################
OF(x0, data)
OF(cbind(x0, x0), data)


###################################################
### code chunk number 24: NMOFex.Rnw:482-492
###################################################
repair <- function(x, data) {
    myFun <- function(x) x/sum(x)
    if (is.null(dim(x)[2L]))
        myFun(x) else apply(x, 2L, myFun)
}
repair2 <- function(x, data) {
    myFun <- function(x) x + (1 - sum(x))/data$na
    if (is.null(dim(x)[2L]))
        myFun(x) else apply(x, 2L, myFun)
}


###################################################
### code chunk number 25: NMOFex.Rnw:496-501
###################################################
sum(x0)
sum(repair(x0, data))

colSums(repair( cbind(x0, x0), data))
colSums(repair2(cbind(x0, x0), data))


###################################################
### code chunk number 26: NMOFex.Rnw:504-515
###################################################
penalty <- function(x, data) {
    up <- data$max
    lo <- data$min
    xadjU <- x - up
    xadjU <- xadjU + abs(xadjU)
    xadjL <- lo - x
    xadjL <- xadjL + abs(xadjL)
    if (is.null(dim(x)[2L]))
        data$w * (sum(xadjU) + sum(xadjL)) else
        data$w * (colSums(xadjU) + colSums(xadjL))
}


###################################################
### code chunk number 27: NMOFex.Rnw:521-527
###################################################
x0[1L] <- 0.30
penalty(x0, data)
penalty(cbind(x0, x0), data)
x0[1L] <- 0
penalty(x0, data)
penalty(cbind(x0, x0), data)


###################################################
### code chunk number 28: NMOFex.Rnw:531-541
###################################################
algo <- list(nP = 100,        ### population size
             nG = 1000,       ### number of generations
             F = 0.25,        ### step size
             CR = 0.9,
             min = data$min,
             max = data$max,
             repair = repair,
             pen = penalty,
             printBar = FALSE, printDetail = TRUE,
             loopOF = TRUE, loopPen = TRUE, loopRepair = TRUE)


###################################################
### code chunk number 29: NMOFex.Rnw:546-553
###################################################
system.time(sol <- DEopt(OF = OF,algo = algo,data = data))
16 * 100 * sqrt(sol$OFvalue)   ### solution quality

## check constraints
all(all.equal(sum(sol$xbest), 1),  ### budget constraint
    sol$xbest <= data$max,         ### holding size constraints
    sol$xbest >= data$min)


###################################################
### code chunk number 30: NMOFex.Rnw:562-569
###################################################
## looping over the population
algo$loopOF <- TRUE; algo$loopPen <- TRUE; algo$loopRepair <- TRUE
system.time(sol <- DEopt(OF = OF,algo = algo, data = data))

## evaluating the population in one step
algo$loopOF <- FALSE; algo$loopPen <- FALSE; algo$loopRepair <- FALSE
system.time(sol <- DEopt(OF = OF,algo = algo, data = data))


###################################################
### code chunk number 31: NMOFex.Rnw:579-592
###################################################
algo$printDetail <- FALSE
restartsDE <- restartOpt(fun = DEopt,      ### what function
                         n = 20L,          ### how many restarts
                         OF = OF,
                         algo = algo,
                         data = data,
                         method = "snow",  ### using package snow
                         cl = 2)           ### 2 cores

## extract best solution
OFvaluesDE <- sapply(restartsDE, `[[`, "OFvalue")
OFvaluesDE <- 16 * 100 * sqrt(OFvaluesDE)
weightsDE <- sapply(restartsDE, `[[`, "xbest")


###################################################
### code chunk number 32: NMOFex.Rnw:595-600
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
plot(sort(OFvaluesDE), (seq_len(length(OFvaluesDE))) / length(OFvaluesDE),
     type = "S", ylim = c(0, 1), xlab = "", ylab = "")
mtext("OF value",  side = 1, line = 1)


###################################################
### code chunk number 33: NMOFex.Rnw:604-610
###################################################
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
boxplot(t(weightsDE),
        outline = FALSE, boxwex = 0.4, ylim = c(-0.06,0.06))
mtext("assets",  side = 1, line = 1)
mtext("weights", side = 2, line = 1.3, las = 1, padj = -5)


###################################################
### code chunk number 34: NMOFex.Rnw:622-660
###################################################
algo$printDetail <- FALSE;  algo$nP <- 200L; restarts <- 20L
nGs <- c(500L, 1000L, 2500L)
lstOFvaluesDE <- list()
for (i in 1:3) {
    algo$nG <- nGs[i]
    restartsDE <- restartOpt(fun = DEopt,
                             n = restarts,
                             OF = OF,  algo = algo, data = data,
                             method = "snow", cl = 2)
    OFvaluesDE <- sapply(restartsDE, `[[`, "OFvalue") ### extract best solution
    OFvaluesDE <- 16 * 100 * sqrt(OFvaluesDE)
    lstOFvaluesDE[[i]] <- OFvaluesDE
}
res <- simplify2array(lstOFvaluesDE)

## now with 'repair2'
algo$repair <- repair2
lstOFvaluesDE <- list()
for (i in 1:3) {
    algo$nG <- nGs[i]
    restartsDE <- restartOpt(fun = DEopt,
                             n = restarts,
                             OF = OF,  algo = algo, data = data,
                             method = "snow", cl = 2)
    OFvaluesDE <- sapply(restartsDE, `[[`, "OFvalue") ### extract best solution
    OFvaluesDE <- 16 * 100 * sqrt(OFvaluesDE)
    lstOFvaluesDE[[i]] <- OFvaluesDE
}
res2 <- simplify2array(lstOFvaluesDE)

## plot results
allres <- as.vector(rbind(res,res2))
xlims <- pretty(allres); xlims <- c(min(xlims), max(xlims))
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
plot(ecdf(res[ ,3L]), xlim = xlims, cex = 0.4, main = "", ylab = "", xlab = "")
for (i in 1:2) lines(ecdf(res[  ,i]), cex = 0.4)
for (i in 1:3) lines(ecdf(res2[ ,i]), col = "blue", cex = 0.4)


###################################################
### code chunk number 35: NMOFex.Rnw:669-676
###################################################
weightsDE <- sapply(restartsDE, `[[`, "xbest")
par(bty = "n", las = 1, mar = c(3, 4, 0, 0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
boxplot(t(weightsDE),
        outline = FALSE, boxwex = 0.4, ylim = c(-0.06, 0.06))
mtext("assets",  side = 1, line = 1)
mtext("weights", side = 2, line = 1.3, las = 1, padj = -5)


###################################################
### code chunk number 36: NMOFex.Rnw:682-699
###################################################
algo <- list(nP = 100L,        ### population size
    nG = 1000L,                ### number of generations
    c1 = 0.5,                  ### weight for individually best solution
    c2 = 1.5,                  ### weight for overall best solution
    min = data$min,
    max = data$max,
    repair = repair, pen = penalty,
    iner = 0.7, initV = 1, maxV = 0.2,
    printBar = FALSE, printDetail = TRUE)

system.time(sol <- PSopt(OF = OF,algo = algo,data = data))
16 * 100 * sqrt(sol$OFvalue)      ### solution quality

## check constraints
all(all.equal(sum(sol$xbest),1),  ### budget constraint
    sol$xbest <= data$max,
    sol$xbest >= data$min)


###################################################
### code chunk number 37: NMOFex.Rnw:705-737
###################################################
## adjusting velocity
changeV <- function(x, data) {
    myFun <- function(x) x - (sum(x))/data$na
    if (is.null(dim(x)[2L]))
        myFun(x) else apply(x, 2L, myFun)
}
sum(changeV(x0, data))
colSums(changeV(cbind(x0, x0), data))

## initial population that meets budget constraint
initP <- data$min + diag(data$max - data$min) %*%
         array(runif(length(data$min) * algo$nP),
               dim = c(length(data$min),  algo$nP))
colSums(initP <- repair(initP,data))[1:10]

## add to 'algo'
algo$changeV <- changeV        ### function to adjust velocity
algo$initP <- initP            ### initial population
algo$repair <- NULL            ### not needed anymore

system.time(sol <- PSopt(OF = OF,algo = algo, data = data))
16 * 100 * sqrt(sol$OFvalue)   ### solution quality

## check constraints
all(all.equal(sum(sol$xbest), 1), ### budget constraint
    sol$xbest <= data$max,
    sol$xbest >= data$min)

## vectorised
algo$loopOF <- FALSE; algo$loopPen <- FALSE
algo$loopRepair <- FALSE; algo$loopChangeV <- FALSE
system.time(sol <- PSopt(OF = OF, algo = algo, data = data))


###################################################
### code chunk number 38: NMOFex.Rnw:740-755
###################################################
algo$printDetail <- FALSE
restartsPS <- restartOpt(fun = PSopt,
                         n = 20L,
                         OF = OF,
                         algo = algo, data = data,
                         method = "snow", cl = 2)

## extract best solution
OFvaluesPS <- sapply(restartsPS, `[[`, "OFvalue")
OFvaluesPS <- 16 * 100 * sqrt(OFvaluesPS)
par(bty = "n", las = 1,mar = c(3,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
plot(sort(OFvaluesPS), (seq_len(length(OFvaluesPS))) / length(OFvaluesPS),
     type = "S", ylim = c(0, 1), xlab = "", ylab = "")
mtext("OF value",  side = 1, line = 1)


###################################################
### code chunk number 39: NMOFex.Rnw:763-786
###################################################
data$R <- R  ## not transposed any more

neighbourU <- function(sol, data){
    resample <- function(x, ...)
        x[sample.int(length(x), ...)]
    wn <- sol$w
    toSell <- wn > data$winf
    toBuy  <- wn < data$wsup
    i <- resample(which(toSell), size = 1L)
    j <- resample(which(toBuy), size = 1L)
    eps <- runif(1) * data$eps
    eps <- min(wn[i] - data$winf, data$wsup - wn[j], eps)
    wn[i] <- wn[i] - eps
    wn[j] <- wn[j] + eps
    Rw <- sol$Rw + data$R[,c(i,j)] %*% c(-eps,eps)
    list(w = wn, Rw = Rw)
}
OF <- function(x, data) {
    Rw <- x$Rw - data$theta
    Rw <- Rw - abs(Rw)
    sum(Rw*Rw) / (4*data$ns)
}



###################################################
### code chunk number 40: NMOFex.Rnw:789-802
###################################################
## a random initial weights
w0 <- runif(data$na); w0 <- w0/sum(w0)
x0 <- list(w = w0, Rw = R %*% w0)
algo <- list(x0 = x0,
             neighbour = neighbourU,
             nS = 2000L,
             nT = 10L,
             nD = 5000L,
             q = 0.20,
             printBar = FALSE,
             printDetail = FALSE)
system.time(sol2 <- TAopt(OF,algo,data))
16 * 100 * sqrt(sol2$OFvalue)


###################################################
### code chunk number 41: NMOFex.Rnw:806-826
###################################################
restartsTA <- restartOpt(fun = TAopt,
                         n = 20L,
                         OF = OF,
                         algo = algo, data = data,
                         method = "snow", cl = 2)

OFvaluesTA <- sapply(restartsTA, `[[`, "OFvalue") # extract best solution
OFvaluesTA <- 16 * 100 * sqrt(OFvaluesTA)
weightsTA <- sapply(restartsTA, `[[`, "xbest")
par(bty = "n", las = 1,mar = c(3,4,0,0), ps = 8,
    tck = 0.001, mgp = c(3, 0.2, 0))

## blue: DE solution with nP = 200 and nG = 2000
xlims <- pretty(c(res2[,3], OFvaluesTA))
plot(ecdf(res2[,3]), col = "blue", cex = 0.4,
     main = "", ylab = "", xlab = "",
xlim = c(min(xlims), max(xlims)) )

## black: TA
lines(ecdf(OFvaluesTA), cex = 0.4)


###################################################
### code chunk number 42: NMOFex.Rnw:843-863
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
### code chunk number 43: NMOFex.Rnw:869-884
###################################################
cfList <- list(cf1,cf2,cf3,cf4,cf5,cf6,cf7,cf8,cf9,cf10,
               cf11,cf12,cf13,cf14,cf15,cf16,cf17,cf18,cf19,cf20)
tmList <- list(tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,
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
### code chunk number 44: NMOFex.Rnw:891-895
###################################################
betaTRUE <- c(5,-2,1,10,1,3)
yM <- NSS(betaTRUE,tm)
diFa <- 1 / ( (1 + yM/100)^tm )
bM <- diFa %*% cfMatrix


###################################################
### code chunk number 45: NMOFex.Rnw:899-903
###################################################
data <- list(bM = bM, tm = tm, cfMatrix = cfMatrix, model = NSS,
             ww = 1,
            min = c( 0,-15,-30,-30,0  ,2.5),
            max = c(15, 30, 30, 30,2.5,5  ))


###################################################
### code chunk number 46: NMOFex.Rnw:913-923
###################################################
OF2 <- function(param, data) {
    tm <- data$tm
    bM <- data$bM
    cfMatrix <- data$cfMatrix
    diFa  <- 1 / ((1 + data$model(param, tm)/100)^tm)
    b <- diFa %*% cfMatrix
    aux <- b - bM; aux <- max(abs(aux))
    if (is.na(aux)) aux <- 1e10
    aux
}


###################################################
### code chunk number 47: NMOFex.Rnw:927-942
###################################################
penalty <- function(mP, data) {
    minV <- data$min
    maxV <- data$max
    ww <- data$ww
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
### code chunk number 48: NMOFex.Rnw:947-963
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

sol <- DEopt(OF = OF2, algo = algo, data = data)


###################################################
### code chunk number 49: NMOFex.Rnw:969-971
###################################################
max( abs(data$model(sol$xbest, tm) - data$model(betaTRUE, tm)))
sol$OFvalue


###################################################
### code chunk number 50: NMOFex.Rnw:975-992
###################################################
s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
system.time(sol2 <- nlminb(s0,OF2,data = data,
                                 lower = data$min,
                                 upper = data$max,
                               control = list(eval.max = 50000,
                                              iter.max = 50000)))
max(abs(data$model(sol2$par,tm) - data$model(betaTRUE,tm)))
sol2$objective

par(ps = 11, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.2, 0), mar = c(4, 4, 1, 1))
plot(tm, yM, xlab = "maturities in years", ylab = "yields in %")
lines(tm,data$model(sol$xbest,tm), col = "blue")
lines(tm,data$model(sol2$par,tm), col = "darkgreen", lty = 2)
legend(x = "bottom", legend = c("true yields", "DE", "nlminb"),
       col = c("black", "blue", "darkgreen"),
       pch = c(1, NA, NA), lty = c(0, 1, 2))


###################################################
### code chunk number 51: NMOFex.Rnw:996-999
###################################################
diFa <- 1 / ((1 + NSS(sol$xbest,tm)/100)^tm)
b <- diFa %*% cfMatrix
b - bM


###################################################
### code chunk number 52: NMOFex.Rnw:1003-1007
###################################################
par(ps = 11, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.2, 0), mar = c(4, 4, 1, 1))
plot(tm, NSS(sol$xbest,tm) - NSS(betaTRUE,tm),
     xlab = "maturities in years", ylab = "yield error in %")


###################################################
### code chunk number 53: NMOFex.Rnw:1013-1017
###################################################
par(ps = 11, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.2, 0), mar = c(4, 4, 1, 1))
plot(as.numeric(unlist(lapply(tmList, max))), as.vector(b - bM),
     xlab = "maturities in years", ylab = "price error in %")


###################################################
### code chunk number 54: NMOFex.Rnw:1028-1047
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
### code chunk number 55: NMOFex.Rnw:1051-1077
###################################################
OF3 <- function(param, data) {
    tm <- data$tm
    rM <- data$rM
    cfMatrix<- data$cfMatrix
    nB <- dim(cfMatrix)[2L]
    zrates <- data$model(param,tm); aux <- 1e10
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
### code chunk number 56: NMOFex.Rnw:1086-1091
###################################################
betaTRUE <- c(5,-2,1,10,1,3)
yM <- NSS(betaTRUE, tm)
diFa <- 1 / ( (1 + yM/100)^tm )
bM <- diFa %*% cfMatrix
rM <- apply(rbind(-bM, cfMatrix), 2, compYield, c(0, tm))


###################################################
### code chunk number 57: NMOFex.Rnw:1095-1115
###################################################
data <- list(rM = rM, tm = tm,
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
### code chunk number 58: NMOFex.Rnw:1118-1121
###################################################
sol <- DEopt(OF = OF3, algo = algo, data = data)
max(abs(data$model(sol$xbest,tm) - data$model(betaTRUE,tm)))
sol$OFvalue


###################################################
### code chunk number 59: NMOFex.Rnw:1125-1133
###################################################
s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
sol2 <- nlminb(s0, OF3, data = data,
                                   lower = algo$min,
                                   upper = algo$max,
                                 control = list(eval.max = 50000L,
                                                iter.max = 50000L))
max(abs(data$model(sol2$par,tm) - data$model(betaTRUE,tm)))
sol2$objective


###################################################
### code chunk number 60: NMOFex.Rnw:1136-1145
###################################################
par(ps = 11, bty = "n", las = 1, tck = 0.01,
    mgp = c(3, 0.2, 0), mar = c(4, 4, 1, 1))
plot(tm, yM, xlab = "maturities in years", ylab = "yields in %")
lines(tm,data$model(sol$xbest,tm), col = "blue")
lines(tm,data$model(sol2$par,tm), col = "darkgreen", lty = 2)

legend(x = "bottom", legend = c("true yields","DE","nlminb"),
       col = c("black", "blue", "darkgreen"),
       pch = c(1, NA, NA), lty = c(0,1,2))


###################################################
### code chunk number 61: NMOFex.Rnw:1149-1151
###################################################
betaTRUE
round(sol$xbest,3)


###################################################
### code chunk number 62: NMOFex.Rnw:1180-1209
###################################################
randomData <- function(
    p = 200L,      ### number of available regressors
    n = 200L,      ### number of observations
    maxReg = 10L,  ### max. number of included regressors
    s = 1,         ### standard deviation of residuals
    constant = TRUE ) {

    X <- rnorm(n * p); dim(X) <- c(n, p)  # regressor matrix X
    if (constant) X[ ,1L] <- 1

    k <- sample.int(maxReg, 1L)   ### the number of true regressors
    K <- sort(sample.int(p, k))   ### the set of true regressors
    betatrue <- rnorm(k)          ### the true coefficients

    ## the response variable y
    y <- X[ ,K] %*% as.matrix(betatrue) + rnorm(n, sd = s)

    list(X = X, y = y, betatrue = betatrue, K = K, n = n, p = p)
}

rD <- randomData(p = 100L, n = 200L, s = 1,
                 constant = TRUE, maxReg = 15L)

data <- list(X = rD$X,
             y = rD$y,
             n = rD$n,
             p = rD$p,
             maxk  = 30L,  ### maximum number of regressors included in model
             lognn = log(rD$n)/rD$n)


###################################################
### code chunk number 63: NMOFex.Rnw:1215-1219
###################################################
x0 <- logical(data$p)
temp <- sample.int(data$maxk, 1L)
temp <- sample.int(data$p, temp)
x0[temp] <- TRUE


###################################################
### code chunk number 64: NMOFex.Rnw:1228-1239
###################################################
require(rbenchmark)
benchmark(lm(data$y ~ -1 + data$X[ ,x0]),
          qr.solve(data$X[ ,x0], data$y),
          columns = c("test", "elapsed", "relative"),
          order = "test",
          replications = 1000L)

## ... should give the same coefficients
ignore1 <- lm(data$y ~ -1 + data$X[ ,x0])
ignore2 <- qr.solve(data$X[ ,x0], data$y)
all.equal(as.numeric(coef(ignore1)), as.numeric(ignore2))


###################################################
### code chunk number 65: NMOFex.Rnw:1248-1254
###################################################
OF <- function(x, data) {
    q <- qr(data$X[ ,x])
    e <- qr.resid(q, data$y)
    log(crossprod(e)/data$n) + sum(x) * data$lognn
}
OF(x0, data)


###################################################
### code chunk number 66: NMOFex.Rnw:1259-1268
###################################################
neighbour <- function(xc, data) {
    xn <- xc
    ex <- sample.int(data$p, 1L)
    xn[ex] <- !xn[ex]
    sumx <- sum(xn)
    if ( sumx < 1L || (sumx > data$maxk) )
        xc else xn
}
neighbour(x0, data)[1:10]


###################################################
### code chunk number 67: NMOFex.Rnw:1273-1281
###################################################
algo <- list(
    nT = 10L,    ### number of thresholds
    nS = 200L,   ### number of steps per threshold
    nD = 1000L,  ### number of random steps to compute thresholds
    neighbour = neighbour,
    x0 = x0,
    printBar = FALSE)
system.time(sol1 <- TAopt(OF, algo = algo, data = data))


###################################################
### code chunk number 68: NMOFex.Rnw:1286-1289
###################################################
sol1$OFvalue
which(sol1$xbest)  ### the selected regressors
rD$K               ### the true regressors


###################################################
### code chunk number 69: NMOFex.Rnw:1296-1300
###################################################
xtrue <- logical(data$p)
xtrue[rD$K] <- TRUE
OF(sol1$xbest, data)
OF(xtrue, data)


###################################################
### code chunk number 70: NMOFex.Rnw:1307-1316
###################################################
restarts <- 50L
algo$printDetail <- FALSE
res <- restartOpt(TAopt, n = restarts,
                  OF = OF, algo = algo, data = data,
                  method = "snow", cl = 2)
par(bty = "n", las = 1,mar = c(3,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
plot(ecdf(sapply(res, `[[`, "OFvalue")),  ### extract solution quality
     cex = 0.4, main = "", ylab = "", xlab = "")


###################################################
### code chunk number 71: NMOFex.Rnw:1321-1329
###################################################
xbestAll <- sapply(res, `[[`, "xbest")    ### extract all solutions
inclReg  <- which(rowSums(xbestAll) > 0L) ### get included regressors
inclReg  <- sort(union(rD$K, inclReg))
data.frame(
    regressor = inclReg,
    times_included = paste(rowSums(xbestAll)[inclReg], "/",
                           restarts, sep = ""),
    true_regressor = inclReg %in% rD$K)


###################################################
### code chunk number 72: NMOFex.Rnw:1360-1367
###################################################
size <- 20L
x <- logical(size)
x[runif(size) > 0.5] <- TRUE

## store information
Data <- list()
Data$size <- size


###################################################
### code chunk number 73: NMOFex.Rnw:1371-1379
###################################################
compareLogicals <- function(x, y, ...) {
    argsL <- list(...)
    if (!("sep" %in% names(argsL))) argsL$sep <- ""
    do.call("cat",
            c(list("\n",as.integer(x), "\n", as.integer(y), "\n",
                   ifelse(x == y, " ", "^"), "\n"), argsL)
            )
        }


###################################################
### code chunk number 74: NMOFex.Rnw:1383-1389
###################################################
## there should be no difference
compareLogicals(x, x)

## change the second element
z <- x; z[2L] <- !z[2L]
compareLogicals(x, z)


###################################################
### code chunk number 75: NMOFex.Rnw:1395-1402
###################################################
Data$n <- 5L  ## how many elements to change
neighbour <- function(x, Data) {
    ii <- sample.int(Data$size, Data$n)
    x[ii] <- !x[ii]
    x
}
compareLogicals(x, neighbour(x, Data))


###################################################
### code chunk number 76: NMOFex.Rnw:1411-1422
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
### code chunk number 77: NMOFex.Rnw:1434-1448
###################################################
## match a binary (logical) string y
size <- 20L             ## the length of the string
OF <- function(x, y)    ## the objective function
    sum(x != y)
y <- runif(size) > 0.5  ## the true solution
OF(y, y)                ## the optimum value is zero
algo <- list(nB = size, nP = 50L, nG = 150L, prob = 0.002,
             printBar = FALSE, methodOF = "loop")
system.time(sol <- GAopt(OF, algo = algo, y = y))
OF(sol$xbest, y)
algo <- list(nB = size, nP = 50L, nG = 150L, prob = 0.002,
             printBar = FALSE, methodOF = "snow", cl = 2L)
system.time(sol <- GAopt(OF, algo = algo, y = y))
OF(sol$xbest, y)


###################################################
### code chunk number 78: NMOFex.Rnw:1453-1466
###################################################
OF <- function(x, y) {
    Sys.sleep(0.01)
    sum(x != y)
}
algo <- list(nB = size, nP = 50L, nG = 10L, prob = 0.002,
             printBar = FALSE, methodOF = "loop")
system.time(sol <- GAopt(OF, algo = algo, y = y))
OF(sol$xbest, y)

algo <- list(nB = size, nP = 50L, nG = 10L, prob = 0.002,
             printBar = FALSE, methodOF = "snow", cl = 2L)
system.time(sol <- GAopt(OF, algo = algo, y = y))
OF(sol$xbest, y)


###################################################
### code chunk number 79: NMOFex.Rnw:1506-1513
###################################################
OF <- tfTrefethen
n <- 100L
surf <- matrix(NA, n, n)
x1 <- seq(from = -10, to = 10, length.out = n)
for (i in seq_len(n))
    for (j in seq_len(n))
         surf[i, j] <- tfTrefethen(c(x1[i], x1[j]))


###################################################
### code chunk number 80: NMOFex.Rnw:1520-1526
###################################################
par(bty = "n", las = 1, mar = c(3,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))

## the actual minimum
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))


###################################################
### code chunk number 81: NMOFex.Rnw:1531-1537
###################################################
algo <- list(nP = 50L, nG = 300L,
             F = 0.6, CR = 0.9,
             min = c(-10,-10), max = c(10,10),
             printDetail = FALSE, printBar = FALSE,
             storeF = TRUE, storeSolutions = TRUE)
sol <- DEopt(OF = OF, algo = algo)


###################################################
### code chunk number 82: NMOFex.Rnw:1540-1546
###################################################
names(sol)
sd(sol$popF)
ts.plot(sol$Fmat, xlab = "generations", ylab = "OF")

length(sol$xlist)
xlist <- sol$xlist[[1L]]


###################################################
### code chunk number 83: NMOFex.Rnw:1554-1569
###################################################
## show solution 1 (column 1) in population over time
xlist[[  1L]][ ,1L]  ### at the end of generation 1
## ...
xlist[[ 10L]][ ,1L]  ### at the end of generation 10
## ...
xlist[[300L]][ ,1L]  ### at the end of generation 300

res  <- sapply(xlist, `[`, 1:2, 1)  ### get row 1 and 2 from column 1
res2 <- sapply(xlist, `[`, TRUE, 1) ### simpler
all.equal(res, res2)

dim(res)
res[ ,1L]
res[ ,2L]
res[ ,300L]


###################################################
### code chunk number 84: NMOFex.Rnw:1573-1587
###################################################
## show parameter 2 (row 2) in population over time
xlist[[  1L]][2L, ]  ### at the end of generation 1
## ...
xlist[[ 10L]][2L, ]  ### at the end of generation 10
## ...
xlist[[300L]][2L, ]  ### at the end of generation 300

res <- sapply(xlist, `[`, 2, 1:50)
res <- sapply(xlist, `[`, 2, TRUE)  ### simpler

dim(res)
res[ ,1L]
res[ ,2L]
res[ ,300L]


###################################################
### code chunk number 85: NMOFex.Rnw:1592-1601 (eval = FALSE)
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
### code chunk number 86: NMOFex.Rnw:1604-1632
###################################################
setEPS()
postscript(file = "figures/c1.eps", width = 2, height = 2)
par(bty="n", las = 1,mar = c(2,2,0,0),
     ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))
points(t(xlist[[1L]]), pch = 21, bg=grey(0.9), col = grey(.2))
invisible(dev.off())

setEPS()
postscript(file = "figures/c2.eps", width = 2, height = 2)
par(bty="n", las = 1,mar = c(2,2,0,0),
     ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))
points(t(xlist[[100L]]), pch = 21, bg=grey(0.9), col = grey(.2))
invisible(dev.off())

setEPS()
postscript(file = "figures/c3.eps", width = 2, height = 2)
par(bty="n", las = 1,mar = c(2,2,0,0),
     ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
contour(x1, x1, surf, nlevels=5, col = grey(0.6))
abline(v = -0.02440308, h = 0.21061243, col = grey(0.6))
points(t(xlist[[300L]]), pch = 21, bg=grey(0.9), col = grey(.2))
invisible(dev.off())
cat("\\includegraphics{figures/c1.eps}\\includegraphics{figures/c2.eps}\\includegraphics{figures/c3.eps}",
sep = "")


###################################################
### code chunk number 87: NMOFex.Rnw:1643-1658
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
### code chunk number 88: NMOFex.Rnw:1663-1680
###################################################
algo <- list(nP = 200L, nG = 100L,
             F = 0.50, CR = 0.99,
             min = c( 0,-10,-10,  0),
             max = c( 1, 10, 10, 10),
             storeSolutions = TRUE, printBar = FALSE)


## set up yield curve and put information in Data
tm <- 1:20                ### times to maturity
parTRUE <- c(5, 3, 2, 1)  ### true parameters
yM <- NS(parTRUE, tm)     ### true market yields
Data <- list(yM = yM, tm = tm, model = NS, ww = 0.1, maxb1 = 4)

## solve with DEopt
sol <- DEopt(OF = OF, algo = algo, Data = Data)
P <- sol$xlist[[1L]] ### all population matrices
p1 <- sapply(P, `[`, 1L, TRUE)


###################################################
### code chunk number 89: NMOFex.Rnw:1687-1696
###################################################
par(bty = "n", las = 1, mar = c(4,4,0,0),
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
plot(jitter(rep(seq_len(algo$nG), each = algo$nP), factor = 5),
     p1,
     pch = 21, cex = 0.01, ylim = c(-5,10),
     xlab = "", ylab = "")

mtext("generation", 1, line = 2)
mtext("parameter\nvalue", 2, line = 1)


###################################################
### code chunk number 90: NMOFex.Rnw:1703-1738
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
    ps = 8, tck = 0.001, mgp = c(3, 0.2, 0))
plot(jitter(rep(seq_len(algo$nG), each = algo$nP), factor = 5),
     p1,
     pch = 21, cex = 0.01, ylim = c(-5,10),
     xlab = "", ylab = "")
abline(h = 4, col=grey(0.5))
mtext("generation", 1, line = 2)
mtext("parameter\nvalue", 2, line = 1)


###################################################
### code chunk number 91: NMOFex.Rnw:1750-1760
###################################################
testFun <- function(x) {
    Sys.sleep(0.1) ## wasting time :-)
    cos(1/x^2)
}
system.time(sol1 <- bracketing(testFun, interval = c(0.3, 0.9),
                           n = 100L))
system.time(sol2 <- bracketing(testFun, interval = c(0.3, 0.9),
                           n = 100L, method = "snow", cl = 2))

all.equal(sol1, sol2)


###################################################
### code chunk number 92: NMOFex.Rnw:1791-1792
###################################################
cfBSM


###################################################
### code chunk number 93: NMOFex.Rnw:1796-1815
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
### code chunk number 94: NMOFex.Rnw:1840-1841
###################################################
toLatex(sessionInfo())


