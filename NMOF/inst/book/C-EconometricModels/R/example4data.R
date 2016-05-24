# example4data.R -- version 2011-01-03
# ... fit Nelson-Siegel-Svensson to bond yields
# ... continues example1.R

compYield <- function(cf, tm, guess = NULL) {
    fy <- function(ytm, cf, tm) sum(cf / ((1 + ytm)^tm))
    logik <- cf != 0
    cf <- cf[logik]
    tm <- tm[logik]
    if (is.null(guess))
        ytm <- 0.05 else ytm <- guess
    h <- 1e-8; dF <- 1; ci <- 0L
    while (abs(dF) > 1e-5) {
        ci <- ci + 1L
        if (ci > 5L) break
        FF <- fy(ytm, cf, tm)
        dFF <- (fy(ytm + h, cf, tm) - FF) / h
        dF <- FF / dFF
        ytm <- ytm - dF
    }
    ytm
}
OF3 <- function(param, data) {       
    tm <- data$tm; rM <- data$rM; model <- data$model
    cfMatrix<- data$cfMatrix; nB <- dim(cfMatrix)[2L]
    zrates <- model(param, tm)
    aux <- 1e8
    diFa    <- 1/((1 + zrates/100)^tm)
    b <- diFa %*% cfMatrix
    r <- numeric(nB)
    if (all(!is.na(b), diFa < 1, diFa > 0, b > 1)) { 
        # compute theoretical yields
        for (bb in 1L:nB) {
            if (bb == 1L) 
                guess <- 0.05 else guess <- r[bb - 1L]
            r[bb] <- compYield(c(-b[bb], cfMatrix[ ,bb]), 
                c(0, tm), guess)  
        }
        aux <- abs(r - rM)
        aux <- sum(aux)
    }
    aux
}

cfList <- bundData$cfList
tmList <- bundData$tmList
mats   <- unlist(tmList, use.names = FALSE)
mats   <- sort(unique(mats))
ISIN   <- names(bundData[[1]])

# set up cash flow matrix
nR <- length(mats); nC <- length(cfList)
cfMatrix <- array(0, dim = c(nR, nC))
for(j in seq(nC)) 
    cfMatrix[mats %in% tmList[[j]], j] <- cfList[[j]]
rownames(cfMatrix) <- mats
colnames(cfMatrix) <- ISIN

# compute artificial market prices
today <- as.Date("2010-05-31")
tm <- as.numeric((as.Date(mats) - today))/365
betaTRUE <- c(5,-2,1,10,1,3); yM <- NSS(betaTRUE, tm)
diFa <- 1 / ((1 + yM/100)^tm)
bM <- diFa %*% cfMatrix
rM <- apply(rbind(-bM, cfMatrix), 2, compYield, c(0,tm))
plot(tm, yM, 
    xlab = "maturities in years", 
    ylab = "yields in %")

# collect all in dataList
data <- list(
    rM = rM, 
    tm = tm, 
    cfMatrix = cfMatrix,
    model = NSS,
    min = c( 0,-15,-30,-30,0  ,2.5),
    max = c(15, 30, 30, 30,2.5,  5), 
    ww = 0.1
)

# set parameters for de
algo <- list(
    nP = 50L, 
    nG = 500L,
    F = 0.50, 
    CR = 0.99,
    min = c( 0,-15,-30,-30,0  ,2.5),
    max = c(15, 30, 30, 30,2.5,5  ),
    pen = penalty, repair = NULL,
    loopOF = TRUE, loopPen = FALSE, 
    loopRepair = FALSE, printDetail = TRUE)

system.time(sol <- DEopt(OF = OF3, algo = algo, data = data))

# maximum error
max(abs(data$model(sol$xbest,tm) - data$model(betaTRUE,tm)))
# maximum abs. yield error and objective function
diFa <- 1 / ((1 + NSS(sol$xbest,tm)/100)^tm)
b <- diFa %*% cfMatrix
r <- apply(rbind(-b, cfMatrix), 2, compYield,c(0, tm))
sum(abs(r - rM))
sol$OFvalue
lines(tm,data$model(sol$xbest, tm), col = "blue")

s0 <- algo$min + (algo$max - algo$min)*runif(length(algo$min))
system.time(sol2 <- nlminb(s0, OF3, data = data,
        lower = algo$min, 
        upper = algo$max, 
        control = list(eval.max = 50000L, iter.max = 50000L)))
# maximum error
max(abs(data$model(sol2$par,tm) - data$model(betaTRUE,tm)))
# maximum abs. yield error and objective function
diFa <- 1 / ((1 + NSS(sol2$par,tm)/100)^tm)
b <- diFa %*% cfMatrix
r <- apply(rbind(-b, cfMatrix), 2, compYield,c(0, tm))
sum(abs(r - rM))
sol2$objective
lines(tm, data$model(sol2$par,tm), col = "green", lty = 2)

legend(x = "bottom", 
    legend = c("true yields", "DE", "nlminb"),
    col = c("black", "blue", "green"),
    pch = c(1, NA, NA), lty = c(0, 1, 2))
