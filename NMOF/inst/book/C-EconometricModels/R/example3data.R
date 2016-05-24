# example3data.R -- version 2011-01-03
# ... fit Nelson-Siegel-Svensson to bond prices
# ... continues example1.R

# objective function
OF2 <- function(param, data) {
    tm  <- data$tm; bM  <- data$bM; model <- data$model
    cfMatrix <- data$cfMatrix
    diFa <- 1 / ( (1 + model(param, tm)/100)^tm )
    b <- diFa %*% cfMatrix
    aux <- b - bM; aux <- max(abs(aux))
    if(is.na(aux)) 1e5 else aux
}

cfList <- bundData$cfList
tmList <- bundData$tmList
mats <- unlist(tmList, use.names = FALSE)
mats <- sort(unique(mats))
ISIN <- names(bundData[[1]])

# set up cash flow matrix
nR <- length(mats); nC <- length(cfList)
cfMatrix <- array(0, dim = c(nR, nC))
for(j in seq(nC)) 
    cfMatrix[mats %in% tmList[[j]], j] <- cfList[[j]]
rownames(cfMatrix) <- mats
colnames(cfMatrix) <- ISIN

# reprice bonds with known yield curve
today <- as.Date("2010-05-31")
tm <- as.numeric((as.Date(mats) - today))/365
betaTRUE <- c(5, -2, 1, 10, 1, 5); yM <- NSS(betaTRUE, tm)
diFa <- 1 / ((1 + yM/100)^tm)
bM <- diFa %*% cfMatrix
#bM <- bundData$bM

plot(tm, yM, 
    xlab = "maturities in years",
    ylab = "yields in %")

# collect all in dataList
data <- list(
    bM = bM, 
    tm = tm, 
    cfMatrix = cfMatrix, 
    model = NSS, 
    ww = 1,
    min = c( 0,-15,-30,-30,0,3),
    max = c(15, 30, 30, 30,3,6))

# set parameters for de
algo <- list(
    nP = 100, 
    nG = 600, 
    F = 0.5, 
    CR = 0.9,
    min	= c( 0,-15,-30,-30,0,3),
    max	= c(15, 30, 30, 30,3,6),
    pen = penalty, 
    repair = NULL,
    loopOF = TRUE, 
    loopPen = FALSE, 
    loopRepair = FALSE)

system.time(sol <- DEopt(OF = OF2, algo = algo, data = data))
# maximum yield error
max(abs(data$model(sol$xbest, tm) - data$model(betaTRUE, tm)))
# max. abs. price error and obj. function: should be the same
diFa <- 1 / ((1 + NSS(sol$xbest,tm)/100)^tm)
b <- diFa %*% cfMatrix
max(abs(b - bM))
sol$OFvalue
lines(tm,data$model(sol$xbest, tm), col = "blue")

s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
system.time(sol2 <- nlminb(s0, OF2, data = data,
        lower = data$min, 
        upper = data$max, 
        control = list(eval.max = 50000, iter.max = 50000)))
# maximum yield error
max(abs(data$model(sol2$par, tm) - data$model(betaTRUE, tm)))
# max. abs. price error and obj. function: should be the same
diFa <- 1 / ((1 + NSS(sol2$par, tm)/100)^tm)
b <- diFa %*% cfMatrix
max(abs(b - bM))
sol2$objective	
lines(tm,data$model(sol2$par,tm), col = "green", lty = 2)

legend(x = "bottom", 
    legend = c("true yields", "DE", "nlminb"),
    col = c("black", "blue", "green"),
    pch = c(1, NA, NA), lty = c(0, 1, 2))