GWsignif <-
function(pvalues, ngroup = 5, ntest.genome, alpha = 0.05, plot.figure = TRUE)
{
if (is.character(pvalues)) {
pvalues.files <- pvalues
pvaluesWasImported <- FALSE
ntest.region <- 0
for (f in pvalues.files){
pvalues <- as.matrix(read.table(f)) #no header in files! header = FALSE!
ntest.region <- ntest.region + ncol(pvalues) - sum(apply(is.na(pvalues), 2, any))  
}
} else {
pvalues <- as.matrix(pvalues)
ntest.region <- sum(apply(!is.na(pvalues), 2, all))  
pvalues.files <- "pvaluesWasImported"
pvaluesWasImported <- TRUE
}

#ntest.genome
if (missing(ntest.genome)) ntest.genome <- ntest.region
ntest.genome <- max(ntest.genome, ntest.region)

#modify ngroup
ngroup <- min(ngroup, floor(1 + log(ntest.region)/log(4))) 
m1 <- ceiling(ntest.region/(2^(ngroup-1))) 
stopifnot(ngroup > 1)

#minimum pvalues in group 1 
minp <- matrix(NA, nrow(pvalues), 2^(ngroup-1))
pval0 <- NULL
testj <- 0
for (f in pvalues.files){
if (!pvaluesWasImported) pvalues <- as.matrix(read.table(f)) #no header in files! header = FALSE!
pvalues <- pvalues[, !apply(is.na(pvalues), 2, any)] 
pvalues <- cbind(pval0, pvalues)
nregion <- floor(ncol(pvalues)/m1)

if (nregion > 0){
for (j in 1:nregion) {
testj <- testj + 1
minp[,testj] <- apply(pvalues[, (1 + (j-1)*m1):(j*m1)], 1, min)
}
}

if (ncol(pvalues) - nregion*m1 > 0){
pval0 <- pvalues[, (1 + nregion*m1):ncol(pvalues), drop=FALSE]
} else {pval0 <- NULL}
}
if (ncol(pvalues) - nregion*m1 > 0){
testj <- testj + 1
minp[,testj] <- apply(pvalues[, max(1, ncol(pvalues) - m1 +1):ncol(pvalues), drop=FALSE], 1, min)
}

stopifnot(testj == 2^(ngroup-1))

#quantile of minimum pvalues: significance threshold
qminp <- list()
qminp[[1]] <- apply(minp, 2, quantile, prob = alpha)
k <- 2
while (k < ngroup){
minp <- pmin(minp[, seq(1, ncol(minp), by=2)], minp[, seq(2, ncol(minp), by=2)])
qminp[[k]] <- apply(minp, 2, quantile, prob = alpha)
k <- k + 1
}
minp <- pmin(minp[,1], minp[,2])
qminp[[k]] <- quantile(minp, prob = alpha)

#number of tests, mean and standard error of -log10(qminp)
#bonf: –log10 of Bonferroni correction
ntest <- c(m1*2^(0:(ngroup-2)), ntest.region)
x <-  - log10(alpha/ntest) 
logqminp <- lapply(qminp, log10)
y <-  - unlist(lapply(logqminp, mean)) 
ysd <- unlist(lapply(logqminp, sd)) 

mlogq <- cbind(ntest =  ntest,  bonf = x,mean = y, sd = ysd)

#linear model
fit<- lm(y ~ x)
xgw <- - log10(alpha/ntest.genome)
ygw <- predict(fit, new=data.frame(x = xgw))
signifgw <- as.numeric(10^(-ygw))

cat(paste0("\nThe number of tests in a large genome-wide region of interest: ntest.genome=", ntest.genome, 
".\nThe significance threshold in the large genome-wide region: GWsignif.threshold=", signif(signifgw), ".\n"))

#plot model fitting
if (plot.figure) 
{
qk <- unlist(qminp)
nk <- NULL
for (k in 1:ngroup) nk <- c(nk, rep(ntest[k], 2^(ngroup-k)))
xk <- - log10(alpha/nk)
yk <- -log10(qk)

plot(xk, yk, xlim=range(c(x, xk, xgw)), ylim=range(c(y, y+ysd, y-ysd, yk, ygw), na.rm = TRUE), col = "gray", cex=0.8,
xlab=paste0("-log10(FWER / number of tests), FWER=", alpha), ylab="-log10(significance threshold)")
abline(0,1, col="gray")
points(x,y, pch = 16, col="blue")
points(xgw, ygw, pch=8, col="red", cex=1.2)
ypre <- predict(fit)
lines(c(x, xgw), c(ypre, ygw), lty="dashed")
co<- round(coefficients(fit), digits=4)
lege=paste("y = ", co[1], " + ", co[2], "x", sep="")
legend("topleft", lege, lty = "dashed", lwd=1, cex = 1, text.col = "black",  merge = TRUE, bg = 'gray90')

cat("\n--------------
In the figure:
Gray circles represent sub-regions significance thresholds.
Blue dots are the mean of -log10(significance thresholds) in the sub-regions of same size.
Red star represents the extrapolated significance threshold in a large genome-wide region of interest.
The dashed line is fit to the blue solid dots. The grey line is the line of equality, y=x.     
\n")
}

list(qminp = qminp, mlogq = mlogq, alpha = alpha, ntest.region = ntest.region, ntest.genome = ntest.genome, GWsignif.threshold = signifgw)
}
