#Example E: local(ly) constant (nonparametric) regression quantiles
#illustrating the use of compContourM1u/compContourM2u and evalContour
# on local(ly) constant (nonparametric) regression quantiles













set.seed(8416)
N <- 1999
M <- 2
Tau <- 0.35
XVec <- matrix(runif(N, -1, 1), N, 1)
YMat <- cbind(XVec, XVec^2) + matrix(runif(N*M, -1, 1), N, M)
X0Vec <- (-4:4)/5




ColPals <- colorRampPalette(c('black', 'red'))
plot(YMat[, 1], YMat[, 2], xlab='y_1', ylab='y_2', type='p', pch='.', cex=3, col=ColPals(N)[as.numeric(cut(XVec, breaks = N))])

for (i in 1:length(X0Vec)){
    XRef <- X0Vec[i]
    Sigma <- 0.4
    WVec <- (2*pi*Sigma^2)^(-1/2)*exp(-((XVec-XRef)/Sigma)^2/2)
    WVec <- N*WVec/sum(WVec)
    COutST <- compContourM1u(Tau, matrix(WVec, N, M)*YMat, WVec)
    if (length(COutST$CharST$HypMat)){
        CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, M+1])
        if (CST$Status == 0){
            RGBColor <- rgb(0, 0.5+(XRef/2.5), 0)
            CnvOrd <- sort(atan2(CST$TVVMat[, 1]-mean(CST$TVVMat[, 1]), CST$TVVMat[, 2]-mean(CST$TVVMat[, 2])), index.return=TRUE)
            ContourVec <- CST$TVVMat[CnvOrd$i, ]; ContourVec <- rbind(ContourVec, ContourVec[1, ])
            lines(ContourVec, type='l', col=RGBColor, lwd=2)
        }
    }
}




