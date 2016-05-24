#Example D: parametric regression quantiles
#illustrating the use of compContourM1u/compContourM2u and evalContour
# on parametric regression quantiles













set.seed(8416)
N <- 1999
M <- 2
Tau <- 0.35
XVec <- matrix(runif(N, -1, 1), N, 1)
YMat <- cbind(XVec, XVec^2) + matrix(runif(N*M, -1, 1), N, M)
X0Vec <- (-4:4)/5




ColPal <- colorRampPalette(c('black', 'red'))
plot(YMat[, 1], YMat[, 2], xlab='y_1', ylab='y_2', type='p', pch='.', cex=3, col=ColPal(N)[as.numeric(cut(XVec, breaks = N))])

COutST = compContourM1u(Tau, YMat, cbind(matrix(1, N, 1), XVec, XVec^2))
if (length(COutST$CharST$HypMat)){
    for (i in 1:length(X0Vec)){
        XRef <- X0Vec[i]
        CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, (M+1):(M+3)]%*%rbind(1, XRef, XRef^2))
        if (CST$Status == 0){
            RGBColor <- rgb(0, 0.5+(XRef/2.5), 0)
            CnvOrd <- sort(atan2(CST$TVVMat[, 1]-mean(CST$TVVMat[, 1]), CST$TVVMat[, 2]-mean(CST$TVVMat[, 2])), index.return=TRUE)
            ContourVec <- CST$TVVMat[CnvOrd$i, ]; ContourVec <- rbind(ContourVec, ContourVec[1, ])
            lines(ContourVec, type='l', col=RGBColor, lwd=2)
        }
    }
}




