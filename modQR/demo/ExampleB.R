#Example B: 2D location quantiles
#illustrating the use of compContourM1u/compContourM2u and evalContour
# on 2D location quantiles













N <- 1357; M <- 2; TauVec <- c(0.3579, 0.1357, 0.0135)




set.seed(8416); YMat <- matrix(runif(N*M, -1, 1), N, M); XVec <- matrix(1, N, 1)

plot(YMat[, 1], YMat[, 2], xlab='y_1', ylab='y_2', type='p', pch='.', cex=3, col='red')

#Way 1:
for (i in 1:length(TauVec)){
    Tau <- TauVec[i]
    COutST <- compContourM1u(Tau, YMat, XVec)
    if (length(COutST$CharST$HypMat)){
        CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, (M+1)])
        if (CST$Status == 0){
            CnvOrd <- sort(atan2(CST$TVVMat[, 1]-mean(CST$TVVMat[, 1]), CST$TVVMat[, 2]-mean(CST$TVVMat[, 2])), index.return=TRUE)
            ContourVec <- CST$TVVMat[CnvOrd$i, ]; ContourVec <- rbind(ContourVec, ContourVec[1, ])
            lines(ContourVec, type='l', col='#008000', lwd=8) #color green, line width 8
        }
    }
}

set.seed(8416); YMat <- matrix(runif(N*M, -1, 1), N, M); XVec <- matrix(1, N, 1)

#Way 2:
for (i in 1:length(TauVec)){
    Tau <- TauVec[i]
    COutST <- compContourM1u(Tau, YMat, XVec)
    if (length(COutST$CharST$HypMat)){
        CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, (M+1)])
        if (CST$Status == 0){
            CnvOrd <- sort(atan2(CST$TVVMat[, 1]-mean(CST$TVVMat[, 1]), CST$TVVMat[, 2]-mean(CST$TVVMat[, 2])), index.return=TRUE)
            ContourVec <- CST$TVVMat[CnvOrd$i, ]; ContourVec <- rbind(ContourVec, ContourVec[1, ])
            lines(ContourVec, type='l', col='#0000a0', lwd=4) #color blue, line width 4
        }
        IsInVec <- (COutST$PosVec == 0)
        NIn <- N - sum(COutST$PosVec > 0)
        if (NIn > 0){
            PseudoIPRow <- colSums(YMat[IsInVec, ])
            YMat <- YMat[!IsInVec, ]
            XVec <- as.matrix(XVec[!IsInVec, ])
            XVec <- rbind(XVec, NIn)
            YMat <- rbind(YMat, PseudoIPRow)
        }
    }
}

set.seed(8416); YMat <- matrix(runif(N*M, -1, 1), N, M); XVec <- matrix(1, N, 1)

#Way 3:
for (i in 1:length(TauVec)){
    Tau <- TauVec[i]
    Tau <- Tau*N/dim(YMat)[1]
    COutST <- compContourM1u(Tau, YMat, XVec)
    if (length(COutST$CharST$HypMat)){
        CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, (M+1)])
        if (CST$Status == 0){
            CnvOrd <- sort(atan2(CST$TVVMat[, 1]-mean(CST$TVVMat[, 1]), CST$TVVMat[, 2]-mean(CST$TVVMat[, 2])), index.return=TRUE)
            ContourVec <- CST$TVVMat[CnvOrd$i, ]; ContourVec <- rbind(ContourVec, ContourVec[1, ])
            lines(ContourVec, type='l', col='#ffffff', lwd=2) #color white, line width 2
        }
        IsInVec <- (COutST$PosVec == 0)
        YMat <- YMat[!IsInVec, ]
        XVec <- as.matrix(XVec[!IsInVec, ])
    }
}




