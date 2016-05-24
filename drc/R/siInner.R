siInner <- function(indPair, pVec, compMatch, object, indexMat, parmMat, varMat, level, reference, type, sifct, interval, degfree, logBase)
{
    jInd <- indPair[1]
    kInd <- indPair[2]

#    indexMat <- object$"indexMat"
    parmInd1 <- indexMat[, jInd]
    parmInd2 <- indexMat[, kInd]
#    parmMat <- matrix(coef(object)[indexMat], ncol = ncol(indexMat))

    parmChosen1 <- parmMat[, jInd]
    parmChosen2 <- parmMat[, kInd]

#    SIeval <- sifct(parmChosen1, parmChosen2, pVec, 1,2,1,2, reference, type, jInd, kInd)
    SIeval <- sifct(parmChosen1, parmChosen2, pVec, jInd, kInd, reference, type)

    SIval <- SIeval$"val"  # SIeval[[1]]
    dSIval <- SIeval$"der"  # SIeval[[2]]
#    print(dSIval)

#    print(varMat)
    oriMatRow <- c(SIval, sqrt(t(dSIval) %*% varMat %*% dSIval))
    siMatRow <- matrix(NA, 1, 4)  # four is the maximum number of columns
    siMatRow[1, 1] <- SIval
    
    ## Using t-distribution for continuous data
    ##  only under the normality assumption
    if (identical(object$"type", "continuous"))
    {
        qFct <- function(x) {qt(x, degfree)}
        pFct <- function(x) {pt(x, degfree)}
    } else {
        qFct <- qnorm
        pFct <- pnorm
    }

    if (identical(interval, "none"))
    {
        siMatRow[2] <- oriMatRow[2]  # sqrt(dSIval%*%varCov%*%dSIval)

        ## Testing SI equal to 1
        tempStat <- (siMatRow[1] - 1)/siMatRow[2]
        siMatRow[3] <- tempStat
        siMatRow[4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))
    }
    if ( (identical(interval, "delta")) || (identical(interval, "fls")) )
    {
        stErr <- oriMatRow[2]  # sqrt(derEval%*%varCov%*%derEval)
        tquan <- qFct(1 - (1 - level)/2)

        siMatRow[2] <- siMatRow[1] - tquan * stErr
        siMatRow[3] <- siMatRow[1] + tquan * stErr
        ciLabel <- "Delta method"
    }
    if (identical(interval, "tfls"))
    {
        lsVal <- log(oriMatRow[1])
        lsdVal <- oriMatRow[2] / oriMatRow[1]
        tquan <- qFct(1 - (1 - level)/2)

        siMatRow[2] <- exp(lsVal - tquan * lsdVal)
        siMatRow[3] <- exp(lsVal + tquan * lsdVal)
        ciLabel <- "To and from log scale"
    }
    if ((!is.null(logBase)) && (identical(interval, "fls")))
    {
        siMatRow[1] <- logBase^(siMatRow[1])
        siMatRow[2] <- logBase^(siMatRow[2])
        siMatRow[3] <- logBase^(siMatRow[3])
        ciLabel <- "From log scale"
    }
    if (identical(interval, "fieller"))  # using t-distribution
    {
        vcMat <- matrix(NA, 2, 2)
#        vcMat[1, 1] <- SIeval$"der1" %*% varMat[parmInd1, parmInd1] %*% SIeval$"der1"
#        vcMat[2, 2] <- SIeval$"der2" %*% varMat[parmInd2, parmInd2] %*% SIeval$"der2"
#        vcMat[1, 2] <- SIeval$"der1" %*% varMat[parmInd1, parmInd2] %*% SIeval$"der2"        
        vcMat[1, 1] <- SIeval$"der1" %*% varMat %*% SIeval$"der1"
        vcMat[2, 2] <- SIeval$"der2" %*% varMat %*% SIeval$"der2"
        vcMat[1, 2] <- SIeval$"der1" %*% varMat %*% SIeval$"der2"
        vcMat[2, 1] <- vcMat[1, 2]
        muVec <- c(SIeval$"valnum", SIeval$"valden")

        siMatRow[2:3] <- fieller(muVec, degfree, vcMat, level = level)
        ciLabel <- "Fieller"
    }
    c(siMatRow, dSIval)
}