### Unit tests of function qskewhyp

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

level2test.qskewhyp <- function()
{
    ## Purpose: Level 2 test of qskewhyp and pskewhyp
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: David Scott, Date: 29 Mar 2010, 15:27

    ## select a random parameter value for testing
    ## eliminate problem values first
    paramSampleSize <- 1
    sampleSize <- 1
    ps <- c(-Inf,-1,0,1,Inf,NA,NA) # 2 values to come later
    qs <- c(0,0.001,0.025,0.3,0.5,0.7,0.975,0.999,1)
    nps <- length(ps)
    nqs <- length(qs)

    data(skewhypParam)
    testParam <- skewhypLargeParam
    ## prune difficult values
    difficult <- matrix(c(-5,1,
                          -5,2,
                          -5,5,
                          0,1,
                          0,2,
                          1,1,
                          1,2,
                          2,1,
                          2,2,
                          2,5,
                          5,1,
                          5,2,
                          5,5,
                          5,10,
                          10,1,
                          10,2,
                          10,5,
                          20,1,
                          20,2,
                          20,5,
                          20,10),
                        ncol = 2, byrow = TRUE)
    vTest <- function(x, lookFor) {
        which(apply(x, 1, function(v) all(v == lookFor)))
    }
    delete <- integer()
    nd <- NROW(difficult)
    for (i in 1:nd){
         indices <- vTest(testParam[,3:4], difficult[i,])
         delete <- c(delete, indices)
     }
    delete <- unique(delete) # probably not necessary
    testParam <- testParam[-delete,]

    ## sample parameter values
    np <- NROW(testParam)
    paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
    paramVals <- testParam[paramNum,]
    ##nv <- dim(paramVals)[1]
    nv <- NROW(paramVals)

    ## Test using default values
    ## Should all be < 10^(-4)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 4 + nqs)
    pqResult[,1:4] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 4 + nps)
    qpResult[,1:4] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv){
        param <- as.numeric(paramVals[i,])
        result <- pskewhyp(qskewhyp(qs, param = param),
                           param = param) - qs
        pqResult[i, 4 + (1:nqs)] <- result


        ## choose sample values from distribution
        x <- rskewhyp(sampleSize, param = param)

        ## add some extreme values
        extreme <- qskewhyp(c(10^(-5),1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x)*c(-1,1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qskewhyp(pskewhyp(ps, param = param),
                           param = param)
        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 4 + (1:nps)] <- result
    }
    ## first test: p then q
    pqMaxDiff <- max(abs(pqResult[,4 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[,4 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-4), msg =
              paste("pqMaxDiff = ", pqMaxDiff,
                    "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4]))

    ## second test: 1 then p
    qpMaxDiff <- max(abs(qpResult[,4 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[,4 + (1:nqs)])))
    }
    checkTrue(qpMaxDiff < 10^(-3), msg =
              paste("qpMaxDiff = ", qpMaxDiff,
                    "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4]))


    ## Now test with higher accuracy
    ## With spline will still only be to 10^(-4)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 4 + nqs)
    pqResult[,1:4] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 4 + nps)
    qpResult[,1:4] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv){
        param <- as.numeric(paramVals[i,])
        result <- pskewhyp(qskewhyp(qs, param = param, uniTol = 10^(-12)),
                           param = param, intTol = 10^(-12)) - qs
        pqResult[i, 4 + (1:nqs)] <- result


        ## choose sample values from distribution
        x <- rskewhyp(sampleSize, param = param)

        ## add some extreme values
        extreme <- qskewhyp(c(10^(-5),1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x)*c(-1,1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qskewhyp(pskewhyp(ps, param = param, intTol = 10^(-12)),
                           param = param, uniTol = 10^(-12))
        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 4 + (1:nps)] <- result
    }
    ## third test: p then q
    pqMaxDiff <- max(abs(pqResult[,4 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[,4 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-6), msg =
              paste("pqMaxDiff = ", pqMaxDiff,
                    "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4]))

    ## fourth test: 1 then p
    qpMaxDiff <- max(abs(qpResult[,4 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[,4 + (1:nqs)])))
    }
    checkTrue(qpMaxDiff < 10^(-3), msg =
              paste("qpMaxDiff = ", qpMaxDiff,
                    "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4]))


    ## Now try with integration instead of spline
    ## Reset parameter set
    data(skewhypParam)
    testParam <- skewhypLargeParam
    ## prune difficult values
    difficult <- matrix(c(-5,1,
                          -5,2,
                          0,1,
                          0,2,
                          1,1,
                          1,2,
                          2,1,
                          2,2,
                          2,5,
                          5,1,
                          5,2,
                          5,5,
                          5,10,
                          10,1,
                          10,2,
                          10,5,
                          20,1,
                          20,2,
                          20,5,
                          20,10),
                        ncol = 2, byrow = TRUE)

    vTest <- function(x, lookFor) {
        which(apply(x, 1, function(v) all(v == lookFor)))
    }
    delete <- integer()
    nd <- NROW(difficult)
    for (i in 1:nd){
         indices <- vTest(testParam[,3:4], difficult[i,])
         delete <- c(delete, indices)
     }
    delete <- unique(delete) # probably not necessary
    testParam <- testParam[-delete,]

    ## sample parameter values
    np <- NROW(testParam)
    paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
    paramVals <- testParam[paramNum,]
    nv <- NROW(paramVals)

    ## Test using default values with integrate
    ## Should all be < 10^(-4)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 4 + nqs)
    pqResult[,1:4] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 4 + nps)
    qpResult[,1:4] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv){
        param <- as.numeric(paramVals[i,])
        result <- pskewhyp(qskewhyp(qs, param = param, method = "integrate"),
                           param = param) - qs
        pqResult[i, 4 + (1:nqs)] <- result


        ## choose sample values from distribution
        x <- rskewhyp(sampleSize, param = param)

        ## add some extreme values
        extreme <- qskewhyp(c(10^(-5),1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x)*c(-1,1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qskewhyp(pskewhyp(ps, param = param),
                           param = param, method = "integrate")
        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 4 + (1:nps)] <- result
    }
    ## fifth test: p then q
    pqMaxDiff <- max(abs(pqResult[,4 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[,4 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-4), msg =
              paste("pqMaxDiff = ", pqMaxDiff,
                    "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4]))

    ## sixth test: 1 then p
    qpMaxDiff <- max(abs(qpResult[,4 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[,4 + (1:nqs)])))
    }
    checkTrue(qpMaxDiff < 10^(-3), msg =
              paste("qpMaxDiff = ", qpMaxDiff,
                    "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4]))


    ## Now test with higher accuracy
    ## With integrate should be to 10^(-10)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 4 + nqs)
    pqResult[,1:4] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 4 + nps)
    qpResult[,1:4] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv){
        param <- as.numeric(paramVals[i,])
        result <- pskewhyp(qskewhyp(qs, param = param, method = "integrate",
                                    uniTol = 10^(-12)),
                           param = param, intTol = 10^(-12)) - qs
        pqResult[i, 4 + (1:nqs)] <- result


        ## choose sample values from distribution
        x <- rskewhyp(sampleSize, param = param)

        ## add some extreme values
        extreme <- qskewhyp(c(10^(-5),1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x)*c(-1,1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qskewhyp(pskewhyp(ps, param = param, intTol = 10^(-10)),
                           param = param, method = "integrate",
                           uniTol = 10^(-10))
        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 4 + (1:nps)] <- result
    }
    ## seventh test: p then q
    pqMaxDiff <- max(abs(pqResult[,4 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[,4 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-6), msg =
              paste("pqMaxDiff = ", pqMaxDiff,
                    "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4]))

    ## eighth test: 1 then p
    qpMaxDiff <- max(abs(qpResult[,4 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[,4 + (1:nqs)])))
    }
    checkTrue(qpMaxDiff < 10^(-3), msg =
              paste("qpMaxDiff = ", qpMaxDiff,
                    "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4]))

    return()
}
