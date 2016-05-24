require(GeneralizedHyperbolic)
### Unit tests of function qgigIBF

### Functions with name test.* run by R CMD check or by make if
### LEVEL = 1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL = n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL = graphics in call to make


test.qgigIBF <- function()
{
    ## select a random parameter value for testing

    paramSampleSize <- 1
    sampleSize <- 1
    ps <- c(-Inf, -1, 0, 1, Inf, NA, NA)
    qs <- c(0, 0.001, 0.025, 0.3, 0.5, 0.7, 0.975, 0.999, 1)
    nps <- length(ps)
    nqs <- length(qs)

    data(gigParam)
    testParam <- gigLargeParam

    ## sample parameter values
    np <- nrow(testParam)[1]
    paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
    paramVals <- testParam[paramNum,,drop = FALSE]
    nv <- nrow(paramVals)[1]

    ## Test using default values
    ## Should all be < 10^(-4)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 3 + nqs)
    pqResult[,1:3] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 3 + nps)
    qpResult[,1:3] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv)
    {
        param <- as.numeric(paramVals[i,])
        result <- pgigIBF(qgigIBF(qs, param = param, method = "integrate"),
                          param = param) - qs
        pqResult[i, 3 + (1:nqs)] <- result

        ## choose sample values from distribution
        x <- rgig(sampleSize, param = param)

        ## add some extreme values
        extreme <- qgigIBF(c(10^(-5), 1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x) * c(-1, 1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qgigIBF(pgigIBF(ps, param = param), param = param)

        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 3 + (1:nps)] <- result
    }
    ## first test: p then q
    pqMaxDiff <- max(abs(pqResult[, 3 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pgResult[, 3 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-4), msg =
              paste("pgMaxDiff = ", pgMaxDiff, "for param = ",
                    pqResult[pgMaxInd, 1], pqResult[pgMaxInd, 2],
                    pqResult[pgMaxInd, 3], pqResult[pgMaxInd, 4],
                    pqResult[pgMaxInd, 5]))

    ## second test: 1 then p
    qpMaxDiff <- max(abs(qpResult[, 3 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[, 3 + (1:nps)])))
    }
    checkTrue(qpMaxDiff < 10^(-4), msg =
              paste("qpMaxDiff =", qpMaxDiff, "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4],
                    qpResult[qpMaxInd, 5]))
    ## Now test with higher accuracy
    ## With spline will still only be 10^(-4)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 3 + nqs)
    pqResult[, 1:3] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 3 + nps)
    qpResult[, 1:3] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv)
    {
        param <- as.numeric(paramVals[i,])
        result <- pgigIBF(qgigIBF(qs, param = param, uniTol = 10^(-12)),
                        param = param, ibfTol = 10^(-12)) - qs
        pqResult[i, 3 + (1:nqs)] <- result

        ## choose sample values from distribution
        x <- rgig(sampleSize, param = param)

        ## add some extreme values
        extreme <- qgigIBF(c(10^(-5), 1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x) * c(-1, 1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qgigIBF(pgigIBF(ps, param = param, ibfTol = 10^(-12)),
                        param = param, uniTol = 10^(-12))
        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 3 + (1:nps)] <- result
    }

    ## third test: p then q
    pqMaxDiff <- max(abs(pqResult[, 3 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[, 3 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-6), msg =
              paste("pgMaxDiff = ", pgMaxDiff, "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4],
                    pqResult[pqMaxInd, 5]))

    ## fourth test: 1 then p
    qpMaxDiff <- max(abs(pqResult[, 3 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[, 3 + (1:nps)])))
    }
    checkTrue(qpMaxDiff < 10^(-4), msg =
              paste("qpMaxDiff = ", qpMaxDiff, "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4],
                    qpResult[qpMaxInd, 5]))

    ## Now try with integration instead of spline
    ## Reset parameter set
    data(gigParam)
    testParam <- gigLargeShape

    ## sample parameter values
    ## sample parameter values
    np <- nrow(testParam)[1]
    paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
    paramVals <- testParam[paramNum,,drop = FALSE]
    nv <- nrow(paramVals)[1]


    ## Test using default values
    ## Should all be < 10^(-4)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 3 + nqs)
    pqResult[, 1:3] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 3 + nps)
    qpResult[, 1:3] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv)
    {
        param <- as.numeric(paramVals[i,])
        result <- pgigIBF(qgigIBF(qs, param = param, method = "integrate"),
                        param = param) - qs
        pqResult[i, 3 + (1:nqs)] <- result

        ## choose sample values from distribution
        x <- rgig(sampleSize, param = param)

        ## add some extreme values
        extreme <- qgigIBF(c(10^(-5), 1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x) * c(-1, 1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qgigIBF(pgigIBF(ps, param = param),
                        param = param, method = "integrate")

        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 3 + (1:nps)] <- result
    }

    ## fifth test: p then q
     pqMaxDiff <- max(abs(pqResult[, 3 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[, 3 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-4), msg =
              paste("pqMaxDiff = ", pqMaxDiff,
                    "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4],
                    pqResult[pqMaxInd, 5]))

    ## sixth test: 1 then p
    qpMaxDiff <- max(abs(qpResult[, 3 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[, 3 + (1:nps)])))
    }
    checkTrue(qpMaxDiff < 10^(-4), msg =
              paste("qpMaxDiff = ", qpMaxDiff,
                    "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4],
                    qpResult[qpMaxInd, 5]))


    ## Now test with higher accuracy
    ## With integrate should be to 10^(-10)

    ## initialize result matrices
    pqResult <- matrix(nrow = nv, ncol = 3 + nqs)
    pqResult[, 1:3] <- as.matrix(paramVals)
    qpResult <- matrix(nrow = nv, ncol = 3 + nps)
    qpResult[, 1:3] <- as.matrix(paramVals)

    ## loop over param values
    for (i in 1:nv){
        param <- as.numeric(paramVals[i,])
        result <- pgigIBF(qskewhyp(qs, param = param, method = "integrate",
                                    uniTol = 10^(-12)),
                           param = param, ibfTol = 10^(-12)) - qs
        pqResult[i, 3 + (1:nqs)] <- result


        ## choose sample values from distribution
        x <- rgig(sampleSize, param = param)

        ## add some extreme values
        extreme <- qgigIBF(c(10^(-5),1 - 10^(-5)), param = param)
        extreme <- extreme + abs(x)*c(-1,1)
        ps[c(nps - 1, nps)] <- extreme
        result <- qgigIBF(pgigIBF(ps, param = param, ibfTol = 10^(-10)),
                           param = param, method = "integrate",
                           uniTol = 10^(-10))

        ## special treatment for +/- Inf
        result <- ifelse(result == ps, 0, result - ps)
        qpResult[i, 3 + (1:nps)] <- result
    }
    ## seventh test: p then q
    pqMaxDiff <- max(abs(pqResult[, 3 + (1:nqs)]))
    if (nv == 1){
        pqMaxInd <- 1
    } else {
        pqMaxInd <- which.max(pmax(abs(pqResult[, 3 + (1:nqs)])))
    }
    checkTrue(pqMaxDiff < 10^(-6), msg =
              paste("pqMaxDiff = ", pqMaxDiff,
                    "for param = ",
                    pqResult[pqMaxInd, 1], pqResult[pqMaxInd, 2],
                    pqResult[pqMaxInd, 3], pqResult[pqMaxInd, 4],
                    pqResult[pqMaxInd, 5]))

    ## eighth test: 1 then p
    qpMaxDiff <- max(abs(qpResult[, 3 + (1:nps)]))
    if (nv == 1){
        qpMaxInd <- 1
    } else {
        qpMaxInd <- which.max(pmax(abs(qpResult[, 3 + (1:nps)])))
    }
    checkTrue(qpMaxDiff < 10^(-4), msg =
              paste("qpMaxDiff = ", qpMaxDiff,
                    "for param = ",
                    qpResult[qpMaxInd, 1], qpResult[qpMaxInd, 2],
                    qpResult[qpMaxInd, 3], qpResult[qpMaxInd, 4],
                    qpResult[qpMaxInd, 5]))

    return()
}



