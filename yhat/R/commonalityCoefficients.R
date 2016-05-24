
commonalityCoefficients<-
function (dataMatrix, dv, ivlist, imat = FALSE) 
{
    ivlist <- unlist(ivlist)
    ilist<-match(ivlist,colnames(dataMatrix))
    ilist<-na.omit(ilist)
    ilist<-colnames(dataMatrix)[ilist]
    dataMatrix<-na.omit(dataMatrix[,c(dv,ilist)])
    nvar = length(ivlist)
    if (nvar < 2) 
        return("Commonality analysis not conducted. Insufficient number of regressors.")
    ivID <- matrix(nrow = nvar, ncol = 1)
    for (i in 0:nvar - 1) {
        ivID[i + 1] = 2^i
    }
    if (imat) 
        print(ivID)
    numcc = 2^nvar - 1
    effectBitMap <- matrix(0, nvar, numcc)
    for (i in 1:numcc) {
        effectBitMap <- setBits(i, effectBitMap)
    }
    if (imat) 
        print(effectBitMap)
    commonalityMatrix <- matrix(nrow = numcc, ncol = 3)
    for (i in 1:numcc) {
        formula = paste(dv, "~", sep = "")
        for (j in 1:nvar) {
            bit = effectBitMap[j, i]
            if (bit == 1) {
                formula = paste(formula, paste("+", ivlist[[j]], 
                  sep = ""), sep = "")
            }
        }
        commonalityMatrix[i, 2] <- summary(lm(formula, dataMatrix))$r.squared
    }
    if (imat) 
        print(commonalityMatrix)
    commonalityList <- vector("list", numcc)
    for (i in 1:numcc) {
        bit = effectBitMap[1, i]
        if (bit == 1) 
            ilist <- c(0, -ivID[1])
        else ilist <- ivID[1]
        for (j in 2:nvar) {
            bit = effectBitMap[j, i]
            if (bit == 1) {
                alist <- ilist
                blist <- genList(ilist, -ivID[j])
                ilist <- c(alist, blist)
            }
            else ilist <- genList(ilist, ivID[j])
        }
        ilist <- ilist * -1
        commonalityList[[i]] <- ilist
    }
    if (imat) 
        print(commonalityList)
    for (i in 1:numcc) {
        r2list <- unlist(commonalityList[i])
        numlist = length(r2list)
        ccsum = 0
        for (j in 1:numlist) {
            indexs = r2list[[j]]
            indexu = abs(indexs)
            if (indexu != 0) {
                ccvalue = commonalityMatrix[indexu, 2]
                if (indexs < 0) 
                  ccvalue = ccvalue * -1
                ccsum = ccsum + ccvalue
            }
        }
        commonalityMatrix[i, 3] = ccsum
    }
    if (imat) 
        print(commonalityMatrix)
    orderList <- vector("list", numcc)
    index = 0
    for (i in 1:nvar) {
        for (j in 1:numcc) {
            nbits = sum(effectBitMap[, j])
            if (nbits == i) {
                index = index + 1
                commonalityMatrix[index, 1] <- j
            }
        }
    }
    if (imat) 
        print(commonalityMatrix)
    outputCommonalityMatrix <- matrix(nrow = numcc + 1, ncol = 2)
    totalRSquare <- sum(commonalityMatrix[, 3])
    for (i in 1:numcc) {
        outputCommonalityMatrix[i, 1] <- round(commonalityMatrix[commonalityMatrix[i, 
            1], 3], digits = 4)
        outputCommonalityMatrix[i, 2] <- round((commonalityMatrix[commonalityMatrix[i, 
            1], 3]/totalRSquare) * 100, digits = 2)
    }
    outputCommonalityMatrix[numcc + 1, 1] <- round(totalRSquare, 
        digits = 4)
    outputCommonalityMatrix[numcc + 1, 2] <- round(100, digits = 4)
    rowNames = NULL
    for (i in 1:numcc) {
        ii = commonalityMatrix[i, 1]
        nbits = sum(effectBitMap[, ii])
        cbits = 0
        if (nbits == 1) 
            rowName = "Unique to "
        else rowName = "Common to "
        for (j in 1:nvar) {
            if (effectBitMap[j, ii] == 1) {
                if (nbits == 1) 
                  rowName = paste(rowName, ivlist[[j]], sep = "")
                else {
                  cbits = cbits + 1
                  if (cbits == nbits) {
                    rowName = paste(rowName, "and ", sep = "")
                    rowName = paste(rowName, ivlist[[j]], sep = "")
                  }
                  else {
                    rowName = paste(rowName, ivlist[[j]], sep = "")
                    rowName = paste(rowName, ", ", sep = "")
                  }
                }
            }
        }
        rowNames = c(rowNames, rowName)
    }
    rowNames = c(rowNames, "Total")
    rowNames <- format.default(rowNames, justify = "left")
    colNames <- format.default(c("Coefficient", " % Total"), 
        justify = "right")
    dimnames(outputCommonalityMatrix) <- list(rowNames, colNames)
    if (imat) 
        print(outputCommonalityMatrix)
    outputCCbyVar <- matrix(nrow = nvar, ncol = 3)
    for (i in 1:nvar) {
        outputCCbyVar[i, 1] = outputCommonalityMatrix[i, 1]
        outputCCbyVar[i, 3] = round(sum(effectBitMap[i, ] * commonalityMatrix[, 
            3]), digits = 4)
        outputCCbyVar[i, 2] = outputCCbyVar[i, 3] - outputCCbyVar[i, 
            1]
    }
    dimnames(outputCCbyVar) <- list(ivlist, c("Unique", "Common", 
        "Total"))
    outputList <- list(CC = outputCommonalityMatrix, CCTotalbyVar = outputCCbyVar)
    return(outputList)
}
 
