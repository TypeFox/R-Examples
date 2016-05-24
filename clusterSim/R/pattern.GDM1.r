pattern.GDM1<-function (data, performanceVariable, scaleType = "i", nomOptValues = NULL, 
    weightsType = "equal", weights = NULL, normalization = "n0", 
    patternType = "upper", patternCoordinates = "dataBounds", 
    patternManual = NULL, nominalTransfMethod = NULL) 
{
    data <- as.matrix(data)
    pdata <- data
    normalizations <-  c("n0", "n1", "n2", "n3", "n3a", "n4", "n5", "n5a", "n6",
           "n6a", "n7", "n8", "n9", "n9a", "n10", "n11", "n12", "n12a", "n13")
    for (v in performanceVariable) {
        if (sum(v == c("s", "n", "d")) == 0) 
            stop("performance variable should be one of the following: s-stimuli,d-destimuli or n-nominant variable")
    }
    if (length(performanceVariable) != ncol(data)) 
        stop("performance variable vector should have the size equal to numberof variables")
    if (length(scaleType) == 1) 
        scaleType <- rep(scaleType, ncol(data))
    if (length(scaleType) != ncol(data)) 
        stop("scale types vector should have the size equal to number of variables")
    for (v in scaleType) {
        if (sum(v == c("r", "i")) == 0) 
            stop("scale type should be one of the following: r-ratio i - interval")
    }
    if (sum(performanceVariable == "n") == 0) 
        nomOptValues <- rep(0, ncol(data))
    if (length(nomOptValues) != ncol(data)) 
        stop("vector of optimal values for nominant variables should have the size equal to numberof variables")
    if (sum(patternType == c("upper", "lower")) == 0) 
        stop("pattern should be one of the following:  \"upper\",\"lower\"")
    if (sum(patternCoordinates == c("dataBounds", "bounds", "manual")) == 
        0) 
        stop("pattern should be one of the following:  \"dataBounds\",\"bounds\",\"manual\"")
    if (patternCoordinates == "manual" && length(patternManual) != 
        ncol(data)) 
        stop("pattern value vector should have the size equal to number of variables")
    if (sum(performanceVariable == "n") == 0) 
        nominalTransfMethod <- rep("q", ncol(data))
    if (length(nominalTransfMethod) == 1) 
        nominalTransfMethod = rep(nominalTransfMethod, ncol(data))
    if (length(nominalTransfMethod) != ncol(data)) 
        stop("vector of transfer methods for nominant variables should have the size equal to numberof variables")
    for (v in nominalTransfMethod) {
        if (sum(v == c("q", "d")) == 0) 
            stop("transfer methods should be one of the following: q-quotient,d-difference")
    }
    if (sum(normalization == normalizations) != 1) 
        stop(paste("Normalization method should be one of the following:", 
            paste(normalizations, collapse = ",")))
    for (i in 1:ncol(data)) {
        if (performanceVariable[i] == "n") {
            if (nominalTransfMethod[i] == "q") {
                if (scaleType[i] == "i") 
                  stop("Quotient transfer method should be used only for ratio data")
                for (j in 1:nrow(data)) {
                  data[j, i] <- min(data[j, i], nomOptValues[i])/max(data[j, 
                    i], nomOptValues[i])
                }
            }
            else {
                data[, i] <- -abs(data[, i] - nomOptValues[i])
                scaleType[i] <- "i"
            }
            performanceVariable[i] <- "s"
        }
    }
    if (patternCoordinates == "manual") {
        pattern <- .manualP(patternManual, data, nomOptValues)
    }
    if (patternCoordinates == "dataBounds") {
        pattern <- rep(0, ncol(data))
        vTypes <- performanceVariable
        for (i in 1:length(performanceVariable)) {
            if (patternType == "lower") {
                if (performanceVariable[i] == "s") {
                  vTypes[i] <- "d"
                }
                if (performanceVariable[i] == "d") {
                  vTypes[i] <- "s"
                }
            }
        }
        for (i in 1:length(vTypes)) {
            if (patternCoordinates == "dataBounds" && vTypes[i] == 
                "s") {
                pattern[i] <- max(data[, i])
            }
            if (patternCoordinates == "dataBounds" && vTypes[i] == 
                "d") {
                pattern[i] <- min(data[, i])
            }
        }
    }
    tdata <- rbind(data, pattern)
    row.names(tdata)[nrow(tdata)] <- "pattern"
    if (sum(scaleType == "i") != 0 && sum(normalization == paste("n", 
        c("6","6a","7","8","9","9a","10","11"), sep = "")) != 0) 
        stop(paste("Normalization ", normalization, " method should not be used for interval data"))
    dataAndPattern <- data.Normalization(rbind(data, pattern), 
        normalization)
    data <- dataAndPattern[-nrow(dataAndPattern), ]
    pattern <- dataAndPattern[nrow(dataAndPattern), ]
    ndata <- rbind(data, pattern)
    row.names(ndata)[nrow(ndata)] <- "pattern"
    cdata <- rbind(pattern, data)
    gdm <- GDM1(cdata, weightsType = weightsType, weights = weights)
    gdm_p <<- as.matrix(gdm)[1, ][-1]
    names(gdm_p) <- row.names(data)
    if (patternType == "upper") {
        sortedgdm_p <- sort(gdm_p)
    }
    if (patternType == "lower") {
        sortedgdm_p <- sort(gdm_p, decreasing = TRUE)
    }
    resul <- list(pdata = pdata, tdata = tdata, data = ndata, 
        distances = gdm_p, sortedDistances = sortedgdm_p)
    resul
}
