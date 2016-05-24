jitterSteps <-
function (dataFrame) 
{
    iSample   <- as.numeric(unlist(dataFrame["Sample"]))
    Z         <- unlist(dataFrame["Z"])
    delta     <- max(Z) - min(Z)
    nrows     <- length(Z)
    Response  <- unlist(dataFrame["Response"])
    o         <- order(iSample, Z, Response)
    do        <- dataFrame[o, ]
    do1       <- do[-1, ]
    do2       <- do[-nrows, ]
    idResp    <- c(FALSE, do1["Response"] == do2["Response"])
    idSamp    <- c(FALSE, do1["Dilution"] == do2["Dilution"])
    x         <- NULL
    counter   <- 0
    for (i in 1:nrows) {
        counter <- ifelse((idSamp & idResp)[i], counter + 1, 0)
        x     <- c(x, counter)
    }
    indexOfTreatment <- iSample[o]
    meanIndex        <- mean(indexOfTreatment)
    Zjitter   <- unlist(dataFrame["Z"])[o] +
        delta * ((indexOfTreatment - meanIndex) * 0.05 + (x - 1) * 0.025)
    dataOrder <- cbind(do, indexOfTreatment, dX = x, Zjitter)
    q         <- order(dataOrder["I"])
    dataQ     <- dataOrder[q, ]
    return(dataQ)
}
