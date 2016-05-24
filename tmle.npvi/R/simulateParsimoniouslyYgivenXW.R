simulateParsimoniouslyYgivenXW <- function(T, Y) {
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Validate arguments
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Argument 'T':
    T <- Arguments$getNumerics(T);
    
    ## Argument 'Y':
    Y <- Arguments$getNumerics(Y);

    test <- (min(Y)<=min(T) & max(T)<=max(Y))
    if (!test) {## if parsimonious method fails
        throw("Parsimonious conditional simulation of Y given (X, W) only works when 'theta' takes its values in 'range(Y)'\nPlease use options 'thetamin' and 'thetamax' in the construction of the 'NPVI' object")
    } else {
        simulationSchemes <- getSimulationScheme(T, Y)
        V <- (runif(length(T)) >= simulationSchemes[, "p1"]) + 1
        out <- simulationSchemes[cbind(1:length(T), V)]
    }
    
    return(out)
}

identifyUniqueEntries <- function(T) {
    ## attributes a unique label to each unique entry of T
    sortedT <- sort(T, index.return=TRUE)
    labelST <- cumsum(c(0, diff(sortedT$x)>0))
    labelT <- rep(NA, length(T))
    labelT[sortedT$ix] <- labelST
    return(labelT+1)
}

getSimulationScheme <- function(T, Y) {
    ## Yinf <- sapply(T, function(tt) {
    ##   max(Y[which(Y<tt)])
    ## })
    ## Ysup <- sapply(T, function(tt) {
    ##   min(Y[which(Y>=tt)])
    ## })
    lab <- identifyUniqueEntries(T)
    suT <- sort(unique(T))
    ## the above  labels  are such  that the  i-th  largest value  of 'T'  is
    ## labelled 'i'
    sortedY <- sort(Y, index.return=TRUE)
    index <- findInterval(suT,
                          sortedY$x)
    Yinf.value <- sortedY$x[index]
    Ysup.value <- sortedY$x[index+1]
    Yinf.index <- sortedY$ix[index]
    Ysup.index <- sortedY$ix[index+1]
    prob <- (Ysup.value - suT)/(Ysup.value-Yinf.value)
    out <- cbind(Yinf.index, Ysup.index, prob, 1-prob)
    out <- out[lab, ]
    colnames(out) <- c("Y1", "Y2", "p1", "p2")
    return(out)
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

