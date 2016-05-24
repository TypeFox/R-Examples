makeChangePointModel <- function(cpmType,ARL0=500,startup=20,lambda=NA) {
    args <- verifyArguments(cpmType,ARL0,startup,lambda)
    if (args$success==FALSE) {
        return
    }
    cpmType <- args$cpmType
    
    thresholds <- NA
    if (is.na(ARL0)) {
        thresholds <- rep(999999,10000)
    } else {
        thresholds <- loadThresholds(cpmType, ARL0, 10000, lambda)
    }
    
    
    callfn <- parse(text=sprintf("makeChangePointModel%s",cpmType))
    
    if (cpmType != "FET") {
        return(eval(callfn)(thresholds,startup))
    } else {
        return(eval(callfn)(thresholds,startup,lambda))
    }
}
