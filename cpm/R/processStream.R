processStream <- function(x,cpmType, ARL0=500 ,startup=20,lambda=NA) {
    args <- verifyArguments(cpmType,ARL0,startup,lambda)
    if (args$success==FALSE) {
        return
    }
    cpmType <- args$cpmType
    
    thresholds <- NA
    if (is.na(ARL0)) {
        thresholds <- rep(999999,length(x))
    } else {
        thresholds <- loadThresholds(cpmType, ARL0, length(x), lambda)
    }
    
    if(is.na(lambda)){lambda<-0}
    res <-.C("cpmProcessStream",as.character(cpmType),as.double(x),as.integer(length(x)),as.double(thresholds),as.integer(length(thresholds)),as.integer(startup),Ds=as.double(numeric(length(x))),cps=as.integer(numeric(length(x))),dts=as.integer(numeric(length(x))),numChanges=as.integer(0),lambda=as.double(lambda),PACKAGE="cpm")
    cps <- numeric()
    dts <- numeric()
    numChanges <- res[["numChanges"]]
    if (numChanges > 0) {
        cps <- res[["cps"]][1:numChanges]
        dts <- res[["dts"]][1:numChanges]
    }
    return(list(x=x,changePoints=cps, detectionTimes = dts))
}
