detectChangePoint <- function(x,cpmType, ARL0=500,startup=20,lambda=NA) {
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
    res <-.C("cpmDetectChange",as.character(cpmType),as.double(x),as.integer(length(x)),as.double(thresholds),as.integer(length(thresholds)),as.integer(startup),Ds=as.double(numeric(length(x))),cp=as.integer(0),dt=as.integer(0),lambda=as.double(lambda),PACKAGE="cpm")
    dt <- res[["dt"]] 
    if (dt==0) {dt <- length(x)} #needed if no change points are detected
    return(list(x=x,Ds=res[["Ds"]][1:dt],thresholds=thresholds,changePoint=res[["cp"]], detectionTime = res[["dt"]],changeDetected = res[["dt"]]>0))
}



