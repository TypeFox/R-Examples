detectChangePointBatch <- function(x,cpmType, alpha=0.05, lambda=NA) {
    if (length(x) < 5) {
        print("Error: Batch change detectoin requires a sequence of length 5 or greater")
        return
    }
    args <- verifyArguments(cpmType,500,20,lambda)
    if (args$success==FALSE) {
        return
    }
    
    cpmType <- args$cpmType
       
    if(is.na(lambda)){lambda<-0}    
    
    threshold <- NA
    if (!is.na(alpha)) {
    	threshold <- getBatchThreshold(cpmType,alpha,length(x),lambda)
    }
  
    res <-.C("cpmDetectChangeBatch",as.character(cpmType),as.double(x),as.integer(length(x)),Ds=as.double(numeric(length(x))),lambda=as.double(lambda),PACKAGE="cpm")

	changeDetected <- NA; 
    changePoint <- which.max(res[["Ds"]])
	if (!is.na(alpha)) {
    	changeDetected <- max(res[["Ds"]]) > threshold
        if (changeDetected==FALSE) {changePoint<-0}
    }
    return(list(x=x,Ds=res[["Ds"]],changePoint=changePoint, changeDetected=changeDetected, alpha=alpha,threshold=threshold))
}



