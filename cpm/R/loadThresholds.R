#loads the h_t threshold sequence from a private package data object


loadThresholds <- function(cpmType, ARL0, desiredLength=5000,lambda=NA, startup=20) {
    #if lambda is passed, append to cpmType name
    if (!is.na(lambda)) {
        lambda <- gsub("\\.","",lambda) #remove the full stop
        cpmType <- sprintf("%s%s",cpmType,lambda)
    }
    
    #load thresholds and expand to the desired length
    str <- sprintf("%sARL%s",cpmType,ARL0)
    thresholds <- c(rep(99999,19),cpmthresholds[[str]][,1])
    thresholds <- thresholds[!is.na(thresholds)]
    len <- length(thresholds)
    if (len< desiredLength) {
        thresholds <- c(thresholds,rep(mean(thresholds[(len-100):len]), desiredLength-len))
    }
    return(thresholds)
}
