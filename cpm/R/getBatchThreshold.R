getBatchThreshold <- function(cpmType, alpha, n, lambda=0.3) {
    #verify function arguments
    allowableAlphas <- c(0.1,0.05,0.01,0.005,0.001)
    maxN <- 10000
    
    args <- verifyArguments(cpmType,500,20,lambda)
    if (args$success==FALSE) {return(NA)}
    
    if (!is.element(alpha, allowableAlphas)) {
        print("Error: No pre-computed thresholds available for selected alpha, please see function documentation for supported values")
        return(NA)
    } else {
        alpha <- which(alpha==allowableAlphas)    
    }
    
    if (n > maxN) {
        print(sprintf("Warning: Pre-computed thresholds only available for n<=10000. Using threshold for n=10000 instead"))
        n <- maxN
    }
    
    cpmnames <- c("Student","Bartlett","FET01","FET03","Joint","JointAdjusted","Exponential","ExponentialAdjusted","MW","Mood","LP", "KS","CVM")
    if (args$cpmType=="FET" && lambda==0.1) {args$cpmType <- "FET01"}
    if (args$cpmType=="FET" && lambda==0.3) {args$cpmType <- "FET03"}

    cpmindex <- which(cpmnames==args$cpmType)
    return(approxfun(cpmpackagens, cpmpackageresults[[alpha]][cpmindex,])(n))
}

