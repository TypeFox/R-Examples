#checks the arguments passed to the CPM
#returns a list where success==TRUE iff checks passed, and cpmType
#is the standard form of the CPM

verifyArguments <- function(cpmType, ARL0, startup, lambda) {
    implementedCpms <- c("MW","Mood","KS","CVM","Student","Bartlett","LP","FET","Joint","JointAdjusted", "JointHawkins", "Poisson", "Exponential", "ExponentialAdjusted")
    implementedARLs <- c(seq(100,900,by=100),370,seq(1000,9000,by=1000),seq(10000,50000,by=10000))
    implementedLambdas <- c(0.1,0.3)
    
    if (cpmType=="FET" && is.na(lambda)) {
        lambda <- 0.1 #default value
    }
    
    #allow cpmType to be specified as either acronym or full name for convenience
    if (cpmType=="T") {cpmType <- "Student"}   
    else if (cpmType=="F") {cpmType <- "Bartlett"} 
    else if (cpmType=="Lepage") {cpmType <- "LP"}
    else if (cpmType=="Mann-Whitney") {cpmType <- "MW"}
    else if (cpmType=="Kolmogorov-Smirnoff") {cpmType <- "KS"}
    else if (cpmType=="Kolmogorov-Smirnov") {cpmType <- "KS"}
    else if (cpmType=="CvM") {cpmType <- "CVM"}
    else if (cpmType=="Cramer-von-Mises") {cpmType <- "CVM"}
    else if (cpmType=="Cramer-Von-Mises") {cpmType <- "CVM"}
    else if (cpmType=="GLR") {cpmType <- "Joint"}

    if (!is.element(cpmType,implementedCpms)) {
        print("Error: cpmType is not a valid ChangePointModel type, please see function documentation for supported types")
        return(list(success=FALSE,cpmType=cpmType))
    }  
    
    if (!is.na(ARL0) && !is.element(ARL0,implementedARLs)) {
        print("Error: No thresholds available for selected ARL0, please see function documentation for supported values")
        return(list(success=FALSE,cpmType=cpmType))
    }
    
    if (cpmType == "FET" && !is.element(lambda,implementedLambdas)) {
        print("Error: No thresholds available for selected lambda, please see function documentation for supported values")
        return(list(success=FALSE,cpmType=cpmType))
    }
    
    return(list(success=TRUE,cpmType=cpmType))
}

