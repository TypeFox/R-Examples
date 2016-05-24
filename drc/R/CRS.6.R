"CRS.6" <- function(
fixed = c(NA, NA, NA, NA, NA, NA), 
names = c("b", "c", "d", "e", "f", "g"),
method = c("1", "2", "3", "4"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 6
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2] + (parmMat[,3] - parmMat[,2] + parmMat[,5]*exp(-1/(dose^parmMat[,6])))/(1 + exp(parmMat[,1]*(log(dose) - log(parmMat[,4]))))
    }

    ## Defining self starter function
if (FALSE)
{    
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)

#        startVal[3]<-max(resp3)+0.001  # the d parameter
#        startVal[2]<-min(resp3)-0.001  # the c parameter
        startVal[3] <- 1.05 * resp3[which.min(dose2)]
        startVal[2] <- 0.95 * min(resp3)
        
#        startVal[!notFixed] <- fixed[!notFixed] 

        if (length(unique(dose2))==1) {return((c(NA,NA,startVal[3],NA,NA))[notFixed])}

        indexT2<-(dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}
        dose3<-dose2[indexT2]
        resp3<-resp3[indexT2]

        logitTrans<-log((startVal[3]-resp3)/(resp3-startVal[2] + 0.001))  # 0.001 to avoid 0 in the denominator
        logitFit<-lm(logitTrans~log(dose3))
        startVal[4]<-exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1]<-coef(logitFit)[2]  # the b parameter

#        startVal[5] <- 0  # the f parameter
        ## Solving equation at x=e
        startVal[6] <- 0        
        startVal[5] <- (2*(median(resp3) - startVal[2]) - (startVal[3] - startVal[2]))*exp(1/(startVal[4]^startVal[6]))

        return(startVal[notFixed])
    }
}
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {   
        ssfct <- function(dframe)
        {
            initval <- c(llogistic()$ssfct(dframe), 0)   
            initval[5] <- (2*(median(dframe[, 2])-initval[2])-(initval[3]-initval[2]))*exp(1/(initval[4]^initval[6]))
    
            return(initval[notFixed])
        }        
    }
 
    ## Defining names
    names <- names[notFixed]

    ## Specifying the derivatives    
    deriv1 <- NULL
    deriv2 <- NULL


    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function    
    edfct <- NULL
    
    ## Finding the maximal hormesis
    maxfct <- NULL
    
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, # lowerc=lowerLimits, upperc=upperLimits, 
    edfct = edfct, maxfct = maxfct,
    name = "CRS.6",
    text = "Generalised Cedergreen-Ritz-Streibig (hormesis)", 
    noParm = sum(is.na(fixed)))

    class(returnList) <- "cedergreen.extended"
    invisible(returnList)
}

