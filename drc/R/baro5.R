"baro5" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b1", "b2", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
#    if (useDer) {stop("Derivatives not available")}

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
    
#        c <- 2*parmMat[, 3]*parmMat[, 5]/abs(parmMat[, 3]+parmMat[, 5])
#        f <- 1/(1+exp(-c*(log(parmMat[, 4]) - log(dose))))
#        g <- exp(parmMat[, 3]*(log(parmMat[, 4]) - log(dose)))
#        h <- exp(parmMat[, 5]*(log(parmMat[, 4]) - log(dose)))
#        parmMat[, 1]+((parmMat[,2]-parmMat[,1])/(1+f*g+(1-f)*h))

        c <- 2*parmMat[, 1]*parmMat[, 2]/abs(parmMat[, 1]+parmMat[, 2])

        tempVal <- log(dose) - log(parmMat[, 5])
        f <- 1/(1+exp(c*tempVal))
        g <- exp(parmMat[, 1]*tempVal)
        h <- exp(parmMat[, 2]*tempVal)
        parmMat[, 3]+((parmMat[,4]-parmMat[,3])/(1+f*g+(1-f)*h))

    }

#    ## Defining value for control measurements (dose=0)
#    confct <- function(drcSign)
#    {
#        if (drcSign>0) {conPos <- 1} else {conPos <- 2}
#        confct2 <- function(parm)
#        { 
#            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
#            parmMat[, notFixed] <- parm
#            parmMat[, conPos]
#        }
#        return(list(pos=conPos, fct=confct2))
#    }
#
#    ## Defining flag to indicate if more general ANOVA model is available as alternative
#    anovaYes <- TRUE

    ## Defining self starter function
if (FALSE)
{    
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)

        startVal[4] <- max(resp3)+0.001  # the d parameter
        startVal[3] <- min(resp3)-0.001  # the c parameter
#        startVal[!notFixed] <- fixed[!notFixed] 
        
        if (length(unique(dose2))==1) {return((c(NA, NA, NA, startVal[4], NA))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        logitTrans <- log((startVal[4]-resp3)/(resp3-startVal[3]+0.001))  # 0.001 to avoid 0 in the denominator
        logitFit <- lm(logitTrans~log(dose3))
        startVal[5] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter
        startVal[2] <- startVal[1]

        return(startVal[notFixed])
    }
}
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {   
        ssfct <- function(dframe)
        {
            initval <- (llogistic()$ssfct(dframe))[c(1, 1:4)]   
    
            return(initval[notFixed])
        }        
    }
   
    ## Defining names
    names <- names[notFixed]

#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[5])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:5]))
#    } else {
#        scaleInd <- NULL
#    }

    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL

#    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}

    ## Defining the ED function
    edfct <- NULL

    ## Defining the SI function
    sifct <- NULL    
    
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct=edfct, sifct=sifct,
    name = "baro5",
    text = "Baroflex", 
    noParm = sum(is.na(fixed)))

    class(returnList) <- "baro5"
    invisible(returnList)
}
