"hewlett" <- function(
fixed = c(NA, NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f", "g"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
eps = 1e-10
)
{
    ## Checking arguments
    numParm <- 6
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        loge <- -parmMat[, 6] * log((1/parmMat[, 4])^(1/parmMat[, 6]) + (1/parmMat[, 5])^(1/parmMat[, 6]))        
#        loge <- -parm[, 6]*log((1/parm[, 4])^(1/parm[, 6]) + (1/parm[, 5])^(1/parm[, 6]))
## old        loge <- -parm[, 6]*log((parm[, 4])^(1/parm[, 6]) + (parm[, 5])^(1/parm[, 6]))

        retVec <- parmMat[, 2] + (parmMat[, 3] - parmMat[, 2]) / (1 + exp(parmMat[, 1] * (log(dose) - loge)))
        ## Handling the case dose=0 where "loge" may become NaN due to the mixture encoding (pct in glymet)
        zeroInd <- dose < eps
        retVec[zeroInd] <- ifelse(parmMat[zeroInd, 1] < 0, parmMat[zeroInd, 2],  parmMat[zeroInd, 3])
        retVec
    }

    ## Defining self starter function
if (FALSE)
{    
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)

        startVal[3] <- max(resp3)+0.001  # the d parameter
        startVal[2] <- min(resp3)-0.001  # the c parameter
        startVal[5] <- 1  # better choice may be possible!        
#        startVal[!notFixed] <- fixed[!notFixed] 
        
        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        logitTrans <- log((startVal[3]-resp3)/(resp3-startVal[2]+0.001))  # 0.001 to avoid 0 in the denominator
        logitFit <- lm(logitTrans~log(dose3))
        startVal[4] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter
        
        startVal[5] <- startVal[4]
        startVal[6] <- 1

        return(startVal[notFixed])
    }
}
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initval <- c((llogistic()$ssfct(dframe))[c(1:4, 4)], 1)
    
            return(initval[notFixed])
        }        
    }    
    
    ## Defining names
    names <- names[notFixed]

    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL

    ## Defining the ED function
    edfct <- NULL

    ## Defining the SI function
    sifct <- NULL

    ## Scale function
    scaleFct <- function(doseScaling, respScaling)
    {        
        c(1, respScaling, respScaling, doseScaling, doseScaling, 1)[notFixed]
    }

    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
    edfct = edfct, sifct = sifct, scaleFct = scaleFct,
    name = "hewlett",
    text = "Hewlett mixture", 
    noParm = sum(is.na(fixed)))
                       
    class(returnList) <- "Hewlett"
    invisible(returnList)
}
