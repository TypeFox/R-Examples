"weibull2x" <- function(
fixed = rep(NA, 5), names = c("b", "c", "d", "e", "t0"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    
    
    if (!is.na(fixed[5])) {stop("The lag time cannot be fixed")}
    
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
    
#        parmMat[,2] + (parmMat[,3] - parmMat[,2]) * (1 - exp(-exp(parmMat[,1] *(log(dose - t0) - log(parmMat[,4])))))        
        ifelse(dose > parmMat[,5], parmMat[,2] + (parmMat[,3] - parmMat[,2]) * (1 - exp(-exp(parmMat[,1] *(log(dose - parmMat[,5]) - log(parmMat[,4]))))), parmMat[,2]) 
    }


    ## Defining the self starter function  
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- function(dataFrame)
        {
            c((weibull2.ssf(method, head(fixed, -1)))(dataFrame), fivenum(dataFrame[, 1])[2])
        }
    }       

   
    ## Defining names
    w2.names <- names[notFixed]

    
    ## Defining derivatives
    deriv1 <- NULL
    deriv2 <- NULL
    derivx <- NULL


    ## Defining the ED function
    edfct <- function(parm, p, reference, type, ...)
    {   
        parmVec[notFixed] <- parm
        
        if ( (parmVec[1] > 0) && (reference == "control") ) 
        {
            p <- 100 - p
            reference <- "upper"  # to avoid resetting of p in 'weibull1'    
        }
        if ( (parmVec[1] < 0) && (reference == "control") ) 
        {
            p <- 100 - p
        }
                
        weibull1(fixed, names)$edfct(head(parm, -1), 100 - p, reference, type, ...) + parmVec[5]
        # not accounting completely for the uncertainty due to estimation of the lag time
    }
    ## using head( , -1) to remove the lag-time parameter


    ## Returning list of functions and values
    returnList <-
    list(fct = fct, ssfct = ssfct, names = w2.names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct, 
#    list(fct=fct, confct=confct, ssfct=ssfct, names=w2.names, deriv1=deriv1, deriv2=deriv2, 
#    lowerc=lowerLimits, upperc=upperLimits, edfct=edfct, anovaYes=anovaYes,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Weibull (type 2)", fctText),     
    noParm = sum(is.na(fixed)),
    fixed = fixed)

    class(returnList) <- "Weibull-2"
    invisible(returnList)
}


"W2x.3" <-
function(fixed = c(NA, NA, NA), names = c("d", "e", "t0"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2x(fixed = c(1, 0, fixed[1:3]), names = c("b", "c", names[1:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Weibull (type 2)"), ...))
}            

"W2x.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("c", "d", "e", "t0"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2x(fixed = c(1, fixed[1:4]), names = c("b", names[1:4]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Weibull (type 2)"), ...))
}
