"weibull2" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    
    
#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}

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
    
        parmMat[,2] + (parmMat[,3] - parmMat[,2]) * (1 - exp(-exp(parmMat[,1] *(log(dose) - log(parmMat[,4])))))
    }


#    ## Defining value for control measurements (dose=0)
#    confct <- function(drcSign)
#    {
#        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
#        confct2 <- function(parm)
#        { 
#            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
#            parmMat[, notFixed] <- parm
#            parmMat[, conPos]
#        }
#        return(list(pos=conPos, fct=confct2))
#    }
#
#
#    ## Defining flag to indicate if more general ANOVA model
#    anovaYes <- TRUE


    ## Defining the self starter function
if (FALSE)
{
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)
        startVal[3] <- max(resp3)  # +0.001  # the upper bound
        startVal[2] <- min(resp3)  # -0.001  # the lower bound
#        startVal[!notFixed] <- fixed[!notFixed] 

        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA))[notFixed])}  
        # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        loglogTrans <- log(-log((startVal[3] - resp3 + 0.001)/(startVal[3]-startVal[2])))  # 0.001 to avoid 0 as argument to log
        
        loglogFit <- lm(loglogTrans~log(dose3))
        startVal[4] <- exp(-coef(loglogFit)[1]/coef(loglogFit)[2])  # the e parameter
        startVal[1] <- coef(loglogFit)[2]  # the b parameter

        ## Avoiding 0 as start value for lower limit (convergence will fail)
        if ( startVal[2] < 1e-12 ) {startVal[2] <- startVal[3]/10}

        return(startVal[notFixed])
    }
}    
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- weibull2.ssf(method, fixed)
    }       

   
    ## Defining names
    w2.names <- names[notFixed]


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }
    
    
    ## Defining derivatives
    deriv1 <- function(dose, parm)
              {
                  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
                  parmMat[, notFixed] <- parm

                  t1 <- parmMat[, 3] - parmMat[, 2]
                  t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
                  t3 <- exp(-t2)

#                  derMat <- as.matrix(cbind( t1*t3*t2*(log(dose) - log(parmMat[, 4])), 
                  derMat <- as.matrix(cbind( t1*xexplogx(dose/parmMat[, 4], parmMat[, 1]), 
                                             1 - (1 - t3), 
                                             1 - t3, 
                                             -t1*xexpx(dose/parmMat[, 4], parmMat[, 1])*parmMat[, 1]/parmMat[, 4] ))                                             
#                                             -t1*t3*t2*parmMat[, 1]/parmMat[, 4] ))
                  return(derMat[, notFixed])
              }
    deriv2 <- NULL

    derivx <- NULL

    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
#    edfct <- function(parm, p, upper=NULL)  # upper argument not used in 'gompertz'
#    {
#        parmVec[notFixed] <- parm
#    
#        tempVal <- log(-log(p/100))
#        EDp <- exp(tempVal/parmVec[1] + log(parmVec[4]))
#
#        EDder <- EDp*c( -tempVal/(parmVec[1]*parmVec[1]), 0, 0, 1/parmVec[4])
#    
#        return(list(EDp, EDder[notFixed]))
#    }
    edfct <- function(parm, p, reference, type, ...)
    {   
        parmVec[notFixed] <- parm

        p <- absToRel(parmVec, p, type)

        ## Reversing p
#        if (identical(type, "absolute"))
#        {
#            p <- 100 - p
#            type <- "relative"
#        }
        
        if (identical(type, "absolute") && (parmVec[1] > 0) && (reference == "control"))
        {
            p <- 100 - p
        }
               
#        if ( (parmVec[1] > 0) && (reference == "control") ) 
#        {
#            p <- 100 - p
#            reference <- "upper"  # to avoid resetting of p in weibull1() called below    
#        }
#        if ( (parmVec[1] < 0) && (reference == "control") ) 
#        {
#            p <- 100 - p
#        }
    
                
#        weibull1(fixed, names)$edfct(parm, 100 - p, reference, type, ...) 
        weibull1(fixed, names)$edfct(parm, p, reference, "relative", ...) 
    }


    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair)
#    {
#        parmVec1[notFixed] <- parm1
#        parmVec2[notFixed] <- parm2
#
#        tempVal1 <- log(-log(pair[1]/100))
#        tempVal2 <- log(-log(pair[2]/100))
#    
#        SIpair <- exp(tempVal1/parmVec1[1] + log(parmVec1[4]))/exp(tempVal2/parmVec2[1] + log(parmVec2[4]))
#    
#        SIder1 <- SIpair*c(-tempVal1/(parmVec1[1]*parmVec1[1]), 0, 0, 1/parmVec1[4])
#        SIder2 <- SIpair*c(tempVal2/(parmVec2[1]*parmVec2[1]), 0, 0, -1/parmVec2[4])
#    
#        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))
#    }

##    sifct <- function(parm1, parm2, pair)
##    {
##        weibull(lowerc, upperc, fixed, names, scaleDose, useDer)$sifct(parm1, parm2, 100-pair)
##    }    


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


"W2.2" <- function(
upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(fixed[1], 0, upper, fixed[2]), names = c(names[1], "c", "d", names[2]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Weibull (type 2)", upper), ...))
}

"W2.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(fixed[1], 0, fixed[2:3]), names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Weibull (type 2)"), ...))
}

"W2.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(fixed[1:2], upper, fixed[3]), 
    names = c(names[1:2], "d", names[3]), 
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Weibull (type 2)", upper), ...))
}

"W2.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!(length(fixed) == numParm)) {stop("Not correct length of 'fixed' argument")}
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}

    return(weibull2(fixed = fixed, names = names,    
    fctName = as.character(match.call()[[1]]),
    fctText = "Weibull (type 2)", ...))
}

"AR.2" <-
function(fixed = c(NA, NA), names = c("d", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(1, 0, fixed[1:2]), 
    names = c("b", "c", names[1:2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Asymptotic regression"), ...))
}

"AR.3" <-
function(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(weibull2(fixed = c(1, fixed[1:3]), 
    names = c("b", names[1:3]), 
    fctName = as.character(match.call()[[1]]),
    fctText = "Shifted asymptotic regression", ...))
}
