"braincousens" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText)
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
        
        parmMat[,2]+(parmMat[,3]+parmMat[,5]*dose-parmMat[,2])/(1+exp(parmMat[,1]*(log(dose)-log(parmMat[,4]))))
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
##    anovaYes <- TRUE
#    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
#    if (is.na(binVar)) {binVar <- FALSE}
#    if (!binVar) {binVar <- NULL}
#    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining the self starter function
if (FALSE)
{
    ssfct <- function(dataFra)
    {
        dose2 <- dataFra[,1]
        resp3 <- dataFra[,2]

        startVal <- rep(0, numParm)
        startVal[3] <- max(resp3)+0.001  # the d parameter
#        startVal[3] <- resp3[which.min(dose2)]        
        startVal[2] <- min(resp3)-0.001  # the c parameter
        startVal[5] <- 0  # better choice may be possible!
#        startVal[5] <- max(resp3) + 0.001 - startVal[3]
#        startVal[!notFixed] <- fixed[!notFixed] 

        if (length(unique(dose2))==1) {return((c(NA, NA, startVal[3], NA, NA))[notFixed])}  # only estimate of upper limit if a single unique dose value 

        indexT2 <- (dose2>0)
        if (!any(indexT2)) {return((rep(NA, numParm))[notFixed])}  # for negative dose value
        dose3 <- dose2[indexT2]
        resp3 <- resp3[indexT2]

        logitTrans <- log((startVal[3]-resp3)/(resp3-startVal[2] + 0.001))  # 0.001 to avoid 0 in the denominator
        logitFit <- lm(logitTrans~log(dose3))
        startVal[4] <- exp((-coef(logitFit)[1]/coef(logitFit)[2]))  # the e parameter
        startVal[1] <- coef(logitFit)[2]  # the b parameter

        return(startVal[notFixed])
    }
} 
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
#        ssfct <- braincousens.ssf(method, fixed)
        ssfct <- function(dframe)
        {
            initval <- llogistic()$ssfct(dframe)   
            initval[5] <- 0
    
            return(initval[notFixed])
        }           
    }
   
    ## Defining names
    names <- names[notFixed]


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }

#    ## Constructing a helper function
#    ##  also used in 'llogistic' and 'llogistic2'
#    xlogx <- function(x, p)
#    {
#        lv <- (x < 1e-12)
#        nlv <- !lv
#        
#        rv <- rep(0, length(x))
#        
#        xlv <- x[lv] 
#        rv[lv] <- log(xlv^(xlv^p[lv]))
#        
#        xnlv <- x[nlv]
#        rv[nlv] <- (xnlv^p[nlv])*log(xnlv)
#    
#        rv
#    }
    
    ## Defining derivatives
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*dose
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1 - 1/t3, 
               1/t3, 
               t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               dose/t3 )[, notFixed]
    }
        
    deriv2 <- NULL


    ## Limits
#    if (length(lowerc)==numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc)==numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}


    ## Defining the ED function
    edfct <- function(parm, p, lower = 1e-3, upper = 1000, ...)
    {
#        if (is.missing(upper)) {upper <- 1000}
        interval <- c(lower, upper)     
     
        parmVec[notFixed] <- parm

        tempVal <- (100-p)/100

        helpEqn <- function(dose) 
        {
            expVal <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            parmVec[5]*(1+expVal*(1-parmVec[1]))-(parmVec[3]-parmVec[2])*expVal*parmVec[1]/dose
        }
        maxAt <- uniroot(helpEqn, interval)$root
    
        eqn <- function(dose) {tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*dose/(parmVec[3]-parmVec[2]))}
        EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root

        EDdose <- EDp
        tempVal1 <- exp(parmVec[1]*(log(EDdose)-log(parmVec[4])))
        tempVal2 <- parmVec[3]-parmVec[2]
        derParm <- c(tempVal*tempVal1*(log(EDdose)-log(parmVec[4])), -parmVec[5]*EDdose/((tempVal2)^2),
                     parmVec[5]*EDdose/((tempVal2)^2), -tempVal*tempVal1*parmVec[1]/parmVec[4],
                     -EDdose/tempVal2)
        derDose <- tempVal*tempVal1*parmVec[1]/EDdose-parmVec[5]/tempVal2 

        EDder <- derParm/derDose
        
        return(list(EDp, EDder[notFixed]))
    }


#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair, upper=1000, interval=c(1e-3, 1000))
#    {
#        parmVec1[notFixed] <- parm1
#        parmVec2[notFixed] <- parm2
#    
#        tempVal1 <- (100-pair[1])/100
#        tempVal2 <- (100-pair[2])/100
#
#        helpEqn1 <- function(dose) 
#                    {
#                        expVal <- exp(parmVec1[1]*(log(dose)-log(parmVec1[4])))
#                        parmVec1[5]*(1+expVal*(1-parmVec1[1]))-(parmVec1[3]-parmVec1[2])*expVal*parmVec1[1]/dose
#                    }
#        maxAt1 <- uniroot(helpEqn1, interval)$root
#        helpEqn2 <- function(dose) 
#                    {
#                        expVal <- exp(parmVec2[1]*(log(dose)-log(parmVec2[4])))
#                        parmVec2[5]*(1+expVal*(1-parmVec2[1]))-(parmVec2[3]-parmVec2[2])*expVal*parmVec2[1]/dose}
#        maxAt2 <- uniroot(helpEqn2, interval)$root
#    
#        eqn1 <- function(dose) {tempVal1*(1+exp(parmVec1[1]*(log(dose)-log(parmVec1[4]))))-(1+parmVec1[5]*dose/(parmVec1[3]-parmVec1[2]))}
#        EDp1 <- uniroot(eqn1, lower=maxAt1, upper=upper)$root
#        eqn2 <- function(dose) {tempVal2*(1+exp(parmVec2[1]*(log(dose)-log(parmVec2[4]))))-(1+parmVec2[5]*dose/(parmVec2[3]-parmVec2[2]))}
#        EDp2 <- uniroot(eqn2, lower=maxAt2, upper=upper)$root
#
#        SIpair <- EDp1/EDp2
#
#        EDdose1 <- EDp1
#        EDdose2 <- EDp2
#        tempVal11 <- exp(parmVec1[1]*(log(EDdose1)-log(parmVec1[4])))
#        tempVal12 <- parmVec1[3]-parmVec1[2]
#        derParm1 <- c(tempVal1*tempVal11*(log(EDdose1)-log(parmVec1[4])), -parmVec1[5]*EDdose1/((tempVal12)^2),
#                      parmVec1[5]*EDdose1/((tempVal12)^2), -tempVal1*tempVal11*parmVec1[1]/parmVec1[4],
#                      -EDdose1/tempVal12)
#        derDose1 <- tempVal1*tempVal11*parmVec1[1]/EDdose1-parmVec1[5]/tempVal12 
#
#        SIder1 <- (derParm1/derDose1)/EDp2
#
#        tempVal21 <- exp(parmVec2[1]*(log(EDdose2)-log(parmVec2[4])))
#        tempVal22 <- parmVec2[3]-parmVec2[2]
#        derParm2 <- c(tempVal2*tempVal21*(log(EDdose2)-log(parmVec2[4])), -parmVec2[5]*EDdose2/((tempVal22)^2),
#                      parmVec2[5]*EDdose2/((tempVal22)^2), -tempVal2*tempVal21*parmVec2[1]/parmVec2[4],
#                      -EDdose2/tempVal22)
#        derDose2 <- tempVal2*tempVal21*parmVec2[1]/EDdose2-parmVec2[5]/tempVal22 
#
#        SIder2 <- (derParm2/derDose2)*(-EDp1/(EDp2*EDp2))
#
#        return(list(SIpair, SIder1[notFixed], SIder2[notFixed]))
#    }


    ## Finding the maximal hormesis
    maxfct <- function(parm, lower = 1e-3, upper = 1000)
    {
#        if (is.null(upper)) {upper <- 1000}
#        if (is.null(interval)) {interval <- c(1e-3, 1000)}     
#        alpha <- 0.5
        parmVec[notFixed] <- parm
        if (parmVec[1]<1) {stop("Brain-Cousens model with b<1 not meaningful")}
        if (parmVec[5]<0) {stop("Brain-Cousens model with f<0 not meaningful")}
        
        optfct <- function(t)
        {
            expTerm1 <- parmVec[5]*t
            expTerm2 <- exp(parmVec[1]*(log(t)-log(parmVec[4])))
            
            return(parmVec[5]*(1+expTerm2)-(parmVec[3]-parmVec[2]+expTerm1)*expTerm2*parmVec[1]/t)
        }
    
        ED1 <- edfct(parm, 1, lower, upper)[[1]]
               
        doseVec <- exp(seq(log(1e-6), log(ED1), length = 100))
#        print((doseVec[optfct(doseVec)>0])[1])

        maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec)>0])[1], ED1))$root
        return(c(maxDose, fct(maxDose, matrix(parm, 1, length(names)))))
    }


    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, 
#    lowerc = lowerLimits, upperc = upperLimits, confct = confct, 
    edfct = edfct, maxfct = maxfct, 
#    scaleInd = scaleInd, anovaYes = anovaYes,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Brain-Cousens (hormesis)", fctText),
    noParm = sum(is.na(fixed)))

#    returnList <- switch(return, "fct+ss" = list(fct,ssfct,names),
#                                 "fct+ss+der" = list(fct,ssfct,names,deriv1,deriv2),
#                                 "ED" = list(edparm, edfct),
#                                 "SI" = list(siparm, sifct))

    class(returnList) <- "braincousens"
    invisible(returnList)
}


"BC.4" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}
 
    return(braincousens(names=c(names[1], "c", names[2:4]), fixed = c(fixed[1], 0, fixed[2:4]),
    fctName = as.character(match.call()[[1]]),
    fctText = "Brain-Cousens (hormesis) with lower limit fixed at 0", ...))
}

bcl3 <- BC.4

"BC.5" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(braincousens(names = names, fixed = fixed,
    fctName = as.character(match.call()[[1]]), ...))
}

bcl4 <- BC.5
