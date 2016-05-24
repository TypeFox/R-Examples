"ucedergreen" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
alpha)
{
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useD)) {stop("Not logical useD argument")}
#    if (useD) {stop("Derivatives not available")}
    
    if (missing(alpha)) {stop("'alpha' argument must be specified")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec

    ## Defining the function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        numTerm <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*exp(-1/dose^alpha)
        denTerm <- 1 + exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        parmMat[, 3] - numTerm/denTerm
    }

    ## Defining self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initval <- llogistic()$ssfct(dframe)   
            initval[1] <- -initval[1]
            initval[5] <- 0  # better solution?
    
            return(initval[notFixed])
        }
    }   
    
#    ## Setting the names of the parameters
#    names <- names[notFixed]


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }


    ## Defining derivatives
    
#    ## Constructing a helper function
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
    
    ## Specifying the derivatives    
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t0 <- exp(-1/(dose^alpha))
        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*t0
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1/t3, 
               1 - 1/t3, 
               -t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               -t0/t3 )[, notFixed]
    }
        
    deriv2 <- NULL

    ## Setting limits
#    if (length(lowerc) == numParm) {lowerLimits <- lowerc[notFixed]} else {lowerLimits <- lowerc}
#    if (length(upperc) == numParm) {upperLimits <- upperc[notFixed]} else {upperLimits <- upperc}

    ## Defining the ED function
    edfct <- function(parm, p, lower = 1e-4, upper = 10000, ...)
    {    
        cedergreen(fixed =  fixed, names = names, alpha = alpha)$edfct(parm, 100 - p, lower, upper, ...) 
    }

#    ## Defining the SI function
#    sifct <- function(parm1, parm2, pair, upper = 10000, interval = c(1e-4, 10000))
#    {
#        cedergreen(alpha = alpha)$sifct(parm1, parm2, 100-pair, upper, interval)
#    }    

    ## Finding the maximal hormesis
    maxfct <- function(parm, upper, interval)
    {
       retVal <- cedergreen(fixed =  fixed, names = names, alpha = alpha)$maxfct(parm, upper, interval)
#       retVal[2] <- (parm[2] + parm[3]) - (retVal[2] - parm[2])
       retVal[2] <- (parm[2] + parm[3]) - retVal[2]
              
       return(retVal)
    }

    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], edfct = edfct, maxfct = maxfct,
    name = "ucedergreen",
    text = "U-shaped Cedergreen-Ritz-Streibig", 
    noParm = sum(is.na(fixed)))
    
    class(returnList) <- "UCRS"
    invisible(returnList)
}


"UCRS.4a" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 1, ...))
}

uml3a <- UCRS.4a

"UCRS.4b" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 0.5, ...))
}

uml3b <- UCRS.4b

"UCRS.4c" <-
function(names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = c(names[1], "c", names[2:4]), fixed = c(NA, 0, NA, NA, NA), alpha = 0.25, ...))
}

uml3c <- UCRS.4c

"UCRS.5a" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 1, ...))
}

uml4a <- UCRS.5a

"UCRS.5b" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 0.5, ...))
}

uml4b <- UCRS.5b

"UCRS.5c" <-
function(names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(ucedergreen(names = names, fixed = c(NA, NA, NA, NA, NA), alpha = 0.25, ...))
}

uml4c <- UCRS.5c
