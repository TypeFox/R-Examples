"gammadr" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
fctName, fctText)
{   
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        
        cParm <- parmMat[, 2]
#        cParm + (parmMat[, 3] - cParm)/((1+exp(parmMat[, 1]*(log(dose)-log(parmMat[, 4]))))^parmMat[, 5])
        cParm + (parmMat[, 3] - cParm) * pgamma(parmMat[, 1] * dose, parmMat[, 4], 1)
    }

    ## Defining the self starter function  
#    if (!is.null(ssfct))
#    {
#        ssfct <- ssfct  # in case it is explicitly provided
#    } else {

        ssfct <- function(dframe)
        {   
            LL4values <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)
            LL4values * c(-1 / (mean(dframe[, 1])), 1, 1, 1)    
        }
#    }     
      
    ##Defining the first and second derivative (in the parameters)
    
    logGamma <- function(x, y)
    {
        lenX <- length(x)
        retVec <- rep(NA, lenX)
        for (i in 1:lenX)
        {
            intFct <- function(t){dgamma(t, y[i], 1) * log(t)}
#            print(ifelse(x[i] < 1e-10, 0, integrate(intFct, 0, x[i])[[1]]))
            retVec[i] <- ifelse(x[i] < 1e-10, 0, integrate(intFct, 0, x[i])[[1]])
        }
        retVec
    }
            
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2]
        t2 <- pgamma(parmMat[, 1] * dose, parmMat[, 4], 1)

        cbind(
        t1 * dgamma(parmMat[, 1] * dose, parmMat[, 4], 1) * parmMat[, 1],
        1 - t2,
        t2,
        t1 * logGamma(parmMat[, 1] * dose, parmMat[, 4])
#        t1 * (parmMat[, 4] - 1) * pgamma(parmMat[, 1] * dose, parmMat[, 4] - 1, 1)
        )[, notFixed]
    }
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        (parmMat[, 3] - parmMat[, 2]) * dgamma(parmMat[, 1] * x, parmMat[, 4], 1) * parmMat[, 1] 
    }

    ## Defining the ED function
    edfct <- NULL
  
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Gamma", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "gamma"
    invisible(returnList)
}
