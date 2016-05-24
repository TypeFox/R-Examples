"gompertzd" <- function(
fixed = c(NA, NA), names = c("a", "b"))
{   
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        innerT1 <- parmMat[, 2]*dose
        innerT2 <- parmMat[, 1]/parmMat[, 2]*(exp(innerT1) - 1)
        parmMat[, 1]*exp(innerT1 - innerT2)
    }

    ## Defining the self starter function
    ssfct <- function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]    
    
        aVal <- max(y)
        bVal <- 1
        
        return(c(aVal, bVal))
    }
   
    ## Defining names
    names <- names[notFixed]
        
    ##Defining the first derivatives (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        help1 <- fct(dose, parm)
        help2 <- exp(parmMat[, 2]*dose)
        help3 <- help2 - 1
        help4 <- parmMat[, 1]/parmMat[, 2]
        
        deriva <- help1*(1/parmMat[, 1] - help3/parmMat[, 2])
        derivb <- help1*(dose + help4*help3/parmMat[, 2] + help4*help2*dose)
        
        cbind(deriva, derivb)[, notFixed]       
    }
        
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
                  
        fct(x, parm)*(parmMat[, 2] - parmMat[, 1]*exp(parmMat[, 2]*x))
    }

    ## Defining the ED function
    edfct <- NULL
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct, 
    name = "gompertzd",
    text = "Gompertz derivative", 
    noParm = sum(is.na(fixed)))
    
    class(returnList) <- "gompertzd"
    invisible(returnList)
}
