"NEC" <- function(
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
        doseDiff <- dose - parmMat[, 4]
        cParm + (parmMat[, 3] - cParm) * exp(-parmMat[, 1] * doseDiff * (doseDiff > 0) )
    }

    ## Defining self starter function
    ssfct <- function(dframe)
    {  
        LLinit <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)  # drc::: not needed      
        c(LLinit[1:3], LLinit[3]/3)[is.na(fixed)]
    }     
      
    ##Defining the first and second derivative (in the parameters) 
    deriv1 <- NULL
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- NULL

    ## Defining the ED function
    edfct <- NULL
  
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "NEC", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "NEC"
    invisible(returnList)
}


"NEC.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = c(fixed[1], 0, upper, fixed[2]), 
    names = c(names[1], "c", "d", names[2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("NEC", upper),
    ...) )
}

"NEC.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = c(fixed[1], 0, fixed[2:3]), 
    names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("NEC"), 
    ...) )
}

"NEC.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( NEC(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), ...) )
}
