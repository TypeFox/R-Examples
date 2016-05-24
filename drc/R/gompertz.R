"gompertz" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
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

    ## Function and derivative obtained using 
    ## deriv(~c+(d-c)*exp(-exp(b*(dose-e))), c("b", "c", "d", "e"), function(dose, b,c,d,e){})   
    fd <- function(dose, b, c, d, e) 
    {
        .expr1 <- d - c
        .expr2 <- dose - e
        .expr4 <- exp(b * .expr2)
        .expr6 <- exp(-.expr4)
        .value <- c + .expr1 * .expr6
        .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
        .grad[, "b"] <- -(.expr1 * (.expr6 * (.expr4 * .expr2)))
        .grad[, "c"] <- 1 - .expr6
        .grad[, "d"] <- .expr6
        .grad[, "e"] <- .expr1 * (.expr6 * (.expr4 * b))
        attr(.value, "gradient") <- .grad
        .value
    }   

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

#        parmMat[,2] + (parmMat[,3] - parmMat[,2]) * exp(-exp(parmMat[,1] *(dose - parmMat[,4])))
        fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4])
    }

    ## Defining self starter function
if (FALSE)
{    
    ssfct <- function(dframe)
    {
        x <- dframe[, 1]
        y <- dframe[, 2]    

        cVal <- ifelse(notFixed[2], 0.99*min(y), fixed[2])    
        dVal <- ifelse(notFixed[3], 1.01*max(y), fixed[3])

        ## Finding b and e based on linear regression
        findbe <- function(x, y, 
        transx = function(x){x},
        transy = function(y) {log(-log((dVal-y)/(dVal-cVal)))})
        {
            transY <- transy(y)  
            transX <- transx(x)

            lmFit <- lm(transY ~ transX)
            coefVec <- coef(lmFit)
#            bVal <- coefVec[2]        
            bVal <- ifelse(notFixed[1], coefVec[2], fixed[1]) 
#            eVal <- -coefVec[1] / bVal    
            eVal <- ifelse(notFixed[4], -coefVec[1] / bVal, fixed[4]) 
    
            return(as.vector(c(bVal, eVal)))
        }
        beVec <- findbe(x, y)
        
        c(beVec[1], cVal, dVal, beVec[2])[notFixed]
    }
}
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- gompertz.ssf(method, fixed)     
    }
   
    ## Defining names
    names <- names[notFixed]
        
    ##Defining the first and second derivative (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        attr(fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4]), "gradient")[, notFixed]
    }    
    deriv2 <- NULL

    ##Defining the first derivative (in the dose)
    derivx <- function(x, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        ## Derivative obtained using 
        ## deriv(~c+(d-c)*exp(-exp(b*(dose-e))), c("dose"), function(dose, b,c,d,e){})
        dFct <- function(x, b, c, d, e) 
        {
            .expr1 <- d - c
            .expr4 <- exp(b * (x - e))
            .expr6 <- exp(-.expr4)
            .value <- c + .expr1 * .expr6
            .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
            .grad[, "x"] <- -(.expr1 * (.expr6 * (.expr4 * b)))
            attr(.value, "gradient") <- .grad
            .value
        }    
        attr(dFct(x, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4]), "gradient")
    }

    ## Defining the ED function
    edfct <- function (parm, respl, reference, type, ...) 
    {
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)
        tempVal <- log(-log((100 - p)/100))
        EDp <- tempVal/parmVec[1] + parmVec[4]
        EDder <- EDp * c(-tempVal/(parmVec[1]^2), 0, 0, 1)
        return(list(EDp, EDder[notFixed]))
    }


    ## Defining functions returning lower and upper limit and monotonicity
    lowerAs <- pickParm(parmVec, notFixed, 2)
    upperAs <- pickParm(parmVec, notFixed, 3)
    monoton <- monoParm(parmVec, notFixed, 1, -1)
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Gompertz", fctText), 
    noParm = sum(is.na(fixed)), lowerAs = lowerAs, upperAs = upperAs, monoton = monoton
    )
    
    class(returnList) <- "gompertz"
    invisible(returnList)
}

"G.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(gompertz(fixed = c(fixed[1], 0, upper, fixed[2]), 
    names = c(names[1], "c", "d", names[2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Gompertz", upper), ...))
}

"G.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(gompertz(fixed = c(fixed[1], 0, fixed[2:3]), 
    names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Gompertz"), ...))
}

"G.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(gompertz(fixed = c(fixed[1:2], upper, fixed[3]), 
    names = c(names[1:2], "d", names[3]),
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Gompertz", upper), ...))
}

"G.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(gompertz(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), ...))
}

