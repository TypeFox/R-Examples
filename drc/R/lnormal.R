"lnormal" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText, loge = FALSE)
{   
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Mean function and first derivative 
    if (!loge)
    {
        ## deriv(~c+(d-c)*pnorm(b*(log(dose)-log(e))), c("b", "c", "d", "e"), function(dose, b,c,d,e){})   
        fd <- function (dose, b, c, d, e) 
        {
            .expr1 <- d - c
            .expr4 <- log(dose) - log(e)
            .expr5 <- b * .expr4
            .expr6 <- pnorm(.expr5)
            .expr9 <- dnorm(.expr5)
            .value <- c + .expr1 * .expr6
            .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
            # needed to handle multiplication of 0 and -infinity
            tempVec <- .expr9 * .expr4
            tempVec[!is.finite(tempVec)] <- 0            
            .grad[, "b"] <- .expr1 * tempVec
#            .grad[, "b"] <- .expr1 * (.expr9 * .expr4)
            .grad[, "c"] <- 1 - .expr6
            .grad[, "d"] <- .expr6
            .grad[, "e"] <- -(.expr1 * (.expr9 * (b * (1/e))))
            attr(.value, "gradient") <- .grad
            .value
        }
        transfe <- log
        backe <- exp
        
    } else {
        ## Modified after
        ## deriv(~c+(d-c)*pnorm(b*(log(dose)-e)), c("b", "c", "d", "e"), function(dose, b,c,d,e){})   
        fd <- function (dose, b, c, d, e) 
        {
            .expr1 <- d - c
            .expr3 <- log(dose) - e
            .expr4 <- b * .expr3
            .expr5 <- pnorm(.expr4)
            .expr8 <- dnorm(.expr4)
            .value <- c + .expr1 * .expr5
            .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
            # needed to handle multiplication of 0 and -infinity            
            tempVec <- .expr8 * .expr3
            tempVec[!is.finite(tempVec)] <- 0
            .grad[, "b"] <- .expr1 * tempVec            
#            .grad[, "b"] <- .expr1 * (.expr8 * .expr3)
            .grad[, "c"] <- 1 - .expr5
            .grad[, "d"] <- .expr5
            .grad[, "e"] <- -(.expr1 * (.expr8 * b))
            attr(.value, "gradient") <- .grad
            .value
        }
        transfe <- function(x){x}
        backe <- function(x){x}                
    }

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

#        parmMat[,2] + (parmMat[,3] - parmMat[,2]) * exp(-exp(parmMat[,1] *(dose - parmMat[,4])))
        fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4])
    }

    ## Defining the self starter function
#if (FALSE)
#{   
#    ssfct <- function(dframe)
#    {
#        x <- dframe[, 1]
#        y <- dframe[, 2]    
#
#        zeroVal <- 1e-12
#        cVal <- 0.99 * ifelse(notFixed[2], min(y), fixed[2])    
#        dVal <- 1.01 * ifelse(notFixed[3], max(y), fixed[3])
#        
#        ## Finding b and e based on linear regression
#        findbe <- function(x, y, 
#        transx = function(x)
#        {
#            xVec <- log(x)
#            xVec[!is.finite(xVec)] <- NA
#            xVec
#       },
#        transy = function(y) 
#        {
#            denomVal <- 1.01 * max(dVal - y)
#            qnorm((dVal - y) / denomVal)             
##            qnorm((dVal - y)/(dVal - cVal))
#        })
#        {
#            transY <- transy(y)  
#            transX <- transx(x)
#
#            lmFit <- lm(transY ~ transX)
#            coefVec <- coef(lmFit)
##            bVal <- coefVec[2]        
#            bVal <- ifelse(notFixed[1], -coefVec[2], fixed[1]) 
##            eVal <- -coefVec[1] / bVal    
#            eVal <- ifelse(notFixed[4], backe(coefVec[1] / bVal), fixed[4]) 
#    
#            return(as.vector(c(bVal, eVal)))
#        }
#        beVec <- findbe(x, y)
#        
#        c(beVec[1], cVal, dVal, beVec[2])[notFixed]
#    }
#}   
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- lnormal.ssf(method, fixed, loge)
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

    ##Defining the first derivative (in the dose) obtained using
    ## deriv(~c+(d-c)*pnorm(b*(log(dose)-log(e))), "dose", function(dose, b,c,d,e){})
    derivx <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        dFct <- function (dose, b, c, d, e) 
        {
            .expr1 <- d - c
            .expr5 <- b * (log(dose) - transfe(e))
            .value <- c + .expr1 * pnorm(.expr5)
            .grad <- array(0, c(length(.value), 1L), list(NULL, c("dose")))
            .grad[, "dose"] <- .expr1 * (dnorm(.expr5) * (b * (1/dose)))
            attr(.value, "gradient") <- .grad
            .value
        }        
        attr(dFct(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4]), "gradient")
    }

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm
#        if (type == "absolute") 
#        {
#            p <- 100*((parmVec[3] - respl)/(parmVec[3] - parmVec[2]))
#        } else {  
#            p <- respl
#        }
#        if ( (parmVec[1] < 0) && (reference == "control") )
#        {
#            p <- 100 - p
#        }
        p <- absToRel(parmVec, respl, type)
    
        ## Reversing p
        if (identical(type, "absolute"))
        {
            p <- 100 - p
        }
        if (identical(type, "relative") && (parmVec[1] < 0) && (reference == "control"))
        {
            p <- 100 - p
        }
    
        pProp <- 1 - (100-p) / 100
#        EDp <- parmVec[4] * exp(qnorm(1-p) / parmVec[1])

        if (!loge)
        {
            ## deriv(~e * exp(22 / b), c("b", "c", "d", "e"), function(b,c,d,e){})
            ## using "22" instead of qnorm(pProp)
            EDfct <- function (b, c, d, e) 
            {
                .expr2 <- exp(qnorm(pProp) / b)
                .value <- e * .expr2
                .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
                .grad[, "b"] <- -(e * (.expr2 * (qnorm(pProp) / (b^2))))
                .grad[, "c"] <- 0
                .grad[, "d"] <- 0
                .grad[, "e"] <- .expr2
                attr(.value, "gradient") <- .grad
                .value
            }
        } else {
#
#            ## Calculating ED on the original scale
#            ## deriv(~exp(e) * exp(22 / b), c("b", "c", "d", "e"), function(b,c,d,e){})
#            ## using "22" instead of qnorm(pProp)
#            EDfct <- function (b, c, d, e) 
#            {
#                .expr1 <- exp(e)
#                .expr3 <- exp(qnorm(pProp) / b)
#                .expr4 <- .expr1 * .expr3
#                .value <- .expr4
#                .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
#                .grad[, "b"] <- -(.expr1 * (.expr3 * (qnorm(pProp) / (b^2))))
#                .grad[, "c"] <- 0
#                .grad[, "d"] <- 0
#                .grad[, "e"] <- .expr4
#                attr(.value, "gradient") <- .grad
#                .value
#            }
            ## Calculating ED on the log scale
            ## deriv(~e + 22 / b, c("b", "c", "d", "e"), function(b,c,d,e){})
            ## using "22" instead of qnorm(pProp)
            EDfct <- function (b, c, d, e) 
            {
                .value <- e + qnorm(pProp)/b
                .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
                .grad[, "b"] <- -(qnorm(pProp)/ (b^2))
                .grad[, "c"] <- 0
                .grad[, "d"] <- 0
                .grad[, "e"] <- 1
                attr(.value, "gradient") <- .grad
                .value
            }
        }
        EDp <- EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4])
        EDder <- attr(EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4]), "gradient")
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
    text = ifelse(missing(fctText), "Log-normal", fctText), 
    noParm = sum(is.na(fixed)), lowerAs = lowerAs, upperAs = upperAs, monoton = monoton,
    fixed = fixed)
    
    class(returnList) <- "log-normal"
    invisible(returnList)
}

"LN.2" <-
function(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
{
    ## Checking arguments
    numParm <- 2
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( lnormal(fixed = c(fixed[1], 0, upper, fixed[2]), 
    names = c(names[1], "c", "d", names[2]), 
    fctName = as.character(match.call()[[1]]), 
    fctText = lowupFixed("Log-normal", upper), ...) )
}

"LN.3" <-
function(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( lnormal(fixed = c(fixed[1], 0, fixed[2:3]), 
    names = c(names[1], "c", names[2:3]),
    fctName = as.character(match.call()[[1]]), 
    fctText = lowFixed("Log-normal"), ...) )
}

"LN.3u" <-
function(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
{
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( lnormal(fixed = c(fixed[1:2], upper, fixed[3]), 
    names = c(names[1:2], "d", names[3]),
    fctName = as.character(match.call()[[1]]),
    fctText = upFixed("Log-normal", upper), 
    ...) )
}

"LN.4" <-
function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return( lnormal(fixed = fixed, names = names,
    fctName = as.character(match.call()[[1]]), ...) )
}

