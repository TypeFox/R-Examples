"fplogistic" <- function(
p1, p2, fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
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

    ## Defining model function and first derivative 
    ## deriv(~c+(d-c)/(1+exp(b*(log(dose+1))^p1+e*(log(dose+1))^p2)), c("b", "c", "d", "e"), function(dose, b,c,d,e){})
    fd <- function (dose, b, c, d, e) 
    {
        .expr1 <- d - c
        .expr3 <- log(dose + 1)
        .expr4 <- .expr3^p1
        .expr6 <- .expr3^p2
        .expr9 <- exp(b * .expr4 + e * .expr6)
        .expr10 <- 1 + .expr9
        .expr15 <- .expr10^2
        .expr18 <- 1/.expr10
        .value <- c + .expr1/.expr10
        .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
        .grad[, "b"] <- -(.expr1 * (.expr9 * .expr4)/.expr15)
        .grad[, "c"] <- 1 - .expr18
        .grad[, "d"] <- .expr18
        .grad[, "e"] <- -(.expr1 * (.expr9 * .expr6)/.expr15)
         attr(.value, "gradient") <- .grad
        .value
    }
    
    ## Defining the nonlinear model function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4])
    }

    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
        ssfct <- function(dframe)
        {
            initval <- llogistic()$ssfct(dframe)[1:4]   
            initval[4] <- initval[1]
            initval[1] <- -initval[1]
    
            return(initval[notFixed])
        }        
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

    ## Defining the first derivative (in the dose) obtained using
    ## deriv(~c+(d-c)/(1+exp(b*(log(dose+1))^p1+e*(log(dose+1))^p2)), "dose", function(dose, b,c,d,e){})
    derivx <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        dFct <- function (dose, b, c, d, e) 
        {
            .expr1 <- d - c
            .expr2 <- dose + 1
            .expr3 <- log(.expr2)
            .expr9 <- exp(b * .expr3^p1 + e * .expr3^p2)
            .expr10 <- 1 + .expr9
            .expr15 <- 1/.expr2
            .value <- c + .expr1/.expr10
            .grad <- array(0, c(length(.value), 1L), list(NULL, c("dose")))
            .grad[, "dose"] <- -(.expr1 * (.expr9 * (b * (.expr3^(p1 - 1) * (p1 * .expr15)) + e * (.expr3^(p2 - 1) * (p2 * .expr15))))/.expr10^2)
            attr(.value, "gradient") <- .grad
            .value
        }
        attr(dFct(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4]), "gradient")
    }

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, loged = FALSE, ...)
    {
        parmVec[notFixed] <- parm
        p <- EDhelper(parmVec, respl, reference, type)

        invfp <- function(resp, b, e)
        {
            fct0 <- function(x){resp - (b*(log(x+1)^p1) + e*(log(x+1)^p2))}
            uniroot(fct0, c(0.001, 1000))$root
        }

        EDfct <- function(b, c, d, e) 
        {
            invfp(log((100-p)/p), b, e)
        }
        
#        ## deriv(~b*(log(dose+1)^p1) + e*(log(dose+1)^p2), c("b", "c", "d", "e"), function(dose, b,c,d,e){})
#        ## note: c and d parameters need not be included
#        derfp <- function (dose, b, c, d, e) 
#        {
#            .expr2 <- log(dose + 1)
#            .expr3 <- .expr2^p1
#            .expr5 <- .expr2^p2
#            .value <- b * .expr3 + e * .expr5
#            .grad <- array(0, c(length(.value), 4L), list(NULL, c("b", "c", "d", "e")))
#            .grad[, "b"] <- .expr3
#            .grad[, "c"] <- 0
#            .grad[, "d"] <- 0
#            .grad[, "e"] <- .expr5
#            attr(.value, "gradient") <- .grad
#            .value
#        }

        EDp <- EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4])
        
        logEDp <- log(EDp+1)
        denVal <- parmVec[1] * p1 * (logEDp)^(p1-1) + parmVec[4] * p2 * (logEDp)^(p2-1)
        derVec <- (EDp+1) * c(logEDp^p1, logEDp^p2) / denVal
        EDder <- c(derVec[1], 0, 0, derVec[2])
        if (loged) 
        {
            EDder <- EDder / EDp
            EDp <- log(EDp)
        }
        
#        EDder <- 1 / attr(derfp(EDp, parmVec[1], parmVec[2], parmVec[3], parmVec[4]), "gradient")
#        EDder <- c(EDder[1], 0, 0, EDder[4])
        return(list(EDp, EDder[notFixed]))
    }
    
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct, 
    name = ifelse(missing(fctName), paste(as.character(match.call()[[1]]), "(", p1,",", p2, ")", sep=""), fctName),
    text = ifelse(missing(fctText), "Fractional polynomial", fctText),     
    noParm = sum(is.na(fixed)),
    fixed = fixed)  
    
    class(returnList) <- "fp-logistic"
    invisible(returnList)    
}

"FPL.4" <-
function(p1, p2, fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
{
    ## Checking arguments
    numParm <- 4
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct names argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct length of 'fixed' argument")}

    return(fplogistic(p1, p2, fixed = fixed, names = names,
    fctName = paste(as.character(match.call()[[1]]), "(", p1,",", p2, ")", sep=""), ...) )
}

