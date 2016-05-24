"gaussian" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"),
method = c("1", "2", "3", "4"), ssfct = NULL,
fctName, fctText, loge = FALSE)
{
    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Mean function and first derivative
#    y0+a*exp(-0.5*abs((x-x0)/b)^c)
    
    ## deriv(~c+(d-c)*exp(-0.5 * ((dose-e)/b)^f), c("b", "c", "d", "e", "f"), function(dose, b,c,d,e,f){})
    ## deriv(~c+(d-c)*exp(-0.5 * (sqrt(((dose-e)/b)^2))^f), c("b", "c", "d", "e", "f"), function(dose, b,c,d,e,f){})
    fd <- function (dose, b, c, d, e, f)
    {
        .expr1 <- d - c
        .expr3 <- dose - e
        .expr4 <- .expr3/b
        .expr5 <- .expr4^2
        .expr6 <- sqrt(.expr5)
        .expr7 <- .expr6^f
        .expr9 <- exp(-0.5 * .expr7)
        .expr13 <- .expr6^(f - 1)
        .expr18 <- .expr5^-0.5
        .value <- c + .expr1 * .expr9
        .grad <- array(0, c(length(.value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
        .grad[, "b"] <- .expr1 * (.expr9 * (0.5 * (.expr13 * (f * (0.5 * (2 * (.expr3/b^2 * .expr4) * .expr18))))))
        .grad[, "c"] <- 1 - .expr9
        .grad[, "d"] <- .expr9
        .grad[, "e"] <- .expr1 * (.expr9 * (0.5 * (.expr13 * (f * (0.5 * (2 * (1/b * .expr4) * .expr18))))))
        .grad[, "f"] <- -(.expr1 * (.expr9 * (0.5 * (.expr7 * log(.expr6)))))
        attr(.value, "gradient") <- .grad
        .value
    }

    ## Defining the dose-response function
    fct <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4], parmMat[, 5])
    }

    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct  # in case it is explicitly provided
    } else {
        ssfct <- gaussian.ssf(method, fixed)
    }

    ## Defining names
    names <- names[notFixed]

    ##Defining the first and second derivative (in the parameters)
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        attr(fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4], parmMat[, 5]), "gradient")[, notFixed]
    }
    deriv2 <- NULL

    ##Defining the first derivative (in the dose) obtained using
    ## deriv(~c+(d-c)*exp(-0.5 * ((dose-e)/b)^f), "dose", function(dose, b,c,d,e,f){})
    ## deriv(~c+(d-c)*exp(-0.5 * (sqrt(((dose-e)/b)^2))^f), "dose", function(dose, b,c,d,e,f){})
    derivx <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm

        dFct <- function (dose, b, c, d, e, f)
        {
            .expr1 <- d - c
            .expr4 <- (dose - e)/b
            .expr5 <- .expr4^2
            .expr6 <- sqrt(.expr5)
            .expr9 <- exp(-0.5 * .expr6^f)
            .value <- c + .expr1 * .expr9
            .grad <- array(0, c(length(.value), 1L), list(NULL, c("dose")))
            .grad[, "dose"] <- -(.expr1 * (.expr9 * (0.5 * (.expr6^(f - 1) * (f * (0.5 * (2 * (1/b * .expr4) * .expr5^-0.5)))))))
            attr(.value, "gradient") <- .grad
            .value
        }
        attr(dFct(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4], parmMat[, 5]), "gradient")
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
        p <- absToRel(parmVec, abs(respl), type)

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

        ## deriv(~b*(-2*22)^(1 / f)+e, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){})
        ## using "22" insted of log(pProp)
        EDfct <- function (b, c, d, e, f) 
        {
#            .expr2 <- -2 * 22
            .expr2 <- -2 * log(pProp)
            .expr4 <- sign(respl) * .expr2^(1/f)
            .value <- b * .expr4 + e
            .grad <- array(0, c(length(.value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
            .grad[, "b"] <- .expr4
            .grad[, "c"] <- 0
            .grad[, "d"] <- 0
            .grad[, "e"] <- 1
            .grad[, "f"] <- -(b * (.expr4 * (log(.expr2) * (1/f^2))))
            attr(.value, "gradient") <- .grad
            .value
        }
        EDp <- EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5])
        EDder <- attr(EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5]), "gradient")
        return(list(EDp, EDder[notFixed]))
    }

    ## Defining functions returning lower and upper limit and monotonicity
    lowerAs <- pickParm(parmVec, notFixed, 2)
    upperAs <- pickParm(parmVec, notFixed, 3)
    monoton <- NA # monoParm(parmVec, notFixed, 1, -1)

    ## Returning the function with self starter and names
    returnList <-
    list(fct = fct, ssfct = ssfct, names = names, deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Gaussian", fctText),
    noParm = sum(is.na(fixed)), lowerAs = lowerAs, upperAs = upperAs, monoton = monoton,
    fixed = fixed)

    class(returnList) <- "gaussian"
    invisible(returnList)
}


if (FALSE)
{

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

}