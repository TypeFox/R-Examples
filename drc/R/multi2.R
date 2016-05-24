"multi2" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b1", "b2", "b3", "c", "d"), 
ssfct = NULL, fctName, fctText)
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
    ## deriv(~c+(d-c)*(1 - exp(-b1-b2*dose-b3*(dose*dose))), c("b1", "b2", "b3", "c", "d"), function(dose, b1, b2, b3, c, d){})           
    fd <- function (dose, b1, b2, b3, c, d) 
    {
        .expr1 <- d - c
        .expr5 <- dose * dose
        .expr8 <- exp(-b1 - b2 * dose - b3 * .expr5)
        .expr9 <- 1 - .expr8
        .value <- c + .expr1 * .expr9
        .grad <- array(0, c(length(.value), 5L), list(NULL, c("b1", "b2", "b3", "c", "d")))
        .grad[, "b1"] <- .expr1 * .expr8
        .grad[, "b2"] <- .expr1 * (.expr8 * dose)
        .grad[, "b3"] <- .expr1 * (.expr8 * .expr5)
        .grad[, "c"] <- 1 - .expr9
        .grad[, "d"] <- .expr9
        attr(.value, "gradient") <- .grad
        .value
    }

    ##Defining the first derivative (in the dose) obtained using
    ## deriv(~c+(d-c)*(1 - exp(-b1-b2*dose-b3*(dose*dose))), "dose", function(dose, b1, b2, b3, c, d){})             
    dFct <- function (dose, b1, b2, b3, c, d) 
    {
        .expr1 <- d - c
        .expr8 <- exp(-b1 - b2 * dose - b3 * (dose * dose))
        .value <- c + .expr1 * (1 - .expr8)
        .grad <- array(0, c(length(.value), 1L), list(NULL, c("dose")))
        .grad[, "dose"] <- .expr1 * (.expr8 * (b2 + b3 * (dose + dose)))
        attr(.value, "gradient") <- .grad
        .value
    }

    ## Defining the non-linear function
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

        ssfct <- function(dframe)
        {  
            first4 <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)    
#            c(0, -first4[1], 0, first4[2:3])[is.na(fixed)]
            c(0, -first4[1] / (mean(dframe[, 1]) * 0.7), 0, first4[2:3])[notFixed]            
        }
    }    
   
    ##Defining the first and second derivative (in the parameters) 
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
        attr(fd(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4], parmMat[, 5]), "gradient")[, notFixed]
    }    
    deriv2 <- NULL

    ##Defining the first derivative in the dose
    derivx <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm    
        attr(dFct(dose, parmMat[, 1], parmMat[, 2], parmMat[, 3], parmMat[, 4], parmMat[, 5]), "gradient")
    }

    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, ...)
    {
        parmVec[notFixed] <- parm
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
#        pProp <- log(1 - (100 - p) / 100)
        pProp <- log((100 - p) / 100)
        
        ## deriv(~ (-b2+sqrt(b2*b2-4*b3*(b1+22)))/(2*b3), c("b1", "b2", "b3", "c", "d"), function(b1, b2, b3, c, d){})
        ## using "22" instead of pProp
        EDfct <- function (b1, b2, b3, c, d) 
        {
            .expr3 <- 4 * b3
            .expr4 <- b1 + pProp
            .expr6 <- b2 * b2 - .expr3 * .expr4
            .expr8 <- -b2 + sqrt(.expr6)
            .expr9 <- 2 * b3
            .expr11 <- .expr6^-0.5
            .value <- .expr8/.expr9
            .grad <- array(0, c(length(.value), 5L), list(NULL, c("b1", "b2", "b3", "c", "d")))
            .grad[, "b1"] <- -(0.5 * (.expr3 * .expr11)/.expr9)
            .grad[, "b2"] <- (0.5 * ((b2 + b2) * .expr11) - 1)/.expr9
            .grad[, "b3"] <- -(0.5 * (4 * .expr4 * .expr11)/.expr9 + .expr8 * 2/.expr9^2)
            .grad[, "c"] <- 0
            .grad[, "d"] <- 0
            attr(.value, "gradient") <- .grad
            .value
        }

        EDp <- EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5])
        EDder <- attr(EDfct(parmVec[1], parmVec[2], parmVec[3], parmVec[4], parmVec[5]), "gradient")
        return(list(EDp, EDder[notFixed]))
    }
   
    ## Returning the function with self starter and names
    returnList <- 
    list(fct = fct, ssfct = ssfct, names = names[notFixed], 
    deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, edfct = edfct,
    name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text = ifelse(missing(fctText), "Multistage", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "multistage"
    invisible(returnList)    
}
