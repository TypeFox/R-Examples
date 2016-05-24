"threephase" <- function(
fixed = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), 
names = c("b1", "c1", "d1", "e1", "b2", "d2", "e2", "b3", "d3", "e3"), 
fctName, fctText)
{   
    ## Checking arguments
    numParm <- 10
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}    
 
    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
#        print("A")
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
#        print(notFixed)
#        print(parm)
#        print(parmMat[, notFixed])
        parmMat[, notFixed] <- parm
#        print("B")
#        LL.4(fixed[1:4])$fct(dose, parmMat[, 1:4]) + LL.3(fixed[5:7])$fct(dose, parmMat[, 5:7])

#        fixed1.4 <- fixed[1:4]
#        fixed5.7 <- fixed[5:7] 
        {LL.4()$fct(dose, parmMat[, 1:4, drop = FALSE]) + LL.3()$fct(dose, parmMat[, 5:7, drop = FALSE]) + LL.3()$fct(dose, parmMat[, 8:10, drop = FALSE])}
    }

    ## Defining self starter function
    ssfct <- function(dframe)
    {  
#        first4 <- drc:::llogistic.ssf(fixed = fixed[1:4])(dframe)  # drc::: not need
#        first4 <- drc:::llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)    
        first4 <- llogistic.ssf(fixed = c(NA, NA, NA, NA, 1))(dframe)    
    
#    print(c(first4[1:2], first4[3]/2, first4[4]/3, first4[1], first4[3]/2, first4[4])[is.na(fixed)])

#        c(first4[1:2], first4[3]/2, first4[4]/3, first4[1], first4[3]/2, first4[4])[is.na(fixed)]
        c(first4[1:2], first4[3]/3, first4[4]/3, first4[1], first4[3]*(2/3), first4[4]*(2/3),
        first4[1], first4[3], first4[4])[is.na(fixed)]
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
    text = ifelse(missing(fctText), "Three-phase", fctText), 
    noParm = sum(is.na(fixed))
    )
    
    class(returnList) <- "three-phase"
    invisible(returnList)
}
