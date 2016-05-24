"rdrm" <- function(nosim, fct, mpar, xerror, xpar = 1, yerror = "rnorm", ypar = c(0, 1), 
onlyY = FALSE)
{        
    ## Constructing the predictor values
    if (is.numeric(xerror))
    {
        x <- xerror
    } else {
        evalStr1 <- paste(xerror, "(", paste(xpar, sep = ",", collapse = ","), ")")
        x <- eval(parse(text = evalStr1))
    }
    lenx <- length(x)
    x <- sort(x)
    x <- rep(x, nosim)
    xMat <- matrix(x, nosim, lenx, byrow = TRUE)
    
    ## Constructing the mean dose-response
    meanVec <- fct$fct(x, matrix(mpar, lenx*nosim, length(mpar), byrow = TRUE))
    
    ## Constructing the simulated response values
    if (yerror == "rbinom")
    {
        if (length(ypar) == 1)
        {
            ypar <- rep(ypar, lenx*nosim)
            wMat <- matrix(ypar, nosim, lenx, byrow = TRUE)
        } else {
            wMat <- matrix(ypar, nosim, lenx, byrow = TRUE)
        }
        evalStr2 <- paste(deparse(substitute(yerror)), "(", lenx*nosim, ", ypar, meanVec)")
        errorVec <- eval(parse(text = evalStr2))
        
        yMat <- matrix(errorVec, nosim, lenx, byrow = TRUE)

        ## Returning the simulated curves
        if (onlyY)
        {
            return(list(y = yMat))
        } else {    
            return(list(x = xMat, w = wMat, y = yMat))
        }
    }  else {
        evalStr2 <- paste(yerror, "(", lenx*nosim, ",",
        paste(ypar, sep = ",", collapse = ","), ")")
        errorVec <- eval(parse(text = evalStr2))
        
        yMat <- matrix(meanVec, nosim, lenx, byrow = TRUE) + errorVec

        ## Returning the simulated curves
        if (onlyY)
        {
            return(list(y = yMat))
        } else {
            return(list(x = xMat, y = yMat))
        }
    }
}
