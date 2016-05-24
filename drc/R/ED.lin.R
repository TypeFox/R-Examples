"ED.lin" <- function(lmObject, respLev)
{    
    parCoef <- coef(lmObject)
    lparco <- length(parCoef)

#    yVal <- lmObject$"model"[, 1]
    xVal <- lmObject$"model"[, 2]
    fittedVal <- fitted(lmObject)
#    maxDose <- max(xVal)

    decreasing <- ((lparco == 2) && (parCoef[lparco] < 0)) || ((lparco == 3) && (parCoef[lparco] > 0))

#    if (parCoef[lparco] < 0)  # decreasing trend
    if (decreasing)
    {
        cVal <- fittedVal[which.max(xVal)]
        dVal <- fittedVal[which.min(xVal)]
    } else {
        cVal <- fittedVal[which.min(xVal)]
        dVal <- fittedVal[which.max(xVal)]
        
#        respLev <- 100 - respLev
    }
    ## Truncating in case the lower limit is negative
    cVal <- pmax(0, cVal)  
    
#    
#    if (cVal < 0)
#    {
#        cVal <- 0  # as.numeric(polyroot(coef(lmObject)))
#    }
#    print(c(cVal, dVal))

    ## Defining apply() function to handle vector "respLev" arguments
    
    if (lparco == 2)
    {
        if (!decreasing) {respLev <- 100 - respLev}
    
        appFct <- function(respLev)
        {
#            deltaMethod(lmObject, paste("(", cVal, "-b0+", (100 - respLev)/(100), "*(", dVal - cVal, "))/b1", collapse = ""))
            deltaMethod(lmObject, paste("(", cVal, "-b0+", (100 - respLev)/(100), "*(", dVal - cVal, "))/b1", collapse = ""),        
            parameterNames=c("b0", "b1"))
        }
    } 
    if (lparco == 3)
    {
        if (parCoef[3] < 0) {respLev <- 100 - respLev}
     
        print(c(max(xVal), (-parCoef[2] / (2*parCoef[3]))))

        ## Deciding which leg of parabola 
        if ((-parCoef[2] / (2*parCoef[3])) > max(xVal) && (parCoef[3] < 0))
        {
            signVal <- 1
        }

        if ((-parCoef[2] / (2*parCoef[3])) > max(xVal) && (parCoef[3] < 0))
        {
            signVal <- 1
        }
        
        
         else {
            signVal <- -1 
        }
        
        ## Deciding whether the parabola is a cap or a cup
        if (parCoef[3] < 0)
        {
            decreasing <- 1
        } else {
            decreasing <- -1
        }
        signVal <- signVal * decreasing
        
        
        
#        print(paste("(-b1+", signVal, "*sqrt(b1*b1 - 4*b2*(b0-", cVal + ((100 - respLev)/100) * (dVal - cVal), ")))/(2*b2)", collapse = ""))
        appFct <- function(respLev)
        {
            deltaMethod(lmObject, paste("(-b1+", signVal, "*sqrt(b1*b1 - 4*b2*(b0-", cVal + ((100 - respLev)/100) * (dVal - cVal), ")))/(2*b2)", 
            collapse = ""))
        }
    }
    t(sapply(respLev, appFct))
}
