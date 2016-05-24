"bmd" <- function(object, bmr, backg = 0, def = c("excess", "additional", "relative", "hybrid"), interval = c("delta"), 
ma = FALSE, maList = NULL, display = FALSE)
{
    def <- match.arg(def)

    respType <- object$"type"

#    if (identical(def, "excess"))
#    {
#        bmrScaled <- bmr + backg
#    }
#    if (identical(def, "additional"))
#    {
#        bmrScaled <- bmr * (1 - backg) + backg
#    }
    
    if (identical(respType, "binomial"))
    {
        bmrScaled <- switch(def,
        "additional" = bmr * (1 - backg) + backg,
        "excess" = bmr + backg)
        typeVal <- "absolute"
    
    }
    if (identical(respType, "continuous"))
    {
        bmrScaled0 <- 100 * switch(def, 
        "relative" = bmr,
        # scaling done automatically relative to the lower and upper limits
        "hybrid" = sqrt(summary(object)$"resVar") * (qnorm(1 - backg) - qnorm(1 - (backg + bmr))) / (abs(diff(predict(object, data.frame(c(0, Inf)))))))
        
#        print(predict(object, data.frame(c(0, Inf))))
#        print(bmrScaled0)
        
        bmrScaled <- as.numeric(format(bmrScaled0, digits = 3))
#        print(bmrScaled)
        typeVal <- "relative"  # as the relative calculation in ED() uses percentages between 0 and 100 (and not values between 0 and 1)
    }    
    
    if (display)
    {
        cat("Effective relative response level: ", bmrScaled, "\n\n")
    }
    
    ## Use default list for model-averaging 
    ## Note: Any list provided as the "maList" argument is over-written
    if (ma)
    {
        if (identical(object$"type", "binomial"))
        {
            maList <- list(LL.2(), LN.2(), W1.2(), W2.2())
        }
        if (identical(object$"type", "continuous"))
        {
            maList <- list(LL.5(), LN.4(), W1.4(), W2.4() , FPL.4(-1,1), FPL.4(-2,3), FPL.4(-0.5,0.5))
        }
    }
    
    ## Model-averaging or a single parametric fit
    if (!is.null(maList))
    {
        resMat <- maED(object, maList, bmrScaled, interval = "buckland", level = 0.90, type = typeVal, display = display)[, c(1, 3), drop = FALSE]  # linreg = TRUE
    } else {
        resMat <- ED(object, bmrScaled, interval = interval, level = 0.90, type = typeVal, display = display)[, c(1, 3), drop = FALSE]
    }
    colnames(resMat) <- c("BMD", "BMDL")

    if (display) 
    {
        cat("\n\n")
    }
    resMat
}