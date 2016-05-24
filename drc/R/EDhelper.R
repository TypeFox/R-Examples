"EDhelper" <- function(parmVec, respl, reference, typeCalc, cond = TRUE)
{
    ## Converting absolute to relative
    if (typeCalc == "absolute") 
    {
        p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
#        typeCalc <- "relative"
    } else {  
        p <- respl
    }
    ## Swapping p for an increasing curve
    if (cond)
    {
        if ((typeCalc == "relative") && (parmVec[1] < 0) && (reference == "control"))
        {
            p <- 100 - p
        }
    } else {
        if ((typeCalc == "relative") && (reference == "control"))
        {
            p <- 100 - p
        }
    }
    p
}