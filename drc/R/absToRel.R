"absToRel" <- function(parmVec, respl, typeCalc)
{
    ## Converting absolute to relative
    if (typeCalc == "absolute") 
    {
        p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
    } else {  
        p <- respl
    }

    p
}