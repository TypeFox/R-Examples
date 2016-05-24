is.bma <-
function (bmao) 
{
    if (any(is.element(class(bmao), c("bma", "bma.fcast", "bma.sar", 
        "oldbma", "bmav0")))) 
        return(TRUE)
    else return(FALSE)
}
