"residuals.mmglm0" <-
function (object, ...) 
{
    object <- as.dthmm(object)
    return(residuals.dthmm(object))
}

