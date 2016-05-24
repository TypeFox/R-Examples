model.frame.bma <-
function (formula, ...) 
{
    if (!is.bma(formula)) 
        stop("argument 'formula' needs to be a bma object")
    return(as.data.frame(formula$X.data))
}
