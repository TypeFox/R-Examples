predict.singleton <-
function(object, ...)
{
    if (!is.null(levels(object))) {
    	return(levels(object))
    } else return(object)
}
