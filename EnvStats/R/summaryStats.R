summaryStats <-
function (object, ...) 
{
    if (is.logical(object)) 
        summaryStats.logical(object, ...)
    else if (is.character(object)) 
        summaryStats.character(object, ...)
    else UseMethod("summaryStats")
}
