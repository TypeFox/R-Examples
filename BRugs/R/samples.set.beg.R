"samplesSetBeg" <-
function(begIt)
#   Set the beg field
{
    if(!is.numeric(begIt))
        stop("begIt ", "must be numeric")
    begIt <- as.integer(begIt)
    options("BRugsSamplesBeg" = begIt)
}
