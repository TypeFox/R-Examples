"samplesSetThin" <-
function(thin)
#   Set the thin field
{
    if(!is.numeric(thin))
        stop("thin ", "must be numeric")
    options("BRugsSamplesThin" = thin)
}
