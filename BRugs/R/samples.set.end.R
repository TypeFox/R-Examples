"samplesSetEnd" <-
function(endIt)
#   Set the end field
{
    if(!is.numeric(endIt))
        stop("endIt ", "must be numeric")
    endIt <- as.integer(endIt)
    options("BRugsSamplesEnd" = endIt)
}
