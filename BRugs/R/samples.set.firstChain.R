"samplesSetFirstChain" <-
function(first)
#   Set the firstChain field
{
    if(!is.numeric(first))
        stop("first ", "must be numeric")
    first <- as.integer(first)
    if(!(first %in% 1:getNumChains()))
        stop("it is required to have 1 <= first <= nchains")
    options("BRugsSamplesFirstChain" = first)
}
