"samplesSetLastChain" <-
function(last)
#   Set the lastChain field
{
    if(!is.numeric(last))
        stop("last ", "must be numeric")
    last <- as.integer(last)
    if(!(last %in% 1:getNumChains()))
        stop("it is required to have 1 <= last <= nchains")
    options("BRugsSamplesLastChain" = last)
}
