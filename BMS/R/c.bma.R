c.bma <-
function (..., recursive = FALSE) 
{
    if (!missing(recursive)) 
        warning("note that argument recursive has no meaning and is retained for compatibility")
    combine_chains(...)
}
