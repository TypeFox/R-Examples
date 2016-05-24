tm_map.DCorpus <-
function(x, FUN, ...)
{
    ## TODO: shouldn't we check provided function for sanity?

    x$content <- DLapply(x$content, FUN, ..., keep = x$keep)
    x
}
