dbase <-
function (env, uniq = T) 
{
    tmp <- switch(bool(uniq) + 1, grep, pmatch)
    tmp <- tmp(env, search())
    if (any(is.na(tmp)) || is.null(tmp)) 
        cat("No pattern", dQuote(env), "could be found in the search list.\n")
    tmp
}
