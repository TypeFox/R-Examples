findplus <-
function (term) 
{
    if (is.numeric(term)) 
        return(0)
    if (!is.language(term)) 
        return(NULL)
    if (length(term) == 1) 
        return(0)
    if (term[[1]] == as.name("|")) 
        return(findplus(term[[2]]))
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("+")) 
        return(1)
    if (term[[1]] == as.name("-")) 
        return(-1)
}

