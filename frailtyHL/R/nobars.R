nobars <-
function (term) 
{
    if (!("|" %in% all.names(term))) 
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

