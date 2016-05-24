subbars <-
function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name("|")) 
        term[[1]] <- as.name("+")
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

