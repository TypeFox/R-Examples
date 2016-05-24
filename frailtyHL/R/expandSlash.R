expandSlash <-
function (bb) 
{
    if (!is.list(bb)) 
        return(expandSlash(list(bb)))
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
            return(lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                bar, list(foo = x[[2]], bar = trm))))
        x
    }))
}

