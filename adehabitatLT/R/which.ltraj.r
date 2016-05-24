##################################################################
##
##
## which.ltr verifie une condition et renvoie les ID, bursts et
## indices de localisations qui la remplisse


which.ltraj <- function(ltraj, expr)
{
    if (!inherits(ltraj, "ltraj"))
        stop("ltraj should be of class \"ltraj\"")
    li <- do.call("rbind",lapply(ltraj, function(x) {
        id <- attr(x, "id")
        bu <- attr(x, "burst")
        if (!is.null(attr(x,"infolocs")))
            x <- cbind(x, attr(x,"infolocs"))
        ex <- parse(text = expr)
        coin <- which(eval(ex, envir = x))
        id <- rep(id,length(coin))
        burst <- rep(bu,length(coin))
        return(data.frame(id=id,burst=burst,results=coin))
    }))
    return(li)
}## OK
