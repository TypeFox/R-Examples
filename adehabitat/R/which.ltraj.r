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
        ex <- parse(text = expr)
        coin <- which(eval(ex, envir = x))
        id <- rep(attr(x,"id"),length(coin))
        burst <- rep(attr(x,"burst"),length(coin))
        return(data.frame(id=id,burst=burst,results=coin))
    }))
    return(li)
}## OK
