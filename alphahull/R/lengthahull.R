lengthahull <-
function (ahull.arcs) 
{
    gamma <- 2 * ahull.arcs[, "theta"]
    length <- gamma * ahull.arcs[, "r"]
    return(length = sum(length))
}
