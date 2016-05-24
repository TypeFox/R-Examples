fitinbounds <-
function (v, lower, upper)
{
    blob <- v
    for (st in names(v)) {
        blob[st] <- max(blob[st], lower[st])
        blob[st] <- min(blob[st], upper[st])
    }
    return(blob)
}
