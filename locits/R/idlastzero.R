idlastzero <-
function (v) 
{
    x <- which(v == 0)
    if (length(x) == 0) 
        stop("Can't find a zero")
    max(x)
}
