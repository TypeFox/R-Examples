setsetdiff <-
function (x, y) 
{
    if (!is.list(x) && !is.list(y)) 
        stop("Arguments must be lists")
    x[match(x, y, 0L) == 0L]
}
