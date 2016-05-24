setsetequal <-
function (x, y) 
{
    if (!is.list(x) && !is.list(y)) 
        stop("Arguments must be lists")
    all(c(setmatch(x, y, 0L) > 0L, setmatch(y, x, 0L) > 0L))
}
