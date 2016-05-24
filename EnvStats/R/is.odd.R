is.odd <-
function (i) 
{
    if (any(i != trunc(i))) 
        stop("All values of 'i' must be integers")
    (i%%2) == 1
}
