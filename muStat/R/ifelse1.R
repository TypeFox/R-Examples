`ifelse1` <-
function(test, x, y, ...)
    if(test) x else if(missing(..1)) y else ifelse1(y, ...)

