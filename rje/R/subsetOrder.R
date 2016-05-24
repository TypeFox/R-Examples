subsetOrder <-
function (x, y) 
{
    if (setequal(x, y)) 
        return(0)
    else if (is.subset(x, y)) 
        return(-1)
    else if (is.subset(y, x)) 
        return(1)
    else return(0)
}
