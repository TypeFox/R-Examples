remove.data <- function(datum, i)
{
    if(is.matrix(datum))
        datum <- datum[-i, ]
    else
        datum <- datum[, , -i]
    return(datum)
}