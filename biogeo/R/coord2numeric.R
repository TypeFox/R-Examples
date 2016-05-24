coord2numeric <-
function (xn) 
{
    if (is.factor(xn)) {
        xn <- as.character(xn)
    }
    x1 <- as.numeric(xn)
    return(x1)
}
