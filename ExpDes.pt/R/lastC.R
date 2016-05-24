lastC <-
function(x)
{
    y <- sub(" +$", "", x)
    p1 <- nchar(y)
    cc <- substr(y, p1, p1)
    return(cc)
}
