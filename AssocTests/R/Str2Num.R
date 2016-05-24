## transite a string into a serial of numbers
Str2Num <- function(S)
{
    y <- as.numeric(strsplit(S, '')[[1]])

    return(y)
}
