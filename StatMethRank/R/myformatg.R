myformatg <- function (x, digits = getOption("digits"), format = paste0("%.", 
    digits, "g"))
{
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector")
    temp <- sprintf(format, x)
    if (is.matrix(x)) 
        matrix(temp, nrow = nrow(x))
    else temp
}