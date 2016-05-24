"matmin" <-
function (...)
{
    x = cbind(...)
    if(!is.numeric(x))
        stop("Input should by numeric.")
    if (!is.matrix(drop(x)))
        x = t(x)
    x[1:nrow(x) + nrow(x) * (max.col(-x) - 1)]
}
