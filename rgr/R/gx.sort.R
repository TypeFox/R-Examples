gx.sort <-
function (x, col = 1, reverse = FALSE) 
{
    ncol <- length(x[1, ])
    if (col > ncol) 
        stop("Column number must be between 1 and", ncol, "\n\n")
    if (reverse) 
        x[rev(order(x[, col])), ]
    else x[order(x[, col]), ]
}
