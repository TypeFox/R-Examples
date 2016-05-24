"mchain" <-
function (x, Pi, delta, nonstat = TRUE)
{
    y <- c(list(x=x, Pi=Pi, delta=delta, nonstat=TRUE))
    class(y) <- "mchain"
    return(y)
}

