var.rwts <-
function (w) 
{
    n <- length(w)
    ans <- mean((n * w - 1)^2)
    return(ans)
}

