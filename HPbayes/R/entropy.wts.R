entropy.wts <-
function(w) 
{
    n <- length(w)
    num <- w * log(w)
    den <- log(n)
    ans <- -sum(num/den, na.rm = TRUE)
    return(ans)
}

