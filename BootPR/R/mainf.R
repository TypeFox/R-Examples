mainf <-
function(b,h,p)
{
    mf <- c(1, b[1], rep(0,h-1))
    for(i in 2:h)
    {
        idx <- 1:min(i,p)
        mf[i+1] <- sum(b[idx] * rev(mf[1:i])[idx])
    }
    return(mf)
}
