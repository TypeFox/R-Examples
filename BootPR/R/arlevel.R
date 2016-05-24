arlevel <-
function(b,p)
{
    a <- numeric(0)
    if(p == 1)
        a <- b[1]
    else
    {
        for(i in 1:p)
        {
            if(i == 1)
                a <- c(a, b[1]+b[2])
            else
                a <- c(a, -b[i]+b[i+1])
            if(i == p)
                a[p] <- -b[i]
        }
    }
    a <- c(a,b[(p+1):length(b)])
    return(as.matrix(a))
}
