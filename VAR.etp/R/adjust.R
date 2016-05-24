adjust <-
function(b,btem1,p)
{
    bias <- b - btem1
    bs1 <- b - bias
    delta3 <- min(Mod(polyroot(c(1,-bs1[1:p]))))
    if(delta3 > 1)
        return(bs1)
    delta1 <- 1
    while(delta3 <= 1)
    {
        delta1 <- delta1-0.01
        bias <- delta1*bias
        bs2 <- b - bias
        if (sum(as.numeric(is.infinite(bs2))) > 0) 
        {bs2 <- bs1; break}
        delta3 <- min(Mod(polyroot(c(1,-bs2[1:p]))))
    }
    return(bs2)
}
