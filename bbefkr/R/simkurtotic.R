simkurtotic <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:length(rdata))
    {
        if(rdata[i]>=0 & rdata[i]<2/3)
        {
            sim[i] = rnorm(1)
        }
        else
        {
            sim[i] = rnorm(1, mean = 0, sd = 1/10)
        }
    }
    return(sim)
}
