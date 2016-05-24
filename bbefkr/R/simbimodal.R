simbimodal <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i] > 0.0 & rdata[i] <= 0.5)
        {
            sim[i] = rnorm(1, mean = -1, sd = 2/3)
        }
        else
        {
            sim[i] = rnorm(1, mean = 1, sd = 2/3)
        }
    }
    return(sim)
}
