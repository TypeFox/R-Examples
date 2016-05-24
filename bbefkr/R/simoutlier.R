simoutlier <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i] > 0.0 & rdata[i] <= 0.1)
        {
            sim[i] = rnorm(1)
        }
        else
        {
            sim[i] = rnorm(1, 0, sd=0.1)
        }
    }
    return(sim)
}
