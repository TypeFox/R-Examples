simsepbimodal <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:length(rdata))
    {
        if(rdata[i]>=0 & rdata[i]<1/2)
        {
            sim[i] = rnorm(1, mean = -3/2, sd=1/2)
        }
        else
        {
            sim[i] = rnorm(1, mean= 3/2, sd=1/2)
        }
    }
    return(sim)
}
