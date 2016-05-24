simtrimodal <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>=0 & rdata[i]<0.45)
        {
            sim[i] = rnorm(1, mean = -6/5, sd = 3/5)
        }
        if(rdata[i]>=0.45 & rdata[i]<0.9)
        {
            sim[i] = rnorm(1, mean = 6/5, sd = 3/5)
        }
        if(rdata[i]>=0.9 & rdata[i]<1)
        {
            sim[i] = rnorm(1, mean = 0, sd = 1/4)
        }
    }	
    return(sim)
}
