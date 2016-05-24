simsmoothcomb <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>=0.0 & rdata[i]<32/63)
        {
            sim[i] = rnorm(1, mean = -31/21, sd = 32/63)
        }
        if(rdata[i]>=32/63 & rdata[i]<48/63)
        {
            sim[i] = rnorm(1, mean = 17/21, sd = 32/126)
        }
        if(rdata[i]>=48/63 & rdata[i]<56/63)
        {
            sim[i] = rnorm(1, mean = 41/21, sd = 32/252)
        }
        if(rdata[i]>=56/63 & rdata[i]<60/63)
        {
            sim[i] = rnorm(1, mean = 53/21, sd = 32/504)
        }
        if(rdata[i]>=60/63 & rdata[i]<62/63)
        {
            sim[i] = rnorm(1, mean = 59/21, sd = 32/1008)
        }
        if(rdata[i]>=62/63 & rdata[i]<1)
        {
            sim[i] = rnorm(1, mean = 62/21, sd = 32/2016)
        }
    }
    return(sim)
}
