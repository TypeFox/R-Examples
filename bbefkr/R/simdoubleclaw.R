simdoubleclaw <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(, samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>=0.0 & rdata[i]<0.49)
        {
            sim[i] = rnorm(1, mean = -1, sd = 2/3)
        }
        if(rdata[i]>=0.49 & rdata[i]<0.98)
        {
            sim[i] = rnorm(1, mean = 1, sd = 2/3)
        }
        if(rdata[i]>=0.98 & rdata[i]<3440/3500)
        {
            sim[i] = rnorm(1, mean = -3/2, sd = 1/100)
        }
        if(rdata[i]>=3440/3500 & rdata[i]<3450/3500)
        {
            sim[i] = rnorm(1, mean = -1, sd = 1/100)
        }
        if(rdata[i]>=3450/3500 & rdata[i]<3460/3500)
        {
            sim[i] = rnorm(1, mean = -0.5, sd = 1/100)
        }
        if(rdata[i]>=3460/3500 & rdata[i]<3470/3500)
        {
            sim[i] = rnorm(1, mean = 0, sd = 1/100)
        }
        if(rdata[i]>=3470/3500 & rdata[i]<3480/3500)
        {
            sim[i] = rnorm(1, mean = 0.5, sd = 1/100)
        }
        if(rdata[i]>=3480/3500 & rdata[i]<3490/3500)
        {
            sim[i] = rnorm(1, mean = 1, sd = 1/100)
        }
        if(rdata[i]>=3490/3500 & rdata[i]<=1)
        {
            sim[i] = rnorm(1, mean = 3/2, sd = 1/100)
        }
    }
    return(sim)
}
