simdiscretecomb <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>=0.0 & rdata[i]<6/21)
        {
            sim[i] = rnorm(1, mean = -15/7, sd = 2/7)
        }
        if(rdata[i]>=6/21 & rdata[i]<12/21)
        {
            sim[i] = rnorm(1, mean = -3/7, sd = 2/7)
        }
        if(rdata[i]>=12/21 & rdata[i]<18/21)
        {
            sim[i] = rnorm(1, mean = 9/7, sd = 2/7)
        }
        if(rdata[i]>=18/21 & rdata[i]<19/21)
        {
            sim[i] = rnorm(1, mean = 16/7, sd = 1/21)
        }
        if(rdata[i]>=19/21 & rdata[i]<20/21)
        {
            sim[i] = rnorm(1, mean = 18/7, sd = 1/21)
        }
        if(rdata[i]>=20/21 & rdata[i]<=1)
        {
            sim[i] = rnorm(1, mean = 20/7, sd = 1/21)
        }
    }
    return(sim)
}
