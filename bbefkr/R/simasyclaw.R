simasyclaw <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(, samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>=0.0 & rdata[i]<31/62)
        {
            sim[i] = rnorm(1)
        }
        if(rdata[i]>=31/62 & rdata[i]<47/62)
        {
            sim[i] = rnorm(1, mean = -3/2, sd = 4/10)
        }
        if(rdata[i]>=47/62 & rdata[i]<55/62)
        {
            sim[i] = rnorm(1, mean = -1/2, sd = 2/10)
        }
        if(rdata[i]>=55/62 & rdata[i]<59/62)
        {
            sim[i] = rnorm(1, mean = 1/2, sd = 1/10)
        }
        if(rdata[i]>=59/62 & rdata[i]<61/62)
        {
            sim[i] = rnorm(1, mean = 3/2, sd = 1/20)
        }
        if(rdata[i]>=61/62 & rdata[i]<=1)
        {
            sim[i] = rnorm(1, mean = 5/2, sd = 1/40)
        }
    }
    return(sim)
}
