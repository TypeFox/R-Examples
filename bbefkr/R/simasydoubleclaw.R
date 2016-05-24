simasydoubleclaw <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>=0.0 & rdata[i]<138/300)
        {
            sim[i] = rnorm(1, mean = -1, sd = 2/3)
        }
        if(rdata[i]>=138/300 & rdata[i]<276/300)
        {
            sim[i] = rnorm(1, mean = 1, sd = 2/3)
        }
        if(rdata[i]>=276/300 & rdata[i]<277/300)
        {
            sim[i] = rnorm(1, mean = -1/2, sd = 1/100)
        }
        if(rdata[i]>=277/300 & rdata[i]<278/300)
        {
            sim[i] = rnorm(1, mean = -1, sd = 1/100)
        }
        if(rdata[i]>=278/300 & rdata[i]<279/300)
        {
            sim[i] = rnorm(1, mean = -3/2, sd = 1/100)
        }
        if(rdata[i]>=279/300 & rdata[i]<286/300)
        {
            sim[i] = rnorm(1, mean = 1/2, sd = 7/100)
        }
        if(rdata[i]>=286/300 & rdata[i]<293/300)
        {
            sim[i] = rnorm(1, mean = 1, sd = 7/100)
        }
        if(rdata[i]>=293/300 & rdata[i]<1)
        {
            sim[i] = rnorm(1, mean = 3/2, sd = 7/100)
        }
    }
    return(sim)
}
