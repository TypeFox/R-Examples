simstrongskew <-
function(samplesize)
{
    rdata = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:length(rdata))
    {
        if(rdata[i]>=0 & rdata[i]<1/8)
        {
            sim[i] = rnorm(1)
        }
        if(rdata[i]>=1/8 & rdata[i]<2/8)
        {
            sim[i] = rnorm(1,-1,2/3)
        }
        if(rdata[i]>=2/8 & rdata[i]<3/8)
        {
            sim[i] = rnorm(1,-5/3,4/9)
        }
        if(rdata[i]>=3/8 & rdata[i]<4/8)
        {
            sim[i] = rnorm(1,-19/9,8/27)
        }
        if(rdata[i]>=4/8 & rdata[i]<5/8)
        {
            sim[i] = rnorm(1,-65/27, 16/81)
        }
        if(rdata[i]>=5/8 & rdata[i]<6/8)
        {
            sim[i] = rnorm(1,-211/81,32/243)
        }
        if(rdata[i]>=6/8 & rdata[i]<7/8)
        {
            sim[i] = rnorm(1,-665/243,64/729)
        }
        if(rdata[i]>=7/8 & rdata[i]<=1)
        {
            sim[i] = rnorm(1,-2059/729,128/2187)
        }
    }
    return(sim)
}
