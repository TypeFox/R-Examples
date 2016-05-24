simclaw <-
function(samplesize)
{
    rdata = runif(samplesize)
    epsilon = vector(, samplesize)
    for(i in 1:samplesize)
    {
        if(rdata[i]>0.0 & rdata[i]<=0.1)
        {
            epsilon[i] = rnorm(1,-1,0.1)
        }
        if(rdata[i]>0.1 & rdata[i]<=0.2)
        {		
            epsilon[i] = rnorm(1,-0.5,0.1)
        }
        if(rdata[i]>0.2 & rdata[i]<=0.3)
        {
            epsilon[i] = rnorm(1,0,0.1)
        }
        if(rdata[i]>0.3 & rdata[i]<=0.4)
        {
            epsilon[i] = rnorm(1,0.5,0.1)
        }
        if(rdata[i]>0.4 & rdata[i]<=0.5)
        {
            epsilon[i] = rnorm(1,1,0.1)
        }
        if(rdata[i]>0.5)
        {
            epsilon[i] = rnorm(1,0,1)
        }
    }
    return(epsilon)
}
