simskewbimodal <- function(samplesize)
{
    uni = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:samplesize)
    {
        if(uni[i] <=0.75)
        {		
            sim[i] = rnorm(1)
        }
        else
        {
            sim[i] = rnorm(1, mean = 1.5, sd = 1/3)
        }
    }
    return(sim)
}
