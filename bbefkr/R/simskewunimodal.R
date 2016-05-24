simskewunimodal <-
function(samplesize)
{ 	
    dum = runif(samplesize)
    sim = vector(,samplesize)
    for(i in 1:length(dum))
    {
        if(dum[i]<=0.2)
        {
            sim[i] = rnorm(1)
        }
        if(dum[i]>0.2 & dum[i]<=0.4)
        {
            sim[i] = rnorm(1, mean = 0.5, sd = 2/3)
        }
        if(dum[i]>0.4 & dum[i]<=1.0)
        {
            sim[i] = rnorm(1, mean = 13/12, sd = 5/9)
        }
    }
    return(sim)
}
