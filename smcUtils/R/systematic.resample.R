
systematic.resample = function(weights, num.samples=length(weights), engine="R")
{
    check.weights(weights, log=F, normalized=T)
    stopifnot(num.samples>0)

    engine=pmatch(engine, c("R","C"))
    n = length(weights)

    switch(engine,
    {
        u = runif(1,0,1/num.samples)+seq(0,by=1/num.samples,length=num.samples)
        return(inverse.cdf.weights(weights,u,engine="R"))
    },
    {
        # C implementation
        out = .C("systematic_resample_R", 
                 as.integer(n),
                 as.double(weights),
                 as.integer(num.samples),
                 id = integer(num.samples))
        return(out$id+1)
    })
}

