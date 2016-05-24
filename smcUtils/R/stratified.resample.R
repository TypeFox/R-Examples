

stratified.resample = function(weights, num.samples=length(weights), engine="R")
{
    engine=pmatch(engine, c("R","C"))
    n = length(weights)

    switch(engine,
    {
        lbs = seq(0, by=1/num.samples, length=num.samples)
        ubs = lbs+1/num.samples
        u = runif(num.samples, lbs, ubs)
        return(inverse.cdf.weights(weights,u,engine="R"))
    },
    {
        # C implementation
        out = .C("stratified_resample_R", 
                 as.integer(n),
                 as.double(weights),
                 as.integer(num.samples),
                 id = integer(num.samples))
        return(out$id+1)
    })
}

