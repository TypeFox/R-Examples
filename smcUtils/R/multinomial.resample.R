
multinomial.resample = function(weights, num.samples=length(weights), engine="R")
{
    check.weights(weights, log=F, normalized=T)
    stopifnot(num.samples>0)

    engine=pmatch(engine, c("R","C"))
    n = length(weights)

    switch(engine,
    {
        # R implementation
        # sample does not perform the same as inverse.cdf method
        # so to be consistent with C, use inverse.cdf
        #return(sort(sample(n, num.samples, replace = TRUE, prob = weights)))

        return(inverse.cdf.weights(weights,runif(num.samples),engine="R"))
    },
    {
        # C implementation
        out = .C("multinomial_resample_R", 
                 as.integer(n),
                 as.double(weights),
                 as.integer(num.samples),
                 id = integer(num.samples))
        return(out$id+1)
    })
}


