
residual.resample = function(weights, num.samples=length(weights), engine="R", rrf="stratified")
{
    rrf = pmatch(rrf,c("stratified","multinomial","systematic"))
    if (is.na(rrf)) stop("No matching residual resampling function.")

    check.weights(weights, log=F, normalized=T)
    stopifnot(num.samples>0)

    engine=pmatch(engine, c("R","C"))
    n = length(weights)

    switch(engine,
    {
        # R implementation
        n.exp.samps = num.samples * weights
        det.reps    = floor(n.exp.samps)
        num.samples   = num.samples - sum(det.reps)
        det.ids     = rep2id(det.reps, engine="R") 
        weights     = renormalize(n.exp.samps-det.reps, engine="R")
        switch(rrf,
            { ran.ids = stratified.resample( weights, num.samples, engine="R") },
            { ran.ids = multinomial.resample(weights, num.samples, engine="R") },
            { ran.ids = systematic.resample( weights, num.samples, engine="R") })
        return(c(det.ids, ran.ids))
    },
    {
        # C implementation
        out = .C("residual_resample_R", 
                 as.integer(n),
                 as.double(weights),
                 as.integer(num.samples),
                 id = integer(num.samples),
                 as.integer(rrf))
        return(out$id+1)
    })
}

