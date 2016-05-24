
renormalize = function(weights, log=FALSE, engine="R")
{   
    if (!log && any(weights<0)) stop("log=F but negative weights exist.")

    engine=pmatch(engine, c("R","C"))

    switch(engine,
    {
        # R implementation
        if (log) weights = exp(weights - max(weights)) # max(weights) for numerical stability
        return(weights/sum(weights))
    },
    {
        # C implementation
        n = length(weights)
        out = .C("renormalize_R", 
                 as.integer(n),
                 as.integer(log),  
                 weights=as.double(weights))
        return(out$weights)
    })   
}

