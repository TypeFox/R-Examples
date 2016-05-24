
inverse.cdf.weights = function(weights, uniforms=runif(length(weights)), engine="R")
{
    check.weights(weights)
    check.weights(uniforms)
    
    engine=pmatch(engine, c("R","C"))
    if (is.unsorted(uniforms)) uniforms = sort(uniforms)
    n.samples = length(uniforms)

    switch(engine,
    {
        # R implementation
        ids       = integer(n.samples)
        cusum     = cumsum(weights)
        index     = 1
        for (i in 1:n.samples) 
        {
            found = FALSE
            while (!found) 
            {
                if (uniforms[i] > cusum[index]) 
                {
                    index = index + 1
                }
                else 
                {
                    found = TRUE
                }
            }
            ids[i] = index
        }
        return(ids)
    },
    {
        # C implementation
        out = .C("inverse_cdf_weights_R", 
                 as.integer(length(weights)),
                 as.double(weights),
                 as.integer(n.samples),
                 as.double(uniforms), 
                 id=integer(n.samples))
        return(out$id)
    })  
}



