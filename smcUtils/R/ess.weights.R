
ess.weights = function(weights, engine="R") 
{
  weights = check.weights(weights, log=F, normalized=T)
  engine=pmatch(engine, c("R","C"))

  switch(engine,
  {
    # R implementation
    return(1/sum(weights^2))
  
    # Could compute this as 
    # return(length(weights)/(1+cov.weights(weights)))
    # as in Liu and Chen 1995.

  },
  {
    # C implementation
    out = .C("ess_R", 
             as.integer(length(weights)),
             as.double(weights), 
             ess = double(1))
    return(out$ess)
  })
}

