
cov.weights = function(weights, engine="R") 
{
  weights = check.weights(weights, log=F, normalized=T)
  engine=pmatch(engine, c("R","C"))

  switch(engine,
  {
    # R implementation
    return(var(weights)/mean(weights)^2)
  
    # Could compute this as
    # return(mean((length(weights)*weights-1)^2))
    # as in Liu and Chen 1995.

  },
  {
    # C implementation
    out = .C("cov2_R", 
             as.integer(length(weights)),
             as.double(weights), 
             cov2 = double(1))
    return(out$cov2) 
  })  
}

