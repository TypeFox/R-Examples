
ent.weights = function(weights, engine="R") 
{
  weights = check.weights(weights, log=F, normalized=T)
  engine=pmatch(engine, c("R","C"))

  switch(engine,
  {
    # R implementation
    return(-sum(weights * log2(weights + .Machine$double.eps)))

    # Could take the maximum of this number and 0 to avoid negative results.
    # Others define this with log (ln) rather than log2.

  },
  {
    # C implementation
    out = .C("entropy_R", 
             as.integer(length(weights)),
             as.double(weights), 
             entropy = double(1))
    return(out$entropy)
  })
}

