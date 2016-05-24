
resample = function(weights, num.samples=length(weights), 
                    method = c("stratified","residual","multinomial","systematic","branching"),
                    nonuniformity = c("none","ess","cov","entropy"), threshold=0.5,
                    rrf = "stratified", engine="R", log=TRUE, normalized=FALSE)
{
  method        <- match.arg(method)
  nonuniformity <- match.arg(nonuniformity)

  if (!normalized) weights = check.weights(weights, log=log, normalized=normalized)

  do.resample = FALSE
  switch(nonuniformity,
    "none"     = { do.resample = TRUE; },
    "ess"      = { if (ess.weights(weights)/num.samples      < threshold) do.resample = TRUE },
    "cov"      = { if (cov.weights(weights)/num.samples      > threshold) do.resample = TRUE },
    "entropy"  = { if (ent.weights(weights)/log2(num.samples)< threshold) do.resample = TRUE }
  )

  if (do.resample) {
    switch(method,
      "multinomial" = { ids = multinomial.resample(weights, num.samples, engine) },
      "residual"    = { ids = residual.resample(   weights, num.samples, engine, rrf) },
      "stratified"  = { ids = stratified.resample( weights, num.samples, engine) },
      "systematic"  = { ids = systematic.resample( weights, num.samples, engine) },
      "branching"   = { ids = branching.resample(  weights, num.samples, engine) }
    )
    weights = rep(1/num.samples,num.samples)
  } else { 
    ids = 1:length(weights)
  }

  return(list(weights=weights,indices=ids))
}

