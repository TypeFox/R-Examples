R2JAGS <-
function(model, data, nchains, inits, burnin, nodes, update, verbose){
  ## Disable JAGS progress bars
  old.pb <- options("jags.pb")
  on.exit(options(old.pb)) 
  options(jags.pb = "none")

  ## Format model
  wrap(model)
  
  ## Define inits
  if (any(inits == "random")) inits <- NULL

  ## Define & Initialize model
  mod <-
    jags.model(file = "modelTempFile.txt",
               data = data,
               #inits = inits,
               n.chains = nchains,
               n.adapt = 1000,
               quiet = !verbose)

  ## Delete 'modelTempFile.txt'
  unlink("modelTempFile.txt")

  ## Burn-in
  update(mod, n.iter = burnin, progress.bar = "none")

  ## Samples
  samples <- coda.samples(mod, nodes, n.iter = update, thin = 1,
                          progress.bar = "none")

  ## Deviance
  dic <- dic.samples(mod, n.iter = update, thin = 1, type = "pD",
                     progress.bar = "none")

  ## Return results
  return(list(mcmc.list = samples, dic = dic))
}