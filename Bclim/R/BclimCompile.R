BclimCompile <-
function(Layers,Mixtures,MCMC,Interpolations,core.name="Core") {
 
  clim.dims <- c("GDD5","MTCO","AET/PET x 1000")
  results <- list(time.grid=Interpolations$time.grid,core.name=core.name,clim.interp=Interpolations$clim.interp,v.interp=Interpolations$v.interp,MDP=Layers,ScMean=Mixtures$ScMean,ScVar=Mixtures$ScVar,clim.dims=clim.dims,n=Mixtures$n,m=Mixtures$m,n.samp=Mixtures$n.samp,Chronsfile=MCMC$chron.loc,nchron=MCMC$nchron,chron.store=MCMC$chron.store)
  class(results) <- "Bclim"
  return(results)
}
