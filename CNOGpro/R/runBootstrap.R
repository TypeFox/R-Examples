runBootstrap <-
function(experiment, replicates=1000, quantiles=c(0.025,0.975)){
  cat("Starting bootstrap analysis for", experiment$accession,"- This might take a while...\n")
  newexp <- copynumberBootstrap(experiment=experiment,replicates=replicates, quantiles=quantiles,windowlength=experiment$windowlength)
  cat("Finished bootstrap analysis.\n")
  experiment <- newexp
  return(experiment)
}
