`diagnosticsJagsMix` <-
function(mcmc.mixture, diagnostics=TRUE, plots=FALSE,
           index=-c( grep("T\\[",varnames(mcmc.mixture$mcmc.list)),
                     grep("b\\[",varnames(mcmc.mixture$mcmc.list)) ),
           trace.plots = FALSE, auto.corrs = FALSE, density.plots = FALSE,
           xy.plots = FALSE, hpd.intervals = FALSE, hdp.prob = 0.95,
           return.results = FALSE)
{
  ## Purpose: print and plot 'coda' diagnostics after fitting dosage
  ##          mixture model using JAGS. The usual parameters of
  ##          interest are selected automatically or may be specified

  ## Arguments:
  ## sumdata:       mcmc summary data produce by 'readJagsMix'
  ## diagnostics:   if TRUE then print several 'coda' dignostic tests
  ## plots:         if TRUE then produce several 'coda' dignostic plots
  ## index:         index of parameters for disgnostic tests/plots
  ##                default: mixture model (and random effects) parameters
  ## trace.plots:   if TRUE plot mcmc traces (default: FALSE)
  ## auto.corrs:    produce autocorrelations of mcmc's (default: FALSE)
  ## density.plots: if TRUE plot parameter densities  (default: FALSE)
  ## xy.plots:      if TRUE plot traces using 'lattice' (default: FALSE)
  ## hpd.intervals: if TRUE print and return highest posterior density
  ##                intervals for parameters specified by 'index'
  ## hdp.prob:          probability for 'hpd.intervals'
  ## return.results:  return results as list

  ## Values:
  ## HPDinterval:   if 'hpd.intervals' return highest posterior density
  ##                intervals for parameters specified by 'index'
  ## autocorr:      if auto.corrs is TRUE return autocorrlations of parameters
  ##                specified by 'index'

  ## require(coda) # obsolete for current version of R since in Depends
  
  if ( !( class(mcmc.mixture) == "segratioMCMC" |
       class(mcmc.mixture) == "runJagsWrapper")) 
    stop("'mcmc.mixture' must be of class 'segratioMCMC' or 'runJagsWrapper'")

  if ( class(mcmc.mixture) == "runJagsWrapper")
    mcmc.mixture <- mcmc.mixture$mcmc.mixture

  res <- list()
  res$raftery <- res$geweke <- res$heidel <- NULL
  
  ##cat("Info: Starting raftery.diag in diagnosticsJagsMix\n")
  options(show.error.messages = FALSE)
  try(res$raftery <- coda::raftery.diag(mcmc.mixture$mcmc.list[,index]), silent=TRUE)
  ##cat("Info: Finished raftery.diag in diagnosticsJagsMix\n")
  options(show.error.messages = TRUE)
  if (length(res$raftery)==0)
     res$raftery  <-
      "Error: Raftery diagnostics not computed due to convergence problems"

  ##cat("Info: Starting geweke.diag in diagnosticsJagsMix\n")
  options(show.error.messages = FALSE)
  try(res$geweke <- coda::geweke.diag(mcmc.mixture$mcmc.list[,index]), silent=TRUE)
  ##cat("Info: Finished geweke.diag in diagnosticsJagsMix\n")
  options(show.error.messages = TRUE)
  if (length(res$geweke)==0)
    res$geweke <-
      "Error: Geweke diagnostics not computed due to convergence problems"

  res$heidel <- NULL
  options(show.error.messages = FALSE)
  ##cat("Info: Starting heidel.diag in diagnosticsJagsMix\n")
  try(res$heidel <- coda::heidel.diag(mcmc.mixture$mcmc.list[,index]), silent=TRUE)
  ##cat("Info: Finished heidel.diag in diagnosticsJagsMix\n")
  options(show.error.messages = TRUE)
  if (length(res$heidel)==0)
    res$heidel <-
      "Error: Heidelberg and Welch's diagnostics not computed due to convergence problems"
  
  if(diagnostics){
    cat("\nRaftery and Lewis's diagnostic\n")
    print(res$raftery)
    
    cat("\nGeweke's convergence diagnostic\n")
    print(res$geweke)

    cat("\nHeidelberger and Welch's convergence diagnostic\n")
    print(res$heidel)
  }

  if(plots) {
    
    coda::acfplot(mcmc.mixture$mcmc.list[,index])

    for (chain in 1:length(mcmc.mixture$mcmc.list)){
      levelplot(mcmc.mixture$mcmc.list[[chain]][,index],
                main=paste("Chain:",chain))
    }

    qqmath(mcmc.mixture$mcmc.list[,index])

    coda::cumuplot(mcmc.mixture$mcmc.list[,index])

  }

  if (trace.plots) {
    coda::traceplot(mcmc.mixture$mcmc.list[,index])
  }

  if (xy.plots) {
    xyplot(mcmc.mixture$mcmc.list[,index])
  }
    
  if (density.plots) {
    densityplot(mcmc.mixture$mcmc.list[,index])
  }
  
   if(auto.corrs) {
    print(acorr <- coda::autocorr(mcmc.mixture$mcmc.list[,index]))
    if (return.results) res$autocorr <- acorr
  }
  

  res$hpd <- coda::HPDinterval(mcmc.mixture$mcmc.list[,index], prob=hdp.prob)
  if(hpd.intervals){
    print(res$hpd)
  }

  if(return.results) {
    return(res)
  } else {
    return(NULL)
  }

}

