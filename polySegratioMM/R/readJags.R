`readJags` <-
function(run.jags, quiet=TRUE, ...)
{
  ## wrapper to read.openbugs so ... can be used for start, end, thin
  ## returns object of class mcmc.list
  
  ## require(coda) # obsolete in recent versions of R as in Depends

  if (class(run.jags) != "runJags")
    stop("'run.jags' must be of class 'runJags'")

  ## superseded now that JAGS Version 1.0 required
  ##if (.Platform$OS.type == "windows"){# to be fixed after JAGS 0.90 superseded
  ##  x.jags <- read.jags(file=paste(run.jags$jags.control$stem,".out",sep=""),
  ##                      quiet=quiet, ...)
  ##  x.jags <- mcmc.list(x.jags)
  ##} else {
  x.jags <- coda::read.openbugs(stem = run.jags$jags.control$stem, quiet=quiet, ...)
  ##}
  
  res <- list(run.jags=run.jags, mcmc.list=x.jags)

  oldClass(res) <- "segratioMCMC"
  return(res)

}

