`print.segratioMCMC` <-
function(x, ..., row.index=c(1:10),
                                  var.index=c(1:6),
                                  marker.index=c(1:8), chain=1)
{

  if (class(x) != "segratioMCMC")
    stop("'x' must be of class 'segratioMCMC'")

  ## print object of type segratioMCMC
  
  cat("MCMC runs for mixture model with",
      x$run.jags$jags.control$model$n.components,
      "dosage classes\n  ploidy level",
      x$run.jags$jags.control$model$ploidy.level,
      "(",x$run.jags$jags.control$model$E.segRatio$ploidy.name,")\n")

  cat("MCMC Burn in:", x$run.jags$jags.control$burn.in, ", Sample:",
      x$run.jags$jags.control$sample,"\n")

  cat("JAGS File:",x$run.jags$jags.control$bugs.file,"\n")

  print(x$run.jags)
  
  T.index <- grep("T\\[",coda::varnames(x$mcmc.list))    # markers
  T.print <- T.index[marker.index]                 # marker selection
  drop <- T.index
  b.index <- grep("b\\[",coda::varnames(x$mcmc.list)) # random effects if set
  if (length(b.index>0)){
    b.print <- T.index[marker.index]
    drop <- c(T.index, b.index)
  }
  var <- c(1:dim(x$mcmc.list[[1]])[2])[-drop]   # variable selection
  var.index <- var[var.index]

  
  for (i in length(chain))
    {
      cat("Chain:",chain,"\nSubset:\n")
      print(x$mcmc.list[ row.index , var.index])
      print(x$mcmc.list[ row.index , T.print] )
      
    }
}

