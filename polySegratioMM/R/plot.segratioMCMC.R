`plot.segratioMCMC` <-
function(x, ..., row.index=c(1:10),
                                  var.index=c(1:6),
                                  marker.index=c(1:8))
{

  if (class(x) != "segratioMCMC")
    stop("'x' must be of class 'segratioMCMC'")
  
  ## based on plot.mcmc - maybe should look at usemethod
  
  T.index <- grep("T\\[",coda::varnames(x$mcmc.list))         # markers
  T.print <- T.index[marker.index]                      # marker selection
  drop <- T.index
  b.index <- grep("b\\[",coda::varnames(x$mcmc.list)) # random effects if set
  if (length(b.index>0)){
    b.print <- b.index[marker.index]
    drop <- c(T.index, b.index)
  }
  var <- c(1:dim(x$mcmc.list[[1]])[2])[-drop]      # variable selection
  var.index <- var[var.index]

  plot(x$mcmc.list [,c(var.index,T.print)], ... )

}

