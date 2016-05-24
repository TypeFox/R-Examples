`summary.segratioMCMC` <-
function(object, ..., row.index=c(1:10),
                                  var.index=NULL, marker.index=c(1:8))
{

  
  if (class(object) != "segratioMCMC")
    stop("'object' must be of class 'segratioMCMC'")

  ## print object summary of type segratioMCMC
  
  T.index <- grep("T\\[",coda::varnames(object$mcmc.list))         # markers
  T.print <- T.index[marker.index]                      # marker selection
  drop <- T.index
  b.index <- grep("b\\[",coda::varnames(object$mcmc.list)) # random effects if set
  if (length(b.index>0)){
    b.print <- T.index[marker.index]
    drop <- c(T.index, b.index)
  }
  var <- c(1:dim(object$mcmc.list[[1]])[2])[-drop]      # variable selection
  if (length( var.index) == 0) {  # take all vars if not set
    var.index <- var
  } else {
    var.index <- var[var.index]
  }

  options(show.error.messages = FALSE)
  ##cat("Info: Starting summary calculations in summary.segratioMCMC\n")
  try(res <- summary(object$mcmc.list[,c(var.index,T.print)]), silent=TRUE)
  ##print(.Last.value)
  ##cat("Info: Finished summary calculations in summary.segratioMCMC\n")
  options(show.error.messages = TRUE)

  class(res) <- "summarySegratioMCMC"
  return(res)    # res is of class "summary.mcmc"

}

