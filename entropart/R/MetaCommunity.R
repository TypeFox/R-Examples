MetaCommunity <-
function(Abundances, Weights = rep(1, ncol(Abundances)))
{
  Nspecies <- length(Abundances[, 1])
  if (is.factor(Abundances[,1])) {
    FirstColumnOfData <- 2
    SpeciesNames <- Abundances[,1]
    Ncommunities <- length(Abundances[1, ])-1
  } else {
    FirstColumnOfData <- 1
    if (is.null(rownames(Abundances))) {
      # Create species names
      SpeciesNames <- as.factor(paste("sp", 1:(nrow(Abundances)), sep=""))
    } else {
      # Read species names
      SpeciesNames <- as.factor(rownames(Abundances))
    }
    Ncommunities <- length(Abundances[1, ])
  }
  
  # Community names
  if (is.null(colnames(Abundances))) {
    # Create community names
    colnames(Abundances) <- paste("P", 1:(ncol(Abundances)), sep="")
  }
  
  # Matrix containing p_si
  Nsi <- as.matrix(Abundances[, FirstColumnOfData:length(Abundances[1, ])])
  dimnames(Nsi)[[1]] <- SpeciesNames
  # Vector of weights
  if (is.vector(Weights)) {
    Wi <- Weights/sum(Weights)
  } else {
    Wi <- Weights$Weights/sum(Weights$Weights)
  }
  
  # Name the weight vector
  names(Wi) <- colnames(Nsi)  
  
  Preprocess.MC(Nsi, Wi)
}


is.MetaCommunity <-
function (x) 
{
  inherits(x, "MetaCommunity")
}


plot.MetaCommunity <- 
function (x, ...) 
{
  graphics::barplot(cbind(x$Psi,
                rep(0, x$Nspecies),
                x$Ps              
  ),
  beside = FALSE,
  width = c(x$Wi, .5, 1),
  names.arg = c(names(x$Wi), "", "Metacommunity"),
  ylab = "Species frequencies",
  ...
  )
}


summary.MetaCommunity <-
function(object, ...) 
{
  
  cat("Meta-community (class 'MetaCommunity') made of", object$N, "individuals in", object$Ncommunities, "communities and", object$Nspecies, "species.", "\n", fill=TRUE)
  cat(paste("Its sample coverage is", object$SampleCoverage, "\n"), fill=TRUE)
  cat("Community weights are:", "\n")
  print(object$Wi)
  cat("Community sample numbers of individuals are:", "\n")
  print(object$Ni)
  cat("Community sample coverages are:", "\n")
  print(object$SampleCoverage.communities)
  
  return(invisible(NULL))
}