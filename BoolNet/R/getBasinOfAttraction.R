# Determine the basin of attraction of attractor <attractorNo> in <attractorInfo>.
getBasinOfAttraction <- function(attractorInfo,attractorNo)
{
  stopifnot(inherits(attractorInfo,"AttractorInfo") || inherits(attractorInfo,"SymbolicSimulation"))

  if (missing(attractorNo) || attractorNo <= 0 || attractorNo > length(attractorInfo$attractors))
    stop("Please provide a valid attractor number!")
  
  table <- getTransitionTable(attractorInfo)
  return(table[which(table$attractorAssignment == attractorNo),,drop=FALSE])
}
