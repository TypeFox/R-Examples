filter.cmap <- function(cm, cutoff.sims = dim(cm)[3]) {

  ## Check input
  if (length(dim(cm)) != 3) {
    stop("Input 'cm' should be a NxNxS 3d array,\n 
      where N is the number of atoms and S is the number of simulations")
  }

  if (!(is.numeric(cutoff.sims) && (cutoff.sims <= dim(cm)[3]))) {
    stop("Input 'cutoff.sims' should a number between 1 and dim(cm)[3],\n 
      i.e. the number of simulations upon which filtering is based")
  }

  if ((is.numeric(cutoff.sims)) <= 0) {
    print("WARNING!! cutoff.sims' should a number between 1 and dim(cm)[3],\n 
      i.e. the number of simulations upon which filtering is based")
  }

  ## Sum across simulations and filter by cutoff.sims de
  cm.sum <- apply(cm, c(1:2), sum, na.rm = TRUE)

  return(array(as.numeric(cm.sum >= cutoff.sims), dim = dim(cm.sum)))
}
