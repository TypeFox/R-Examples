Preprocess.Tree <-
function(Tree)
{
  # The tree may be NULL or already processed
  if (is.null(Tree) | inherits(Tree, "PPtree")) return (Tree)

  # Tree must be either a phylog, phylo or a hclust object
  if (inherits(Tree, "phylog")) {
    # build an hclust object to use cutree later. Distances in $Wdist are actually 2*sqrt(distance)
    hTree <- stats::hclust(Tree$Wdist^2/2, "average")
    # build a phylo object
    phyTree <- ape::as.phylo.hclust(hTree)
    # edge lengths are divided by 2 during the conversion. See ?as.phylo.hclust
    phyTree$edge.length <- 2*phyTree$edge.length
  } else {
    if (inherits(Tree, "phylo")) {
      phyTree <- Tree
      # build an hclust object to use cutree later.
      hTree <- ape::as.hclust.phylo(Tree)
    } else {
      if (inherits(Tree, "hclust")) {
        hTree <- Tree
        # build a phylo object to use $droot later
        phyTree <- ape::as.phylo.hclust(Tree)
        # edge lengths are divided by 2 during the conversion. See ?as.phylo.hclust
        phyTree$edge.length <- 2*phyTree$edge.length
      } else {
        stop("Tree must be an object of class phylo, phylog or hclust")
      }
    }
  }

  # Calculate distances between nodes and leaves
  DistancesFromLeaves <- ape::branching.times(phyTree)
  # Get a sorted list of cuts (eliminate leaves)
  Cuts <- sort(DistancesFromLeaves[setdiff(names(DistancesFromLeaves), names(phyTree$leaves))])
  # Calculate intervals between cuts (add 0 to Cuts to get the first interval)
  Intervals <- diff(c(0, Cuts))
  # Eliminate 0 intervals (when a node contains more than 2 tips), including rounding errors
  RoundingError <- max(DistancesFromLeaves)*10*.Machine$double.eps
  Cuts <- Cuts[Intervals > RoundingError]
  Intervals <- Intervals[Intervals > RoundingError]

  ppTree <- list(
    phyTree   = phyTree,
    hTree     = hTree,
    Height    = Cuts[length(Cuts)],
    Cuts      = Cuts,
    Intervals = Intervals
    )
  class(ppTree) <- "PPtree"
  return (ppTree)
}
