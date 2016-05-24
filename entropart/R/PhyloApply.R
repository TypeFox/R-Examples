PhyloApply <-
function(Tree, FUN, NorP, Normalize = TRUE, ..., CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree: calculate cuts, intervals, and get it in both hclust and phylog formats.
  if(is.null(Tree)){
    stop("Tree cannot be NULL in PhyloApply")
  } else {
    ppTree <- Preprocess.Tree(Tree)
  }
  # The tree must be ultrametric
  if (!ape::is.ultrametric(ppTree$phyTree))
    stop("The tree must be ultrametric to apply a function over it.")

  # NorP may be a vector or a matrix. If it is a vector, double it into a matrix to simplify the code (use rownames)
  if (is.vector(NorP) | is.SpeciesDistribution(NorP)) {
    NorPisVector <- TRUE
  } else {
    if (length(dim(NorP)) == 1) {
      NorPisVector <- TRUE
    } else {
      NorPisVector <- (dim(NorP)[2] == 1)
    } 
  }
  # Save the name of NorP before manipulating it
  NorPName <- deparse(substitute(NorP))
  if (NorPisVector) {
    NorP <- matrix(NorP, nrow = length(NorP), ncol = 2, dimnames = list(names(NorP), c("NorP", "Dummy")))
  }
  # NorP should be named. If it is not, but has the same number of elements as the tree, just warn.
  if (is.null(rownames(NorP))) {
    if (nrow(NorP) == length(ppTree$phyTree$tip.label)) {
      rownames(NorP) <- ppTree$phyTree$tip.label
      warning("The abundance or frequency vector was not named. It was supposed to be order as the tree leaves.")
    } else {
      stop("The abundance or frequency vector is not named and does not have the same number of elements as the tree. Abundances and species could not be related.")
    }
  }

  # Eliminate abundances which are not in the tree (with a warning)
  SpeciesNotFound <- setdiff(rownames(NorP), ppTree$phyTree$tip.label)
  if (length(SpeciesNotFound) > 0) {
    NorP <- NorP[intersect(rownames(NorP), ppTree$phyTree$tip.label), ]
    if (nrow(NorP) > 1) { 
      # Some species have been dropped
      warning(paste("Species not found in the tree: ", SpeciesNotFound, collapse="; "))
    } else {
      # Less than 2 species were kept. Cannot calculate diversity.
      stop("Species cannot be found in the tree")
    }
  }
  
  # Rounding errors in cutree: is.unsorted(hTree$height) may return TRUE even though height is sorted by construction
  # Values are not sorted properly because of rounding errors, e.g. 4e-16 (2 * .Machine$double.eps) in a taxonomy where Cuts contains (1,2,3)
  # Run sort so that is.unsorted returns FALSE.
  OriginalHeights <- ppTree$hTree$height
  ppTree$hTree$height <- sort(OriginalHeights)
  # If there is no rounding error, add one (10 * .Machine$double.eps times the tree height) or cutree will miss some nodes.
  RoundingError <- max(ppTree$hTree$height) * 10 * .Machine$double.eps
  # Cut the tree at each node (eliminate the root). Cut at the values + RoundingError or many nodes will be missed.
  DatedGroups <- stats::cutree(ppTree$hTree, h=c(0, ppTree$Cuts[-length(ppTree$Cuts)]) + RoundingError)
  # DatedGroups is a table. Lines are species, columns are intervals between nodes.
  # Reorder NorP to fit DatedGroups
  NorP <- NorP[intersect(rownames(NorP), dimnames(DatedGroups)[[1]]),]
  # In each column, use the tip number as a factor to sum abundances of both NorP columns
  DatedN <- sapply(colnames(DatedGroups), function(group) apply(NorP, 2, function(n) tapply(n, as.factor(DatedGroups[rownames(NorP),group]), sum)), simplify=FALSE)
  # DatedN is a list of two-column matrices. It must be cleaned up if NorP was a vector
  if (NorPisVector) {
    DatedN <- lapply(DatedN, function(m) m[,1])
  }
  # Apply Function to each slice.
  DatedResult <- unlist(parallel::mclapply(DatedN, FUN, ..., CheckArguments = FALSE))
  # Names of slices should be the cut time, without the rounding error
  names(DatedResult) <- ppTree$Cuts
  # Normalization
  if (Normalize) {
    Normalization <- sum(ppTree$Intervals)
  } else {
    Normalization <- 1
  }
  # Return a weighted average of the results in each slice
  Value <- list(
    Distribution = NorPName,
    Function = deparse(substitute(FUN)),
    Tree = deparse(substitute(Tree)),
    Normalized = Normalize,
    Cuts = DatedResult, 
    Total = sum(DatedResult * ppTree$Intervals / Normalization)
  )
  class(Value) <- "PhyloValue"
  
  return(Value)  
  
}
