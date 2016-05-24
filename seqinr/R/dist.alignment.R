#
# Pairwise Distances from Aligned Protein or DNA/RNA Sequences
#

dist.alignment <- function(x, matrix = c("similarity", "identity"),gap )
{
  #
  # Check arguments:
  #
  if (!inherits(x, "alignment")) stop("Object of class 'alignment' expected")
  #
  # Match arguments:
  #
  matrix <- match.arg(matrix)
  #
  # Compute arguments for the C distance function:
  #
  sequences <- toupper(x$seq)
  nbseq <- x$nb
  matNumber <-ifelse(matrix == "similarity", 1, 2)
    #
    # The following shouldn't be hard encoded, an argument for full
    # user control should be added.
    #
  seqtype <- as.numeric(.Call("is_a_protein_seq", sequences[1], PACKAGE = "seqinr") >= 0.8)
  #
  # Call the C distance function:
  #
  if (missing(gap)) {
  	dist <- .Call("distance", sequences, nbseq, matNumber, seqtype,0, PACKAGE = "seqinr")
  	}
  else {
  	dist <- .Call("distance", sequences, nbseq, matNumber, seqtype,gap, PACKAGE = "seqinr")
  	}
  #
  # Convert the result in a object of class dist:
  #
  mat <- matrix(dist, nbseq, nbseq, byrow = TRUE)
  dimnames(mat) <- list(x$nam, x$nam)
  return( as.dist(mat) )
}
