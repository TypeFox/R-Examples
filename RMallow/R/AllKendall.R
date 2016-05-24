#' All Kendall's distances between two sets of rankings.
#' 
#' Calculates all of the Kendall's distances between two different sets of
#' rankings.
#' 

#' 
#' @param r One set of sequences.
#' @param seqs Another set of sequences.
#' @param data.info Optional argument, a 0/1/NA matrix specifying all of the
#' relevant information to calculate Kendall's difference for "r".  Used for
#' efficiency in "Solve".
#' @return Matrix where output[i, j] represents the distance from sequence "i"
#' in "r" to sequence "j" in "seqs".
#' @author Erik Gregory
#' @keywords Kendall distance
#' @examples
#' 
#' data1 <- do.call("rbind", list(1:5, 5:1, c(3, 2, 1, 4, 5)))
#' data2 <- do.call("rbind", list(1:5, 5:1))
#' # AllKendall(data1, data2)
#' 

AllKendall <-
function(r, seqs, data.info = NULL) {
  N <- nrow(r)
  n.seq <- nrow(seqs)
  dists <- matrix(0, nrow = N, ncol = n.seq)
  inds <- combn(ncol(r), 2)
  if (is.null(data.info)) {
    # Info for the data
    data.info <- KendallInfo(r, inds)
  }
  # Info for the sequences
  seqs.info <- KendallInfo(seqs, inds)
  if (n.seq > 1) {
    for (i in 1:nrow(seqs.info)) {
      dists[, i] <- rowSums(abs(sweep(data.info, 2, seqs.info[i, ], "-")), na.rm = TRUE)
    }
  }
  else {
    dists[, i] <- rowSums(abs(sweep(data.info, 2, seqs.info, "-")), na.rm = TRUE)
  }
  return(dists)
}
