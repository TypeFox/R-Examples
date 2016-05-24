#' Update modal sequences in each cluster.
#' 
#' Maximizes the likelihood of the data by updating the cluster centers of the
#' model.
#' 
#' @param r Matrix of sequences being clustered.
#' @param z Probability of cluster membership for each sequence and each
#' cluster.
#' @param infos The KendallInfo matrix for "r".
#' @return New cluster centers for each cluster.
#' @author Erik Gregory
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords cluster center
UpdateR <-
function(r, z, infos = NULL) {
  # Number of modes
  G <- ncol(z)
  # Number of subjects.
  N <- nrow(r)
  if (is.null(infos)) {
    # Information in the sequences.
    infos <- KendallInfo(r)
  }
  # Number of comparisons for Kendall's distance.
  n <- ncol(infos)
  # Number of abilities
  n.abils <- ncol(r)
  zeros <- list()
  ones <- list()
  zero.sums <- matrix(0, nrow = G, ncol = n)
  one.sums <- zero.sums
  if (G > 1) {
    for (i in 1:n) {
      # Weighted values of orderings.
      zero.sums[, i] <- colSums((z[which(infos[, i] == 0), ]))
      one.sums[, i] <- colSums(z[which(infos[, i] == 1), ])
    }
  }
  else {
    for (i in 1:n) {
      zero.sums[, i] <- sum(z[which(infos[, i] == 0), ])
      one.sums[, i] <- sum(z[which(infos[, i] == 1), ])
    }
  }
  prefs <- zero.sums - one.sums
  # Make sure none of them are zero...
  # print(paste(length(which(prefs == 0)), "zero values..."))
  prefs[prefs > 0] <- 0
  prefs[prefs < 0] <- 1
  # Reverse-engineer the sequences.
  R <- ConstructSeqs(prefs, n.abils)
  return(R)
}
