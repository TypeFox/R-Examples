#' Change the form of ordered sequences.
#' 
#' Simplifies sequences so that each tie group is only of distance 1 to the
#' next tie group.  For example, we would simplify (1, 1, 2, 4, 4, 5) to (1, 1,
#' 2, 3, 3, 4).
#' 
#' @param loss.time Matrix of sequences to be simplified.
#' @return Simplified sequences, as described in Description.
#' @author Erik Gregory
#' @keywords simplify sequence
SimplifySequences <-
function(loss.time) {
  maxs <- apply(loss.time, 1, max) + 1
  for (i in 1:nrow(loss.time)){
    nums <- sort(unique(loss.time[i, ]))
    for(j in 1:length(nums)) {
      loss.time[i, ][loss.time[i, ] == nums[j]] <- j
    }
  }
  return(loss.time)
}
