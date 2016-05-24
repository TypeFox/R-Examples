findKtree <-
function(bb, parent, sons, compBranch, n) {
  kappaTop <- numeric(bb)
  kappaBottom <- numeric(bb)
  for (i in seq_len(bb)) {
    if (is.na(parent[i])) {
      kappaBottom[i] <- NA
      kappaTop[i] <- NA
    } else if (parent[i] == 0) { 
        kappaBottom[i] <- 0 
        if (i <= length(sons)) {
          Ksons <- sons[[i]]
          if (!is.null(Ksons) && !is.na(Ksons)) {
          kappaTop[i] <-
              (length(compBranch[[i]]) - length(unlist(compBranch[Ksons]))) / n
          } else {
            kappaTop[i] <- length(compBranch[[i]]) / n
          }
        } else {
          kappaTop[i] <- length(compBranch[[i]]) / n
        }
      } else {
        kappaBottom[i] <- kappaTop[parent[i]]
        if (i <= length(sons)) {
          Ksons <- sons[[i]]
          if (!is.null(Ksons) && length(Ksons) != 0 && !is.na(Ksons)) {
            kappaTop[i] <- kappaBottom[i] + (length(compBranch[[i]]) -
                length(unlist(compBranch[Ksons]))) / n
          } else {
            kappaTop[i] <- kappaBottom[i] + length(compBranch[[i]]) / n
          }
        } else {
          kappaTop[i] <- kappaBottom[i] + length(compBranch[[i]]) / n
        }
      }
    }

  return(list("kappaTop" = kappaTop, "kappaBottom" = kappaBottom))
}
