################################################################################
# Minimize comparison sizes
################################################################################
# from a set of comparisons and lengths, determine the best possible
# arrangement of maps
minimize_comps <- function(comparisons, xlims, lengths, prel_offsets,
                           fixed_gap_length=FALSE){
  # function to minimize. Calculates the mean of the absolute differences
  # between starts and ends for direct hits only
  mean_w_fixed_gaps <- function(first_offset, fixed_offsets, ...){
    mean_w_gaps(c(first_offset, fixed_offsets), ...)
  }
  mean_w_gaps <- function(offsets, offsets_ref, xlim, xlim_ref, comp, side_ref){
    if (nrow(comp) == 0) return(0)
    side_test <- if (side_ref == 2) 1 else 2
    # recalc for ref side
    comp <- calc_comp_coor(gap=offsets_ref, xlim=xlim_ref, comp=comp,
                           side=side_ref)
    # test
    comp <- calc_comp_coor(gap=offsets, xlim=xlim, comp=comp, side=side_test)
    direction <- sign(comp$end1-comp$start1)*sign(comp$end2-comp$start2)
    dists <- c(abs(comp$start1-comp$start2),
             abs(comp$end1-comp$end2))#[direction > 0]
    lengths <- abs(comp$end1-comp$start1) + abs(comp$end2-comp$start2)
    lengths[lengths==0] <- 1
    weighted.mean(dists, c(lengths, lengths))
  }
  
  n_org <- length(lengths)
  offsets <- prel_offsets
  if (length(comparisons) < 1) return(offsets)
  idx_ref <- which.max(lengths)
  max_len <- max(lengths)

  # go up from ref
  if (idx_ref > 1){
    # comp i is between org i and i+1
    for (i in (idx_ref-1):1){
      # optimise
      # if fixed_gap_length, optimse only on the first gap
      if (fixed_gap_length && length(offsets[[i]]) > 1){
        fixed_offsets <- offsets[[i]][2:length(offsets[[i]])]
        opt <- optim(par=offsets[[i]][1], fn=mean_w_fixed_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i]][1],
                     fixed_offsets=fixed_offsets, offsets_ref=offsets[[i+1]],
                     xlim=xlims[[i]], xlim_ref=xlims[[i+1]],
                     comp=comparisons[[i]], side_ref=2)
        offsets[[i]][1] <- opt$par
      }
      # else optimise on all offsets
      else {
        opt <- optim(par=offsets[[i]], fn=mean_w_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i]],
                     offsets_ref=offsets[[i+1]],
                     xlim=xlims[[i]], xlim_ref=xlims[[i+1]],
                     comp=comparisons[[i]], side_ref=2)
        offsets[[i]] <- opt$par
      }
    }
  }
  # go down
  if (idx_ref < n_org){
    for (i in idx_ref:(n_org-1)){
      # optimise
      # if fixed_gap_length, optimse only on the first gap
      if (fixed_gap_length && length(offsets[[i]]) > 1){
        fixed_offsets <- offsets[[i]][2:length(offsets[[i]])]
        opt <- optim(par=offsets[[i+1]][1], fn=mean_w_fixed_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i+1]][1],
                     fixed_offsets=fixed_offsets, offsets_ref=offsets[[i]],
                     xlim=xlims[[i+1]], xlim_ref=xlims[[i]],
                     comp=comparisons[[i]], side_ref=1)
        offsets[[i+1]][1] <- opt$par
      }
      # else optimise on all offsets
      else {
        opt <- optim(par=offsets[[i+1]], fn=mean_w_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i+1]],
                     offsets_ref=offsets[[i]],
                     xlim=xlims[[i+1]], xlim_ref=xlims[[i]],
                     comp=comparisons[[i]], side_ref=1)
        offsets[[i+1]] <- opt$par
      }
    }
  }
  offsets
}
