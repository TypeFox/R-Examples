alleleEnrichment.perm <- function(
    subsample.size, 
    select.fun, 
    alleles, 
    var.min, 
    num.perms
  ) {
  
  cat(paste("Determining the best selection out of", num.perms, "permutations...\n"))
  
  trials <- sapply(
    1:num.perms, 
    function(x) {
      # set a zero vector to 1 (= selected) at num.samples.select random positions
      sel <- rep(0, ncol(alleles))
      sel[sample(1:ncol(alleles), subsample.size)] <- 1
      return(sel)
    }
  )
  
  # observed allele counts for trials (result vector): rows = SNPs and cols = selection trial
  obs.varsum <- alleles %*% trials
  
  # quality for balanced distribution: minimize the deviation from target frequencies
  qual <- apply(
    obs.varsum, 
    2, 
    function(obs.varsum.single) return(select.fun(var.min, obs.varsum.single, 0))
  )
  
  return(trials[, which.min(qual)])
  
}

### retired code: selects mode: max solution

# we check which trials are valid
# that is, they have at least the frequency of the risk allele as in the original population
# meet minimum allele count, <= (equality) is mandatory for zero alleles!
#    trials.accept <- apply(obs.varsum >= var.min, 2, all, na.rm = TRUE) 
#    if(!any(trials.accept))
#      stop(paste("Did not find a solution using", num.perms, "random trials with minimum allele counts of", paste(var.min, collapse = " / "), "\n", collapse = " "))
#    cat(paste("Found", length(which(trials.accept)), "solutions satisfying the minimum allele frequency constraints\n"))
#
#    qual <- apply(
#      obs.varsum, 
#      2, 
#      function(obs.varsum.single) return(deviation.weighted.maximize(var.min, obs.varsum.single, 0))
#    )
#    
#    # choose solution with best quality score and fulfilling min req
#    sol.idx <- which(qual == max(qual[which(trials.accept)]))[1]
#    return(trials[, sol.idx])
