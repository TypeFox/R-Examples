alleleEnrichment.lp <- function(
    subsample.size, 
    mode, 
    alleles, 
    var.min, 
    weights = weights,
    search.seconds
  ) {
  
  
  require(lpSolveAPI)
  
  # iteratively increase var.min until desired number of samples is exploited
  var.min.incr <- var.min
  # to-be solution, a 0-1 vector where 1 is a selected sample
  sol <- NA
  
  inc.fun <- lapply(
               rownames(alleles), # over snps
               function(x) {
                 if(mode[[x]] == "max") {
                   return(function(old.min) return(old.min + (old.min / subsample.size) * weights[[x]]))
                 } else {
                   return(identity)
                 }
               }
             )
  
  while(TRUE) {

    lp <- make.lp(dim(alleles)[1], dim(alleles)[2])
    lapply(
      1:ncol(alleles), 
      function(idx) {
        set.column(lp, idx, alleles[, idx])
        set.type(lp, idx, "binary")           # the decision variable, i.e. select this sample or not, has to be binary
      }
    )
    set.objfn(lp, rep(1, ncol(alleles))) # objective function is the sum of samples
    set.constr.type(lp, rep(">=", nrow(alleles)), 1:nrow(alleles))
    set.rhs(lp, floor(var.min.incr)) # is much faster with int constraints, also allele sums cannot be real numbers
    
    # remark: setting bb.depthlimit does not seem to work properly, does not terminate in experimental runs
    lp.control(lp, timeout = search.seconds, sense = "min")
    
    #lp.control(lp, bb.rule = c("pseudoratio", "stronginit", "randomize"))
    
    cat(paste("Searching a solution for constraints", paste(floor(var.min.incr), collapse = "/"), "(timeout:", search.seconds, "secs)...\n"))
    status <- solve(lp)
    
    # status == 0 means optimal solution found, 1 suboptimal solution found
    if(status >= 2) {
      # when this is the first run, terminate with error
      if(all(is.na(sol)))
        stop("No solution found / premature abort! Try again with different parameters...\n")
      # otherwise, break the search, this keeps the last solution
      break
    }
    
    sol.new <- get.variables(lp)
    if(sum(sol.new) > subsample.size) {
      if(all(is.na(sol)))
        stop("No solution found / premature abort! Try again with different parameters...\n")
      break
    } else {
      sol <- sol.new
    }
    
    # next af to try
    #      af.target <- af.target + (1-af.target) * weight * 0.01 # approx. 1% increase of min af
    #      var.min.incr <- floor(af.target * 2*apply(alleles, 1, function(x) sum(!is.na(x))) * subsample.size / ncol(alleles))
    var.min.incr <- mapply(do.call, inc.fun, lapply(var.min.incr, as.list))
    
    # this generates strange messages so we do not tidy up... :-(
    # delete.lp(lp)
    
  }

  return(sol)
  
}