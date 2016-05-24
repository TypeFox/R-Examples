#library(FrF2)

# DONE: Ensure that there is 1 or more entries in the list; pick the first one
# DONE: what happens if a full factorial is requested
# DONE: what happens if a "twice" or "four times" is requested

# TODO: add the defining relationship as an output. There are 2^p words in that
#       relationship.
# TODO: allow the user to specify the factor names in a list, and use these instead.
tradeoff <- function(runs=8, factors=7, display=TRUE){
  
  if(as.integer(runs) != runs){
    stop('The "runs" input must be an integer.')
  }
  if(as.integer(factors) != factors){
    stop('The "factors" input must be an integer.')
  }
  
  # Use FrF2 catalog to locate the experiments
  cat.entry <- FrF2::catlg
  subset_runs <- FrF2::nruns(cat.entry) %in% runs
  subset_factors <- FrF2::nfac(cat.entry) %in% factors
  if ((sum(subset_runs * subset_factors)) == 0){
    stop("The number of runs and/or number of factors do not specify a fractional factorial.\n 
         Please type: tradeOffTable()\n         for assistance.")
  }else{
    joint <- cat.entry[as.logical(subset_runs * subset_factors)]  
  }
  joint <- joint[[1]]
  resolution <- as.character(as.roman(joint$res))
  
  if (display){
    pr.header <- paste0("With ", runs, " experiments, and ", factors, " factors:\n")
    pr.res <- paste0("  Resolution: ", resolution, "\n")
    cat(pr.header)
    cat(pr.res)
  }
  # Fractional factorial experiments have 2^{k-p} number of runs, 
  # where factors == k and ngen == p.
  ngen  <- factors - log2(runs)
  frf2.catalog <- FrF2(nruns=runs, nfactors=factors)
  aliasing <- attr(frf2.catalog, "design.info")
  gen <- generators(frf2.catalog)
  if (display){
    if (ngen > 1){
      pr.gen <- paste0("  Generators:\n")
      for (elem in gen$generators){
        pr.gen <- paste0(pr.gen, "      ", elem, "\n")
      }
    }else{
      pr.gen <- paste0("  Generator: ", gen$generators, "\n")
    }
    pr.alias <- paste0("  Aliasing (related ONLY to main effects and 2-factor interactions):\n")
    if (resolution != "III"){
      #length(aliasing$aliased$main)==0) is an alternative test
      pr.alias <- paste0(pr.alias, "      Main effects are not aliased with 2-factor interactions.\n")
    }else{
      #pr.alias <- paste0(pr.alias, "    Main effect aliasing is:\n")
      for (elem in aliasing$aliased$main){
        pr.alias <- paste0(pr.alias, "      ", elem, "\n")
      }
    }
    
    if (length(aliasing$aliased$fi2) > 0){
      #pr.alias <- paste0(pr.alias, "    Two factor aliasing is:\n")
      for (elem in aliasing$aliased$fi2){
        pr.alias <- paste0(pr.alias, "      ", elem, "\n")
      }
    }
   
    pr.alias <- paste0(pr.alias, "\n")
    
    # Show the generator section and the aliasing section.
    cat(pr.gen)
    cat(pr.alias)
  }
  invisible(list(resolution=resolution, 
              generators=gen$generators, 
              aliases=aliasing$aliased))
}
if(FALSE){
  tradeoff(runs=8, factors=4)
  tradeoff(runs=8, factors=5)
  tradeoff(runs=8, factors=6)
  tradeoff(runs=8, factors=7)
  tradeoff(runs=16, factors=5)
  tradeoff(runs=16, factors=6)
  tradeoff(runs=16, factors=7)
  tradeoff(runs=16, factors=8)
  tradeoff(runs=16, factors=9)
}
