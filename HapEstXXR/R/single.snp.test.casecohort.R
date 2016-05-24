single.snp.test.casecohort <-
  function(snps, trait, 
           patid, start.time, stop.time, subcohort,
           stratvar = NA, 
           ties = "efron", 
           robust = FALSE, adj.var = NULL, prt = TRUE) {
    
    snps <- as.matrix(snps)
    if (any (!((trait == 0)|(trait == 1)))) {
      stop (paste("Error in single.snp.test.casecohort: ",
                  "case must be 1 (= case) or 0 (= non-case).",
                  sep = ""))
    }
    
    #if (method == "prentice" ) {
    res <- single.snp.test.casecohort.prentice (
      snps = snps,
      trait = trait, 
      patid = patid, 
      start.time = start.time, 
      stop.time = stop.time, 
      subcohort = subcohort,
      stratvar = stratvar, 
      robust = robust, 
      adj.var = adj.var,
      prt = prt,
      ties = ties)
    #}
    # Barlow is missing ............
    return (res)
  }
