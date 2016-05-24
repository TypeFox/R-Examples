single.snp.test <-
  function(snps, trait, adj.var = NULL, 
           type = c("gaussian", "binomial", "families", "casecohort"), 
           famid, patid, fid, mid, 
           start.time, stop.time, subcohort, stratvar = NA, robust = FALSE, 
           marker.label = NA, 
           prt = TRUE, ties = "efron") {
    
    type <- match.arg(type)
    switch(type, 
           gaussian =  test <- "F", 
           binomial =  test <- "Chisq", 
           families  = test <- "wTDT", 
           casecohort = test <- "case.cohort"
    )
    snps <- as.matrix(snps)
    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    if (!all(is.na(marker.label))) {
      if (ns != length(marker.label)) {
        stop(paste("Error in single.snp.test: ", 
                   "marker.label don't matched nummber of SNPs.",
                   sep = ""))
      }
    }
    if (length(trait) !=  N) {
      stop("length of trait does not match dimension of snps")
    }
    adjusted <- FALSE
    if (!all(is.null(adj.var))) { 
      adjusted <- TRUE
    }
    # handle missing values
    miss.value <- which(is.na(trait))
    if (adjusted) {
      adj.var <- as.matrix(adj.var)
      miss.value <- unique(c(miss.value,
                             which(apply(is.na(adj.var), 1, any))))
    }
    if (length(miss.value) > 0) {
      if (adjusted) {
        adj.var <- adj.var[-miss.value, , drop = FALSE]
      }
      trait <- as.numeric(trait[-miss.value])
      snps  <- snps[-miss.value, , drop = F]
      N <- dim(snps)[1]
    }
    # Adjustment
    if (adjusted) {
      adj.var <- as.matrix(adj.var)
      if (dim(adj.var)[1] !=  N)
        stop(paste("Error in stepwise: ",
                   "length of patid does not match number of rows in adj.cov",
                   sep = ""))
    }
    ######################
    # print    
    if (prt == TRUE) {
      cat("Start procedure: single.snp.test\n")
      cat("Individuals:          ", N, "(", length(miss.value), 
          " excluded)", "\n", sep = "")
      cat("SNPs:                 ", ns, "\n", sep = "")
      cat("Trait type:           ", type, "\n", sep = "")
      cat("Statistic:            ", test, "\n", sep = "")
      if (is.null(dim(adj.var))) {
        cat("Number of covariates: 0\n\n", sep = "")
      } else {
        cat("Number of covariates: ", dim(adj.var)[2], "\n\n", sep = "")
      }
    }
    ######################
    # single SNP analysis
    switch(type, 
           gaussian   = res <- single.snp.test.gaussian(snps,
                                                        trait,
                                                        adj.var = adj.var, 
                                                        prt = prt), 
           binomial   = res <- single.snp.test.binomial(snps, 
                                                        trait,
                                                        adj.var = adj.var,
                                                        prt = prt), 
           families   = res <- single.snp.test.families(snps, 
                                                        trait,
                                                        adj.var = adj.var, 
                                                        famid, 
                                                        patid, 
                                                        fid, 
                                                        mid, 
                                                        prt = prt), 
           casecohort = res <- single.snp.test.casecohort(
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
             ties = "efron"))  
    if (!all(is.na(marker.label))) {
      res[, "SNP"] <- marker.label
    }
    return(res)    
  }
