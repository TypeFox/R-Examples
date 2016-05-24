single.haplotype.test <-
  function(snps, trait, famid, patid, fid, mid, 
           adj.var = NULL, type = c("gaussian", "binomial", "families"),
           prt = TRUE, lim = 0.05, min.count = 10, 
           alpha = 0.05, sort = FALSE) {
    
    type <- match.arg(type)
    snps <- as.matrix(snps)
    switch(type, 
           gaussian = test <- "F", 
           binomial = test <- "Chisq", 
           families = test <- "wTDT", 
           survival = { print("not ready!"); stop("end") })
    
    trait <- as.numeric(trait)
    if (type == "families") {
      famid <- as.character(famid)
      lenx  <- length(famid)
      patid <- as.character(patid)
      fid   <- as.character(fid)
      mid   <- as.character(mid)
      if ((length(patid) !=  lenx) |(length(fid) !=  lenx) |
            (length(mid) !=  lenx) |(length(trait) !=  lenx) |
            (dim(snps)[1] !=  lenx))
        stop(paste("length mismatch of input parameter.\n", sep = ""))
    }
    
    snps <- as.matrix(snps)
    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    if (length(trait) !=  N)
      stop("length of trait does not match dimension of snps")
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
      if (type == "families") {
        famid <- famid[-miss.value]
        patid <- patid[-miss.value]
        fid   <- fid  [-miss.value]
        mid   <- mid  [-miss.value]
        lenx<-length(famid)
      }
      trait <- as.numeric(trait[-miss.value])
      snps  <- snps[-miss.value, , drop = FALSE]
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
      cat("Start procedure: single.haplotype.test\n")
      if (type == "families") {
        cat("Families:              ", 
            length(unique(famid)), "\n", sep = "")
      }
      cat("Individuals:           ", N, "(", 
          N - dim(snps)[1], " individual(s) exluded)\n", sep = "")
      cat("SNPs:                  ", ns, "\n", sep = "")
      cat("Trait type:            ", type, 
          ifelse(adjusted, "(adjusted)", "(unadjusted)"), "\n", sep = "")
      cat("Statistic:             ", test, "\n", sep = "")
      cat("Threshold(lim):        ", lim, "\n", sep = "")
      if (!adjusted) {
        cat("Number of covariates:  0\n\n")
      } else {
        cat("Number of covariates:  ", dim(adj.var)[2], "\n\n", sep = "")
      }
    }
    ############################
    # single haplotype analysis
    switch(type, 
           gaussian = res <- single.haplotype.test.gaussian(
             snps, 
             trait, 
             adj.var = adj.var, 
             lim = lim, 
             baseline.hap = "max", 
             do.hap.specific.test = TRUE, 
             min.count = min.count, 
             alpha = alpha), 
           binomial = res <- single.haplotype.test.binomial(
             snps, 
             trait, 
             adj.var = adj.var, 
             lim = lim, 
             baseline.hap = "max", 
             do.hap.specific.test = TRUE, 
             min.count = min.count, 
             alpha = alpha), 
           families   = res <- single.haplotype.test.families(
             snps, 
             trait, 
             famid, 
             patid, 
             fid, 
             mid, 
             lim = lim), 
           survival = { print("not ready!") ; stop("end") })
    return(res)
  }
