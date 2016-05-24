hapest.gaussian <-
  function(geno, trait, lim = 0.05) {
    n <- length(trait)
    ck <- 2
    # infer haplotypes
    hpool <-  itegeppXXR(geno, des = 0, lim = lim)
    # failed inferring haplotypes ?
    if (all(is.na(hpool$hap))) {
      return(list(haplotypes = NA, desres = NA))
    }
    hapgauss.pool <- colSums(hpool$desres) /(ck * dim(hpool$desres)[1])
    hapgauss.pool <- hapgauss.pool[names(hapgauss.pool) != "R", drop = FALSE]
    hapgauss <- data.frame(
      as.character(names(hapgauss.pool)), 
      as.numeric(hapgauss.pool),  
      stringsAsFactors = FALSE)
    desres <- hpool$desres
    hapgauss[, 1] <- as.character(hapgauss[, 1])
    hapgauss[, 2] <- as.numeric(as.character(hapgauss[, 2]))
    colnames(hapgauss) <- c("Hap", "Pool")
    # Threshold for protecting haplotypes
    bin <-(hapgauss$Pool >= lim)
    hapgauss <- hapgauss[bin, , drop = FALSE]
    if( dim(hapgauss)[1] == 0) {
      print(paste("All inferred haplotypes or haplotype pairs less than ", 
                  "lim = ", lim, sep=""))
      return(list(haplotypes = NA, desres = NA))
    }
    hapsi <- which(!is.na(match(colnames(desres), hapgauss[, 1])))    
    desres <- as.matrix(desres[, hapsi, drop = FALSE])
    #   desres <- desres / 2
    if(!any(colnames(desres) == "R")) {
      desres <- cbind(desres, R = 2 - rowSums(desres))
    }
    desres <- as.matrix(desres)
    # round to 6 decimal places   
    hapgauss[, 2] <-  round(hapgauss[, 2], 6)
    hi <- match(hapgauss[, 1], colnames(desres))
    hapgauss <- hapgauss [!is.na(hi), , drop = FALSE]
    rownames(hapgauss) <- NULL
    return(list(haplotypes = hapgauss,
                desres = desres))
  }
