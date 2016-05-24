hapest.caco <-
  function(geno, trait, lim = 0.05) {
    ca <- trait == 1
    co <- trait == 0
    nca <- sum(ca)
    nco <- sum(co)
    ck <- 2
    # infer haplotypes in pooled sample(cases and controls together)
    hpool <-  itegeppXXR(geno, des = 0, lim = lim)
    # failed inferring haplotypes ?
    if (all(is.na(hpool$hap))) {
      return(list(haplotypes = NA, desres = NA))
    }
    hapcaco.pool <- colSums(hpool$desres) /(ck * dim(hpool$desres)[1])
    hapcaco.ca <- colSums(hpool$desres [ca, ]) / (ck * nca)
    hapcaco.co <- colSums(hpool$desres [co, ]) / (ck * nco)
    hapcaco.pool <- hapcaco.pool[names(hapcaco.pool) != "R", drop = FALSE]
    hapcaco.ca   <- hapcaco.ca  [names(hapcaco.ca)   != "R", drop = FALSE]
    hapcaco.co   <- hapcaco.co  [names(hapcaco.co)   != "R", drop = FALSE]
    hapcaco <- data.frame(
      as.character(names(hapcaco.pool)), 
      as.numeric(hapcaco.pool), 
      as.numeric(hapcaco.ca), 
      as.numeric(hapcaco.co), 
      stringsAsFactors = FALSE)
    desres <- hpool$desres
    hapcaco[, 1] <- as.character(hapcaco[, 1])
    hapcaco[, 2] <- as.numeric(as.character(hapcaco[, 2]))
    hapcaco[, 3] <- as.numeric(as.character(hapcaco[, 3]))
    hapcaco[, 4] <- as.numeric(as.character(hapcaco[, 4]))
    colnames(hapcaco) <- c("Hap", "Pool", "Case", "Control")
    # Threshold for protecting haplotypes
    bin <-(hapcaco$Pool >= lim)
    hapcaco <- hapcaco[bin, , drop = FALSE]
    if( dim(hapcaco)[1] == 0) {
      print(paste("All inferred haplotypes or haplotype pairs less than ", 
                  "lim = ", lim, sep=""))
      return(list(haplotypes = NA, desres = NA))
    }
    hapsi <- which(!is.na(match(colnames(desres), hapcaco[, 1])))    
    desres <- as.matrix(desres[, hapsi, drop = FALSE])
    #   desres <- desres / 2
    if(!any(colnames(desres) == "R")) {
      desres <- cbind(desres, R = 2 - rowSums(desres))
    }
    desres <- as.matrix(desres)
    # round to 6 decimal places   
    hapcaco[, 2] <-  round(hapcaco[, 2], 6)
    hapcaco[, 3] <-  round(hapcaco[, 3], 6)
    hapcaco[, 4] <-  round(hapcaco[, 4], 6)
    hi <- match(hapcaco[, 1], colnames(desres))
    hapcaco <- hapcaco [!is.na(hi), , drop = FALSE]
    rownames(hapcaco) <- NULL
    return(list(haplotypes = hapcaco,
                desres = desres))
  }
