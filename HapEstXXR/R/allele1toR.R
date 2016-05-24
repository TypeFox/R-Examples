allele1toR <-
  function(geno, marker.label = NULL, miss.val = c(-1, NA)){
    geno <- as.matrix(geno)
    if (!all(is.na(miss.val))) {
      geno[geno %in% miss.val ] <- NA  
    }
    if (any(!(geno %in% c(0, 1, 2, miss.val)))) 
      stop ("Error in allele1toR: miss match in genotypes")
    geno[!(geno %in% c(0, 1, 2))] <- NA
    geno <- replace(geno, (geno %in% 1), 3)
    geno <- replace(geno, (geno %in% 0), 1)
    geno[is.na(geno) ] <- 0
    if (is.null(marker.label))
      marker.label <- paste("S", 1:ncol(geno), sep = "")
    colnames(geno) <- marker.label
    rownames (geno) <- rownames(geno)
    return ( geno)
  }
