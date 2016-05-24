allele2to1 <-
  function (geno, marker.label = NULL, miss.val = NA){
    geno <- as.matrix(geno)
    if (!all(is.na(miss.val))) { 
      geno [geno %in% miss.val] <- NA
    }
    ns <- ncol(geno)
    if ((ns %% 2) != 0) { 
      stop("geno must have at least 2 loci.") 
    }
    ns <- ns / 2
    N <-  nrow(geno)
    code <- matrix (0, N, ns)
    for (i in 1:ns)
    {
      g1 <-  geno [, 2 * i - 1]
      g2 <-  geno [, 2 * i]
      tgeno <- sort(table (c(g1, g2)), decreasing = TRUE)
      tgeno <- tgeno[names(tgeno)!=0]
      if ((length(tgeno) != 2)) {
        stop (paste("Marker ", 
                    i, 
                    " is not biallelic.",
                    sep = "")) 
      }
      allele1 <- names(tgeno)[1]
      allele2 <- names(tgeno)[2]
      code[, i] <- rep(NA, N)
      code[, i] <- ifelse (((g1 == g2) & (as.character(g1) == allele1)),
                           0, 
                           ifelse (((g1==g2) & (as.character(g1)==allele2)),
                                   2, 
                                   ifelse ((g1 != g2), 1, code[, i])))
    }
    if (is.null(marker.label)) 
      marker.label <- paste("S", 1:ns, sep = "")
    colnames(code) <- marker.label
    rownames (code) <- rownames(geno)
    return (code)
  }
