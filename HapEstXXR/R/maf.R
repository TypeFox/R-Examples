maf  <- 
  function(geno, marker.label = NA) {
    geno <- as.matrix(geno)
    maf <- matrix(0, ncol = 9, nrow = dim(geno)[2])
    colnames(maf) <- c(0:3, "Total", "call.rate", "minor.allele",
                       "maf", "hwe.chisq.p.value")
    for(i in 1:dim(geno)[2]) {
      tabl <- table(geno[, i])
      maf [i, names(tabl)] <- tabl
      N  <- maf[i, "1"] + maf[i, "2"] + maf[i, "3"]
      maf[i, "Total"] <- N
      maf[i, "call.rate"] <- 1 -(maf[i, "0"] / (N + maf[i, "0"]))
      p  <- (2 * maf[i, "1"] + maf[i, "3"]) / N / 2
      maf[i, "maf"] <- ifelse(p <= 0.5, p, 1-p)
      maf[i, "minor.allele"] <- ifelse(maf[i, "1"] <= maf[i, "2"], 1, 2)
      # hwe
      observed <- maf[i, c("1", "2", "3")]
      n1  <- (2*maf[i, "1"] + maf[i, "3"])^2 / 4 / N
      n2  <- (2*maf[i, "2"] + maf[i, "3"])^2 / 4 / N
      n3  <- (2*maf[i, "1"] + maf[i, "3"]) * 
        (2 * maf[i, "2"] + maf[i, "3"]) / (2 * N) 
      expexted <-  c(n1, n2, n3)
      hwe.chisq <- sum((observed - expexted)^2 / expexted)
      maf[i, "hwe.chisq.p.value"] <- 1-pchisq(hwe.chisq, 1)
    }
    maf <- maf[, c("1", "3", "2", "Total", "0", "call.rate", 
                   "minor.allele", "maf", "hwe.chisq.p.value")]
    colnames(maf) <- c("1/1", "1/2", "2/2", "Total", "NMISS", 
                       "call.rate", "minor.allele",
                       "maf", "hwe.chisq.p.value")
    if(!(all(is.na(marker.label)))) {
      rownames(maf) <- marker.label
    }
    return(maf)
  }
