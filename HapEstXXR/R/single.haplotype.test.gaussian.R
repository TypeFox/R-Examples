single.haplotype.test.gaussian <-
  function(snps, trait, adj.var = NA, lim = 0.05, 
           baseline.hap = "max", 
           do.hap.specific.test = TRUE, 
           min.count = 10, alpha = 0.05) {
    
    # check number of snps    
    snps <- as.matrix(snps)
    # Infer haplotypes
    hapest <- hapest.gaussian(snps, trait, lim = lim)
    desres <- as.matrix(hapest$desres)
    if (all(is.na(hapest$haplotypes))) {
      return(list(haplotypes = NA, nSubj = NA, df = NA, 
                  global.p.value = NA))
    }
    # check, if "rest" < min.count
    if (any(colnames(desres) == "R")) {
      # skipping "rest" from design matrix when "colsum < min.count"
      if (sum((desres[, colnames(desres) == "R"])) < min.count) {
        desres <- desres [, colnames(desres) != "R", drop = FALSE]
      }
    }
    # select baseline haplotype and drop from design matrix
    # if there is only one haplotype then don't drop.
    if (length(hapest$haplotypes[, "Pool"]) > 1) {
      baseline <- which.max(hapest$haplotypes[, "Pool"])
      desres <- desres[, 
                       -c(which(hapest$haplotypes[baseline, "Hap"] == 
                                  colnames(desres))), drop = FALSE]            
    }
    # generalized linear models
    fit <- multi.snp.test(trait, desres, adj.var, type = "gaussian")
    ## find overall p-value(goodness of fit)
    df.model <-(fit$aov.glm)$Df
    nind <-  length((fit$fit.glm1)$residuals)
    df <- as.integer(df.model[!is.na(df.model)])
    pval.model <-(fit$aov.glm)$`Pr(>F)`
    pval <- as.numeric(pval.model[!is.na(pval.model)])
    aic <- AIC(fit$fit.glm1)
    if (do.hap.specific.test == TRUE) {
    # haplotype-secific test; one haplotype or haplotype pair on time
    desres  <- hapest$desres[, colnames(hapest$desres) != "R", drop = FALSE]
    ntest   <- dim(desres)[2]
    pvali   <- rep(-1, ntest)
    aici    <- rep(-1, ntest)
    betai   <- rep(-1, ntest)
    betasdi <- rep(-1, ntest)
    for(i in 1:ntest) {
      fiti <- multi.snp.test(trait, desres[, i, drop = FALSE], adj.var, 
                             type = "gaussian")
      pval.model <-(fiti$aov.glm)$`Pr(>F)`
      pvali[i]   <- as.numeric(pval.model[!is.na(pval.model)])
      aici[i]    <- AIC(fiti$fit.glm1)
      betai[i]   <- as.numeric(coefficients(fiti$fit.glm1)[2])
      betasdi[i] <- as.numeric(
        (coefficients(summary(fiti$fit.glm1)))[2, "Std. Error"])
    }
    names(pvali)<- colnames(desres)
    names(aici) <- names(pvali)
    names(betai) <- names(pvali)
    names(betasdi) <- names(pvali)
    ## order of pvali
    ord <- match(hapest$haplotype[, 1], names(pvali))
    haplotype.i <- data.frame(
      beta             = betai[ord], 
      beta.err         = betasdi[ord], 
      "p value"        = pvali[ord], 
      "aic"            = aici[ord], 
      stringsAsFactors = FALSE)
    } else {
      haplotype.i <- NA
    }
    res.list <- list(haplotypes = hapest$haplotypes, 
                     nSubj = nind, 
                     df = df, 
                     global.p.value = pval,
                     global.aic = aic, 
                     haplotype.i = haplotype.i)
    return(res.list)
  }
