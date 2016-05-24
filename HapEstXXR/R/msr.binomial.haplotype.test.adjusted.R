msr.binomial.haplotype.test.adjusted <-
  function(snps, trait, adj.var, lim = 0.05, 
           baseline.hap="max", min.count = 10){
    adj.var <- as.matrix(adj.var)
    # Infer haplotypes
    hapest <- itegeppXXR(snps, des = 0, lim = lim)
    desres <- as.matrix(hapest$desres)
    if(all(is.na(hapest$hap))) {
      return(list(haplotypes = NA,
                  nSubj = NA,
                  df = NA, 
                  global.p.value = NA,
                  aic = NA,
                  aicc = NA))
    }
    # check, if "rest" < min.count
    if(any(colnames(desres) == "R")) {
      # skipping "rest" from design matrix when "colsum < min.count"  
      if(sum((desres[, colnames(desres) == "R"])) < min.count) {
        desres <- desres [, colnames(desres) != "R", drop = FALSE]  
      }
    }
    # select baseline haplotype and drop from design matrix
    # if there is only one haplotype then don't drop.  
    if(length(hapest$freq) > 1) {
      baseline <- which.max(hapest$freq)
      desres <- desres[, -c(which(hapest$hap[baseline] == colnames(desres))), 
                       drop = FALSE]
    }
    # generalized linear models
    # covariates 
    fit.glm0 <- glm(trait ~ adj.var, family = "binomial")
    fit.glm1 <- glm(trait ~ ., data = as.data.frame(cbind(desres, adj.var)), 
                    family = "binomial")
    aov.glm  <- anova(fit.glm0, fit.glm1, test = "Chisq")
    ## find overall p-value(goodness of fit)
    df.model <-(aov.glm)$Df
    nind <-  length(fit.glm1$residuals)
    df <- as.integer(df.model[!is.na(df.model)])
    pval.model <-(aov.glm)$`Pr(>Chi)`
    pval <- as.numeric(pval.model[!is.na(pval.model)])
    aic  <- AIC(fit.glm1)  
    p    <- length(coef(fit.glm1))
    N    <- length(fit.glm1$residuals)
    aicc <- aic + 2 * (p * (p + 1)) / (N - p - 1)
    return(list(haplotypes = 1, 
                nSubj = nind, 
                df = df,
                global.p.value = pval,
                AIC = aic,
                AICc = aicc))
  }
