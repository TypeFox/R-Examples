single.snp.test.gaussian <-
  function(snps, trait, adj.var = NULL, prt = TRUE) {
    
    snps <- as.matrix(snps)
    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    adjusted <- FALSE
    if (!all(is.null(adj.var))) {
      adjusted <- TRUE 
    }
    pval <- rep(-1, ns)
    nind <- rep(-1, ns)
    aic  <- rep(-1, ns)
    aicc <- rep(-1, ns)
    beta <- rep(-1, ns)
    f    <- rep(-1, ns)
    beta.stderr <- rep(-1, ns)
    y <- trait
    for(i in 1:ns) {
      x <- as.matrix(as.numeric(alleleRto1(snps[, i])), ncol = 1)
      miss.value <- which(is.na(x))
      if (length(miss.value)>0) {
        fit <- multi.snp.test(y[-miss.value], 
                              x[-miss.value, , drop = FALSE], 
                              x.adj = adj.var[-miss.value, , drop = FALSE],
                              type = "gaussian")
      } else {
        fit <- multi.snp.test(y, x, x.adj = adj.var, type = "gaussian")
      }
      nind[i]    <- length((fit$fit.glm1)$residuals)
      pval.model <-(fit$aov.glm)$`Pr(>F)`
      if (all(is.na(pval.model))) {
        pval[i]        <- NA
        aic[i]         <- NA
        aicc[i]        <- NA
        beta[i]        <- NA
        f[i]           <- NA
        beta.stderr[i] <- NA
      } else {
        pval[i]        <- pval.model[!is.na(pval.model)]
        aic[i]         <- AIC(fit$fit.glm1)
        p              <- length(coef(fit$fit.glm1))
        N              <- length(residuals(fit$fit.glm1))
        aicc[i]        <- aic[i] + 2 * (p * (p + 1)) / (N - p -1)
        beta[i]        <- fit$beta[2, 1]
        beta.stderr[i] <- fit$beta[2, 2]
        f[i]           <-(fit$aov.glm)$F[!is.na((fit$aov.glm)$F)]
      }
    }
    res <- data.frame(snp = 1:ns, 
                      N = nind,
                      type = "gaussian", 
                      beta.estimate = beta, 
                      beta.stderr = beta.stderr,
                      test = "F", 
                      "F" = f, 
                      p.value = pval,
                      aic = aic, 
                      aicc = aicc,
                      stringsAsFactors = FALSE)
    colnames(res) <- c("SNP", "N", "type", "beta", "se(beta)", 
                       "Test", "F", "p.value", "AIC", "AICc")
    return(res)
  }
