single.snp.test.binomial <-
  function(snps, trait, adj.var = NULL, prt = TRUE) {
    
    snps <- as.matrix(snps)
    N <- dim(snps)[1]
    ns <- dim(snps)[2]
    adjusted <- FALSE
    if (!all(is.null(adj.var))) {
      adjusted <- TRUE
    }
    if (!all(trait[!is.na(trait)] ==  0 | trait[!is.na(trait)] ==  1))
      stop("trait should be 0 for controls or 1 for cases")
    pval <- rep(-1, ns)
    nind <- rep(-1, ns)
    aic  <- rep(-1, ns)
    aicc <- rep(-1, ns)
    beta <- rep(-1, ns)
    beta.stderr <- rep(-1, ns)
    r2 <- rep(-1, ns)
    r2.adjusted <- rep(-1, ns)
    y <- trait
    for(i in 1:ns) {
      x <- as.matrix(as.numeric(alleleRto1(snps[, i])), ncol = 1)
      miss.value <- which(is.na(x))
      if (length(miss.value) > 0) {
        fit <- multi.snp.test(y[-miss.value], 
                              x[-miss.value, , drop = F], 
                              x.adj = adj.var[-miss.value, ], 
                              type = "binomial")
      } else {
        fit <- multi.snp.test(y, 
                              x, 
                              x.adj = adj.var, 
                              type = "binomial")
      }
      nind[i]        <- length((fit$fit.glm1)$residuals)
      pval.model     <-(fit$aov.glm)$`Pr(>Chi)`
      if (all(is.na(pval.model))) {
        pval[i]        <- NA
        aic[i]         <- NA
        aicc[i]        <- NA
        beta[i]        <- NA
        beta.stderr[i] <- NA
      } else {
        pval[i]        <- pval.model[!is.na(pval.model)]
        aic[i]         <- AIC(fit$fit.glm1)
        p              <- length(coef(fit$fit.glm1))
        N              <- length((fit$fit.glm1)$residuals)
        aicc[i]        <- aic[i] + 2 * (p * (p + 1)) / (N - p -1)
        beta[i]        <- fit$beta[2, 1]
        beta.stderr[i] <- fit$beta[2, 2]
      }
    }
    xalpha <- qnorm(.975)
    upper.95 <- exp(beta + xalpha * beta.stderr)
    lower.95 <- exp(beta - xalpha * beta.stderr)
    res <- data.frame(snp = 1:ns, 
                      N = nind,
                      type = "binomial", 
                      beta = beta, 
                      "se(beta)" = beta.stderr, 
                      "OR" = exp(beta), 
                      lower.95 = lower.95, 
                      upper.95 = upper.95, 
                      test = "Chisq", 
                      p.value = pval,
                      aic = aic, 
                      aicc = aicc,
                      stringsAsFactors = FALSE)
    colnames(res) <- c("SNP", "N", "type", "beta", "se(beta)",
                       "exp(beta)", "lower.95", 
                       "upper.95", "Test", "p.value", "AIC", "AICc")
    return(res)
  }
