single.snp.test.casecohort.prentice <-
  function(snps, trait, 
           patid, start.time, stop.time, subcohort, stratvar = NA, 
           robust, adj.var, prt, ties = "efron") {
    
    eps        <- 0.00005
    patid      <- as.character(patid)
    start.time <- as.numeric(start.time)
    stop.time  <- as.numeric(stop.time)
    N <- nrow(snps)
    adjusted <- FALSE
    if (!(is.null(adj.var))) {
      adjusted <- TRUE
      adj.var <- as.matrix(adj.var)
    }
    if (!all(is.na(stratvar))) {
      stratvar <- as.matrix(stratvar)
    }
    # Prentice:
    start.time <- ifelse((trait == 1) & (subcohort == 0), 
                         stop.time - eps, start.time)
    stop    <- NULL
    start   <- NULL
    subco   <- NULL
    tt      <- NULL
    xx      <- NULL
    id      <- NULL
    adj.mat <- NULL
    strat   <- NULL    
    for(i in 1:length(trait)) {
      # Standard Cox model counting process approach
      if (trait[i] == 1) {  # Cases
        start  <- c(start, stop.time[i] - 0.00005)
        stop   <- c(stop, stop.time[i])
        subco  <- c(subco, subcohort[i])
        tt     <- c(tt, trait[i])
        xx     <- rbind(xx, snps[i, ])
        id     <- c(id, patid[i])
        if (adjusted) {
          adj.mat <- rbind(adj.mat, adj.var[i, ])
        }
        if (!all(is.na(stratvar))) {
          strat <- rbind(strat, stratvar[i, ])
        }
        if (subcohort[i] == 1) {   # Subcohort case treated as control
          start  <- c(start, start.time[i])
          stop   <- c(stop, stop.time[i] - 0.00005)
          subco  <- c(subco, subcohort[i])
          tt     <- c(tt, 0)
          xx     <- rbind(xx, snps[i, ])
          id     <- c(id, patid[i])
          if (adjusted) {
            adj.mat <- rbind(adj.mat, adj.var[i, ])
          }
          if (!all(is.na(stratvar))) {
            strat <- rbind(strat, stratvar[i, ])
          }
        }
      } else {
        if (subcohort[i] == 1) { # Subcohort controls 
          start  <- c(start, start.time[i])
          stop   <- c(stop, stop.time[i])
          subco  <- c(subco, subcohort[i])
          tt     <- c(tt, trait[i])
          xx     <- rbind(xx, snps[i, ])
          id     <- c(id, patid[i])
          if (adjusted) {
            adj.mat <- rbind(adj.mat, adj.var[i, ])
          }
          if (!all(is.na(stratvar))) {
            strat <- rbind(strat, stratvar[i, ])
          }
        }
      }
    }
    ysurv <- Surv(start, stop, tt)
    nloci <- dim(snps)[2]
    if (robust == TRUE) {
      res <-as.data.frame(matrix(0, 
                                 nrow = nloci,
                                 ncol = 11), 
                          stringsAsFactors = FALSE)
    } else {
      res <-as.data.frame(matrix(0,
                                 nrow = nloci,
                                 ncol = 10),
                          stringsAsFactors = FALSE)
    }    
    if (!all(is.na(stratvar))) { 
      fml <- "ysurv  ~  cbind(x, adj.mat) + cluster(as.factor(1:length(id)))"
      for(ii in 1:dim(stratvar)[2]) {
        fml <- paste(fml, "+strata(strat[, ", ii, "])", collapse = "")
      }
      for(j in 1:dim(xx)[2]) {
        x <- as.matrix(as.numeric(alleleRto1(xx[, j])), ncol = 1)
        phfit <- summary(coxph(as.formula(fml), 
                               method = ties, 
                               robust = robust))  # marker and covariates 
        # chisq LR
        LR <-  phfit$logtest[1]
        # pseudo R2 Nagelkerke
        Nagelkerke.R2 <- phfit$rsq[1] / phfit$rsq[2]
        if (robust == FALSE) {
          res [j, ] <-
            c(j, phfit$n, 0, 
              as.numeric(formatC(c(
                (phfit$coefficients)[1, c("coef", "se(coef)")], 
                (phfit$conf.int)[1, c("exp(coef)", "lower .95", "upper .95")], 
                (phfit$coefficients)[1, c("Pr(>|z|)")], LR, Nagelkerke.R2), 
                format = "fg")))
        } else {
          res [j, ] <-
            c(j, phfit$n, 0, 
              as.numeric(formatC(c(
                (phfit$coefficients)[1, c("coef", "robust se")], 
                (phfit$conf.int)[1, c("exp(coef)", "lower .95", "upper .95")], 
                (phfit$coefficients)[1, c("Pr(>|z|)")], LR, Nagelkerke.R2),
                format = "fg")))
        }
      }
    } else {
      # no stratvar
      for(j in 1:dim(xx)[2]) {        
        x <- as.matrix(as.numeric(alleleRto1(xx[, j])), ncol = 1)
        phfit <- summary(coxph(ysurv  ~  cbind(x, adj.mat) + 
                                 cluster(as.factor(1:length(id))), 
                               method = ties, 
                               robust = robust))  # marker and covariates
        # chisq LR
        LR <- phfit$logtest[1]
        # pseudo R2 Nagelkerke
        Nagelkerke.R2 <- phfit$rsq[1]/phfit$rsq[2]
        if (robust == FALSE) {
          res [j, ] <-
            c(j, phfit$n, 0, 
              as.numeric(formatC(c(
                (phfit$coefficients)[1, c("coef", "se(coef)")], 
                (phfit$conf.int)[1, c("exp(coef)", "lower .95", "upper .95")], 
                (phfit$coefficients)[1, c("Pr(>|z|)")], LR, Nagelkerke.R2), 
                format = "fg")))
        } else {
          res [j, ] <-
            c(j, phfit$n, 0, 
              as.numeric(formatC(c(
                (phfit$coefficients)[1, c("coef", "robust se")], 
                (phfit$conf.int)[1, c("exp(coef)", "lower .95", "upper .95")], 
                (phfit$coefficients)[1, c("Pr(>|z|)")], LR, Nagelkerke.R2),
                format = "fg")))
        }
      }
    }
    if (robust == TRUE) {
      colnames(res) <- c("SNP", "N", "type", "beta", "se(beta)", 
                         "exp(beta)", "lower.95", "upper.95",
                         "p.value", "LR", "Nagelkerke.R2")
    } else {
      colnames(res) <- c("SNP", "N", "type", "beta", "robust se(beta)", 
                         "exp(beta)", "lower.95", "upper.95", 
                         "LR", "Nagelkerke.R2")
    }
    res[, 3] <- "case.cohort.prentice"
    return(res)
  }
