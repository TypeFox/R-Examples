##generic
modavgShrink <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                          ...){
  cand.set <- formatCands(cand.set)
  UseMethod("modavgShrink", cand.set)
}



##default
modavgShrink.default <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                  nobs = NULL,  uncond.se = "revised", conf.level = 0.95, 
                                  ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov
modavgShrink.AICaov.lm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                  nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                  ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##betareg
modavgShrink.AICbetareg <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                    nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                    ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract labels
    ##determine if parameter is on mean or phi
    if(regexpr(pattern = "\\(phi\\)_", parm) == "-1") {
      parm.phi <- NULL
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients$mean))
    } else {
      ##replace parm
      parm.phi <- gsub(pattern = "\\(phi\\)_", "", parm)
      if(regexpr(pattern = ":", parm) != "-1") {
        warning(cat("\nthis function does not yet support interaction terms on phi:\n",
                    "use 'modavgCustom' instead\n"))
      }
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients$precision))
    }
  
  
  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) stop("\nTo compute a shrinkage version of model-averaged estimate, each term must appear with the same frequency across models\n")


  ##check whether parm is involved in interaction
  ##if parameters on mean
  if(is.null(parm.phi)) {
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  }
  

  ##if parameters on phi
  if(!is.null(parm.phi)) {
    parm.inter <- c(paste(parm.phi, ":", sep = ""), paste(":", parm.phi, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  }

  
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##clm
modavgShrink.AICsclm.clm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL,
           uncond.se = "revised", conf.level = 0.95, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }


  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
  unique.link <- unique(x = check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                   "from models using different link functions\n")


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$beta))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set=cand.set, modnames=modnames,
                      second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##clmm
modavgShrink.AICclmm <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }


  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
  unique.link <- unique(x = check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                   "from models using different link functions\n")


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$beta))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##coxme
modavgShrink.AICcoxme <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                      nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                      ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

  
    ##check for frequency of each terms    
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) names(fixef(i))) #extract model formula for each model in cand.set

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    terms.freq <- table(pooled.terms)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

  
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


    ##compute table
    new_table <- aictab(cand.set = cand.set, modnames = modnames,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) extractSE(i)[paste(parm)]))
  

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
    
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 

      
    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
    
    }
  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter" =  paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##coxph and clogit
modavgShrink.AICcoxph <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                      nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                      ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  
    ##check for frequency of each terms    
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    terms.freq <- table(pooled.terms)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

  
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


    ##compute table
    new_table <- aictab(cand.set = cand.set, modnames = modnames,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
    
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 

      
    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
    
    }
  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter" =  paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)

  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##glm
modavgShrink.AICglm.lm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, gamdisp = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  #check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$family$link))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model averaged beta estimate\n",
"from models using different link functions\n")

  
  ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
  fam.type <- unlist(lapply(cand.set, FUN = function(i) family(i)$family))
  fam.unique <- unique(fam.type)
  if(identical(fam.unique, "gaussian")) {disp <- NULL} else{disp <- 1}
  ##poisson and binomial defaults to 1 (no separate parameter for variance)

  ##for negative binomial - reset to NULL
  if(any(regexpr("Negative Binomial", fam.type) != -1)) {
    disp <- NULL
    ##check for mixture of negative binomial and other
    ##number of models with negative binomial
    negbin.num <- sum(regexpr("Negative Binomial", fam.type) != -1)
    if(negbin.num < length(fam.type)) {
      stop("Function does not support mixture of negative binomial with other distributions in model set")
    }
  }  
  ##gamma is treated separately

  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
    
    
  ##check for frequency of each terms    
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                   mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE,
                      c.hat = c.hat)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i, dispersion = disp)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }
  
    
  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <-new_table$SE*sqrt(c.hat)
  } 

  gam1 <- unlist(lapply(cand.set, FUN=function(i) family(i)$family[1]=="Gamma")) #check for gamma regression models
  ##correct SE's for estimates of gamma regressions
  if(any(gam1) == TRUE)  {
    ##check for specification of gamdisp argument
    if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
    new_table$SE <- unlist(lapply(cand.set,
                                  FUN = function(i) sqrt(diag(vcov(i, dispersion = gamdisp)))[paste(parm)]))
  } 

  

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

    
    
  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     



  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }



  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE,
                     "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)

}



##gls
modavgShrink.AICgls <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

  
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

  
    ##check for frequency of each terms    
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

    
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


    ##compute table
    new_table <- aictab(cand.set = cand.set, modnames = modnames,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))
    
    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }
    
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
    
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
    
    }
  
    zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter" =  paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
  
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##hurdle
modavgShrink.AIChurdle <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                   "from models using different link functions\n")
  
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
    
    
  ##check for frequency of each terms    
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i))) #extract model formula for each model in cand.set

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "count_(Intercept)" & pooled.terms != "zero_(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                   mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN=function(i) coefficients(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##AIC
  if(second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }



  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE,
                     "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)

  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
  
}



##lm
modavgShrink.AIClm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   ...){


  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##lme
modavgShrink.AIClme <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)


    ##check for frequency of each terms 
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients$fixed))

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


    ##compute table
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }

                      
    }
  
    zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##lmekin
modavgShrink.AIClmekin <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  

  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))


  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) extractSE(i)[paste(parm)]))               
 

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##maxlike
modavgShrink.AICmaxlikeFit.list <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


    ##check that link function is the same for all models
    check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                     "from models using different link functions\n")

    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)


    ##check for frequency of each terms 
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) names(coef(i)))
   

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE,
                      c.hat = 1)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##mer
modavgShrink.AICmer <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


    ##determine families of model
    fam.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$family))
    check.fam <- unique(fam.list)
    if(length(check.fam) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                   "from models using different families of distributions\n")
    ##determine link functions
    link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
    check.link <- unique(link.list)
    if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                    "from models using different link functions\n")
###################       

    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
  
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

  
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")



    ##compute table
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) extractSE(i)[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord==TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
    
    
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower_CL<-Modavg_beta-zcrit*Uncond_SE
  Upper_CL<-Modavg_beta+zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##glmerMod
modavgShrink.AICglmerMod <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      }
      modnames <- names(cand.set)
    }
    

    ##determine families of model
    fam.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$family))
    check.fam <- unique(fam.list)
    if(length(check.fam) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                   "from models using different families of distributions\n")
    ##determine link functions
    link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
    check.link <- unique(link.list)
    if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                    "from models using different link functions\n")
###################       

    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
    
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

  
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")



    ##compute table
    new_table <- aictab(cand.set = cand.set, modnames = modnames,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) extractSE(i)[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord==TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 


    } else {
      Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
    
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
    
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    
    
    }
  
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL<-Modavg_beta-zcrit*Uncond_SE
    Upper_CL<-Modavg_beta+zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
  
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  
  }



##lmerMod
modavgShrink.AIClmerMod <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      }
      modnames <- names(cand.set)
    }
    
###################       

  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
  
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

  
  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                      mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")



  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) extractSE(i)[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord==TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
    
    
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower_CL<-Modavg_beta-zcrit*Uncond_SE
  Upper_CL<-Modavg_beta+zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
  
}



##multinom
modavgShrink.AICmultinom.nnet <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }



    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)
  
 
    ##extract model formula for each model in cand.set    
    mod_formula <- lapply(cand.set, FUN = function(i) colnames(summary(i)$coefficients)) 

    nmods <- length(cand.set)
  
  
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

    
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                      pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


    ##determine number of levels - 1
    mod.levels <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients)) #extract level of response variable 
    check.levels <- unlist(unique(mod.levels))


    ##recompute AIC table and associated measures
    new_table <- aictab(cand.set = cand.set, modnames = modnames,
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat) 

  ##create object to store model-averaged estimate and SE's of k - 1 level of response
  out.est <- matrix(data = NA, nrow = length(check.levels), ncol = 4)
  colnames(out.est) <- c("Mod.avg.est", "Uncond.SE", "Lower.CL", "Upper.CL")
  rownames(out.est) <- check.levels

  ##iterate over levels of response variable
  for (g in 1:length(check.levels)) {

    ##extract coefficients from each model for given level
    coefs.levels <- lapply(cand.set, FUN = function(i) coef(i)[check.levels[g], ])
    ##extract coefficients from each model for all levels
    SE.all.levels <- lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i))))
    id.coef <- lapply(coefs.levels, FUN = function(i) which(names(i) == paste(parm)))
    ##temporary matrix to hold estimates and SE's from models and set to 0 otherwise
    tmp.coef <- matrix(NA, ncol = 2, nrow = nmods)
    for(k in 1:nmods) {
          tmp.coef[k, 1] <- ifelse(length(id.coef[[k]]) != 0, coefs.levels[[k]][paste(parm)], 0)
          tmp.coef[k, 2] <- ifelse(length(id.coef[[k]]) != 0, SE.all.levels[[k]][paste(check.levels[g], ":",
                                    parm, sep="")], 0)
        }
          
    
    ##extract beta estimate for parm
    new_table$Beta_est <- tmp.coef[, 1]
    
    ##extract SE of estimate for parm
    new_table$SE <- tmp.coef[, 2]

    
    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

      
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {new_table$SE <- new_table$SE*sqrt(c.hat)} 

    ##compute model-averaged estimates, unconditional SE, and 95% CL
    #AICc
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    #QAICc
    #if c-hat is estimated compute values accordingly and adjust table names
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 
    }

    
    #AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }
    

    #QAIC
    #if c-hat is estimated compute values accordingly and adjust table names  
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 
    }
    

    out.est[g, 1] <- Modavg_beta 
    out.est[g, 2] <- Uncond_SE 
  }
     
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  out.est[,3] <- out.est[,1] - zcrit*out.est[,2]
  out.est[,4] <- out.est[,1] + zcrit*out.est[,2]
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = out.est[,1],
                     "Uncond.SE" = out.est[,2], "Conf.level" = conf.level, "Lower.CL"= out.est[,3],
                     "Upper.CL" = out.est[,4])

  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)

}



##polr
modavgShrink.AICpolr <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  
  ##extract model formula for each model in cand.set    
  mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) 

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[attr(regexpr(pattern = "\\|", text = pooled.terms), "match.length") == -1 ]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                   mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  
  ##add logical test to distinguish between intercepts and other coefs
  if(attr(regexpr(pattern = "\\|", text = parm), "match.length") == -1) {
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)]))
  } else {new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) (i)$zeta[paste(parm)])) }
        
  ##extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  
  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL based on AICc
  if(second.ord == TRUE) {
    Modavg_beta<-sum(new_table$AICcWt*new_table$Beta_est)
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }   
  }

  
  ##AICc  
  ##compute model-averaged estimates, unconditional SE, and 95% CL based on AIC
  if(second.ord == FALSE) {
    Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL,
                     "Upper.CL" = Upper_CL)
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##rlm
modavgShrink.AICrlm.lm <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##survreg
modavgShrink.AICsurvreg <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                                    nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                    ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  
  ##check that distribution is the same for all models
  check.dist <- sapply(X = cand.set, FUN = function(i) i$dist)
  unique.dist <- unique(x = check.dist)
  if(length(unique.dist) > 1) stop("\nFunction does not support model-averaging estimates from different distributions\n")

  
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)


  ##check for frequency of each terms 
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) names(summary(i)$coefficients))

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }

                      
  }
  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
}



##vglm
modavgShrink.AICvglm <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         c.hat = 1, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i@family@blurb[3]))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model averaged beta estimate\n",
"from models using different link functions\n")

  
  ##check family of vglm to avoid problems when requesting predictions with argument 'dispersion'
  fam.type <- unlist(lapply(cand.set, FUN=function(i) i@family@vfamily))
  fam.unique <- unique(fam.type)
  if(identical(fam.unique, "gaussianff")) {disp <- NULL} else{disp <- 1}
  if(identical(fam.unique, "gammaff")) stop("\nGamma distribution is not supported yet\n")
  ##poisson and binomial defaults to 1 (no separate parameter for variance)
  ##for negative binomial - reset to NULL
  if(identical(fam.unique, "negbinomial")) {disp <- NULL}

    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

    
  ##check for frequency of each terms    
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i))) #extract model formula for each model in cand.set

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)" & pooled.terms != "(Intercept):1" & pooled.terms != "(Intercept):2")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction or if label changes for some models - e.g., ZIP models
  ##if : not already included
  if(regexpr(":", parm, fixed = TRUE) == -1){
    
    ##if : not included
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) warning("\nLabel of parameter of interest seems to change across models:\n",
                                     "check model syntax for possible problems\n")
  } else {
    ##if : already included
    ##remove : from parm
    simple.parm <- unlist(strsplit(parm, split = ":"))[1]

    ##search for simple.parm and parm in model formulae
    no.colon <- sum(ifelse(attr(regexpr(simple.parm, mod_formula, fixed = TRUE), "match.length") != "-1", 1, 0))
    with.colon <- sum(ifelse(attr(regexpr(parm, mod_formula, fixed = TRUE), "match.length") != "-1", 0, 1))
    
    ##check if both are > 0
    if(no.colon > 0 && with.colon > 0) warning("\nLabel of parameter of interest seems to change across models:\n",
                                               "check model syntax for possible problems\n")
  }

  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN=function(i) coefficients(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN=function(i) sqrt(diag(vcov(i, dispersion = disp)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
  ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
  if(c.hat > 1) {
    new_table$SE <-new_table$SE*sqrt(c.hat)
  } 

  ##gam1 <- unlist(lapply(cand.set, FUN=function(i) family(i)$family[1]=="Gamma")) #check for gamma regression models
  ##correct SE's for estimates of gamma regressions
  ##if(any(gam1) == TRUE)  {
  ##check for specification of gamdisp argument
  ## if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
  ## new_table$SE <- unlist(lapply(cand.set,
  ##                             FUN = function(i) sqrt(diag(vcov(i, dispersion = gamdisp)))[paste(parm)]))
  ##} 
    
    

  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(c.hat == 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }
    
    
    
  ##QAICc
  if(c.hat > 1 && second.ord == TRUE) {
    Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {      
      Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }     



  ##AIC
  if(c.hat == 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }
    


  ##QAIC
  if(c.hat > 1 && second.ord == FALSE) {
    Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }  
  }     

  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE,
                     "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
  
}



##zeroinfl
modavgShrink.AICzeroinfl <-
function(cand.set, parm, modnames = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  
  ##check that link function is the same for all models
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                   "from models using different link functions\n")
  
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
    
    
  ##check for frequency of each terms    
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i))) #extract model formula for each model in cand.set

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "count_(Intercept)" & pooled.terms != "zero_(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                   mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                      second.ord = second.ord, nobs = nobs, sort = FALSE)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(cand.set, FUN=function(i) coefficients(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord == TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }


  ##AIC
  if(second.ord == FALSE) {
    Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }



  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE,
                     "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)

  class(out.modavg) <- c("modavgShrink", "list")
  return(out.modavg)
  
}



####added functionality for reversing parameters
##added additional argument parm.type = "psi", "gamma", "epsilon", "lambda", "omega", "detect"
##model type:  parameters labeled in unmarked - parameters labeled in AICcmodavg.unmarked
##single season: state, det - USE psi, detect
##multiseason model:  psi, col, ext, det - USE psi, gamma, epsilon, detect
##RN heterogeneity model: state, det - USE lambda, detect
##N-mixture: state, det - USE lambda, detect
##Open N-mixture: lambda, gamma, omega, det - USE lambda, gamma, omega, detect
##distsamp: state, det - USE lambda, detect
##gdistsamp: state, det, phi - USE lambda, detect, phi
##false-positive occupancy: state, det, fp - USE psi, detect, fp
##gpcount: lambda, phi, det - USE lambda, phi, detect
##gmultmix: lambda, phi, det - USE lambda, phi, detect
##multinomPois: state, det - USE lambda, detect



##occu
modavgShrink.AICunmarkedFitOccu <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    ##single-season occupancy model
    ##psi
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm.unmarked <- "psi"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    

  
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")


    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##colext
modavgShrink.AICunmarkedFitColExt <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    ##multiseason occupancy model
    ##psi - initial occupancy
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$psi)))
      ##create label for parm
      parm.unmarked <- "psi"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##gamma - extinction
    if(identical(parm.type, "gamma")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$col)))
      ##create label for parm
      parm.unmarked <- "col"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##epsilon - extinction
    if(identical(parm.type, "epsilon")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$ext)))
      ##create label for parm
      parm.unmarked <- "ext"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    

    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

    
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##occuRN
modavgShrink.AICunmarkedFitOccuRN <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    ##Royle-Nichols heterogeneity model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm.unmarked <- "lam"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    

    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

    
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##pcount
modavgShrink.AICunmarkedFitPCount <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    ##single season N-mixture model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      ##create label for parm
      parm.unmarked <- "lam"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

    
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##pcountOpen
modavgShrink.AICunmarkedFitPCO <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    
    ##open version of N-mixture model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
      parm.unmarked <- "lam"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##gamma - recruitment
    if(identical(parm.type, "gamma")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$gamma)))
      ##create label for parm
      parm.unmarked <- "gam"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##omega - apparent survival
    if(identical(parm.type, "omega")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$omega)))
      ##create label for parm
      parm.unmarked <- "omega"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")

    
    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##distsamp
modavgShrink.AICunmarkedFitDS <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   

    ##Distance sampling model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm.unmarked <- "lam"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      if(identical(parm.type, "detect")) {
        stop("\nModel-averaging estimates of detection covariates not yet supported for unmarkedFitDS class\n")
      }
    }
    
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")
    

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##gdistsamp
modavgShrink.AICunmarkedFitGDS <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    ##Distance sampling model with availability
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
      parm.unmarked <- "lambda"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      if(identical(parm.type, "detect")) {
        stop("\nModel-averaging estimates of detection covariates not yet supported for unmarkedFitGDS class\n")
      }
    }
    ##availability
    if(identical(parm.type, "phi")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
      parm.unmarked <- "phi"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
  

    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")
    

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##occuFP
modavgShrink.AICunmarkedFitOccuFP <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    
    ##single-season false-positive occupancy model
    ##psi
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm.unmarked <- "psi"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##false positives - fp
    if(identical(parm.type, "fp")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$fp)))
      parm.unmarked <- "fp"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
  
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")
    

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##multinomPois
modavgShrink.AICunmarkedFitMPois <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   

    ##multinomPois model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm.unmarked <- "lambda"
      ##create label for parm
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    
    
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")
    

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##gmultmix
modavgShrink.AICunmarkedFitGMM <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

    
    ##gmultmix model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
      parm.unmarked <- "lambda"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##availability
    if(identical(parm.type, "phi")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
      parm.unmarked <- "phi"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
  
    
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")
    

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }



##gpcount
modavgShrink.AICunmarkedFitGPC <- 
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           c.hat = 1, parm.type = NULL, ...){

    ##note that parameter is referenced differently from unmarked object - see labels( )
    
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgShrink for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    
    ##gpcount model
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
      parm.unmarked <- "lambda"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm.unmarked <- "p"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    ##availability
    if(identical(parm.type, "phi")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
      parm.unmarked <- "phi"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
    
    
    
    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) warning("\nVariables do not appear with same frequency across models, proceed with caution\n")
    

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab(cand.set = cand.set, modnames = modnames, 
                        second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavgShrink", "list")
    return(out.modavg)
  }




print.modavgShrink <-
  function(x, digits = 2, ...) {
    ic <- colnames(x$Mod.avg.table)[3]
    cat("\nMultimodel inference on \"", x$Parameter, "\" based on ", ic, "\n", sep = "")
    cat("\n", ic, " table used to obtain model-averaged estimate with shrinkage:\n", sep = "")
    oldtab <- x$Mod.avg.table
    if (any(names(oldtab) == "c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n", sep = "")}
    cat("\n")
    if (any(names(oldtab)=="c_hat")) {
      nice.tab <- cbind(oldtab[, 2], oldtab[, 3], oldtab[, 4], oldtab[, 6],
                        oldtab[, 9], oldtab[, 10])
    } else {nice.tab <- cbind(oldtab[, 2], oldtab[, 3], oldtab[, 4], oldtab[, 6],
                              oldtab[, 8], oldtab[, 9])
          }

##modify printing style if multinomial model is used  
    if(length(x$Mod.avg.beta) == 1) {
      colnames(nice.tab) <- c(colnames(oldtab)[c(2, 3, 4, 6)], "Estimate", "SE")
      rownames(nice.tab) <- oldtab[, 1]
      print(round(nice.tab, digits = digits))
      cat("\nModel-averaged estimate with shrinkage:", eval(round(x$Mod.avg.beta, digits = digits)), "\n")
      cat("Unconditional SE:", eval(round(x$Uncond.SE, digits = digits)), "\n")
      cat("",x$Conf.level*100, "% Unconditional confidence interval: ", round(x$Lower.CL, digits = digits),
          ", ", round(x$Upper.CL, digits = digits), "\n\n", sep = "")
    } else {
      col.ns <- ncol(nice.tab)
      nice.tab <- nice.tab[,-c(col.ns - 1, col.ns)]
      colnames(nice.tab) <- c(colnames(oldtab)[c(2, 3, 4, 6)])
      rownames(nice.tab) <- oldtab[, 1]
      print(round(nice.tab, digits = digits))
      cat("\n\nModel-averaged estimates with shrinkage for different levels of response variable:", "\n\n")
      resp.labels <- labels(x$Mod.avg.beta)
      mult.out <- matrix(NA, nrow = length(resp.labels), ncol = 4)
      colnames(mult.out) <- c("Model-averaged estimate with shrinkage", "Uncond. SE", paste(x$Conf.level*100,"% lower CL", sep = ""),
                              paste(x$Conf.level*100, "% upper CL", sep = ""))
      rownames(mult.out) <- resp.labels
      mult.out[, 1] <- round(x$Mod.avg.beta, digits = digits)
      mult.out[, 2] <- round(x$Uncond.SE, digits = digits)
      mult.out[, 3] <- round(x$Lower.CL, digits = digits)
      mult.out[, 4] <- round(x$Upper.CL, digits = digits)
      print(mult.out)
      cat("\n")
    }
  }
