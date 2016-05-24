intsvy.log <- 
function(y, x, by, data, export=FALSE, name= "output", folder=getwd(), config) {
  
  regform <- paste(y, "~", paste(x, collapse="+"))

  reg.input <- function(y, x, data, config) {
    # If no variability in y or x, or if all missing print NA
    if (any(sapply(data[c(y, x)], function(i) all(duplicated(i))))) {
      results <- list("replicates"=NA, "residuals"= NA, "reg"=NA)
      return(results)
    }

    # BRR / JK
    if (config$parameters$weights == "BRR") {
      # balanced repeated replication
      # Replicate weighted %s (sampling error)
      # in PISA
      
      # Replicate weighted coefficients, normalised weights
      coef.rp <- suppressWarnings(lapply(1:config$parameters$BRRreps, 
                        function(i) summary(glm(formula=as.formula(regform), family=quasibinomial("logit"), 
                                    weights=nrow(data)*data[[paste0(config$variables$weightBRR, i)]]/sum(data[[paste0(config$variables$weightBRR, i)]]),
                                                      data=data))))

      # Retrieving coefficients
      rp.coef <- sapply(1:config$parameters$BRRreps, function(i) coef.rp[[i]]$coefficients[,1])
      # Total weighted regressions 
      reg.pv <- suppressWarnings(summary(glm(formula=as.formula(regform), family=quasibinomial("logit"), 
                            weights=nrow(data)*data[[config$variables$weightFinal]]/sum(data[[config$variables$weightFinal]]), data=data)))
      
      # Total weighted coefficients
      tot.coef <- reg.pv$coefficients[, 1]
      # Sampling error 
      coef.se <- (0.05*apply((rp.coef-tot.coef)^2, 1, sum))^(1/2)
      t.stat <- tot.coef/coef.se
      
      # Odds ratios and confidence intervals
      OR <- exp(tot.coef)
      
      # OR confidence intervals 
      CI95low <- exp(tot.coef - 1.96*coef.se)
      CI95up <- exp(tot.coef + 1.96*coef.se)
      
      # Table with estimates
      log.tab <- data.frame("Coef."=tot.coef, "Std. Error"=coef.se, "t value"=t.stat, 
                            as.data.frame(cbind(OR, CI95low, CI95up)), check.names=F)
      
      results <- list("replicates"=t(rp.coef), "reg"=log.tab)
      return(results)
      
    } 
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS
      
      # Replicate weights
      rp.wt <- sapply(1:max(data[[config$variables$jackknifeZone]]), function(rp) 
                    ifelse(data[[config$variables$jackknifeZone]] == rp, 
                             2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], 
                                    data[[config$variables$weight]]))
      rp.wt.n <- nrow(data)*rp.wt/apply(rp.wt, 2, sum)
      # Replicate weights coefficients for sampling error
      reg.rp <- suppressWarnings(lapply(1:max(data[[config$variables$jackknifeZone]]), function(rp) 
                      summary(glm(formula=as.formula(regform), 
                                 family=quasibinomial("logit"), weights=rp.wt.n[, rp], data=data))))
      
      # Combine coefficients 
      coef.rp <- do.call("cbind", lapply(1:max(data[[config$variables$jackknifeZone]]), function(rp) 
                        reg.rp[[rp]]$coefficients[,1]))
      # Total weighted coefficient for each PV for imputation (between) error
      reg.tot <- suppressWarnings(summary(glm(formula=as.formula(regform), family=quasibinomial("logit"), 
                             weights=nrow(data)*data[[config$variables$weight]]/sum(data[[config$variables$weight]]), data=data)))
      
      # Total weighted coefficients
      coef.tot <- reg.tot$coefficients[, 1]
      # Sampling error 
      coef.se <- (0.05*apply((coef.rp-coef.tot)^2, 1, sum))^(1/2)
      t.stat <- coef.tot/coef.se
      # Odds ratios and confidence intervals
      OR <- exp(coef.tot)
      # OR confidence intervals 
      CI95low <- exp(coef.tot - 1.96*coef.se)
      CI95up <- exp(coef.tot + 1.96*coef.se)
      
      # Table with estimates
      log.tab <- data.frame("Coef."=coef.tot, "Std. Error"=coef.se, "t value"=t.stat, 
                            as.data.frame(cbind(OR, CI95low, CI95up)), check.names=F)
      
      results <- list("replicates"=coef.rp, "reg"=log.tab)
      return(results)
      
    }
    if (config$parameters$weights == "mixed_piaac") {
      # mixed design, different for different coutnries
      # PIAAC
      stop("Not implemented yet")
    } 
  }
  
  # If by no supplied, calculate for the complete sample    
  if (missing(by)) { 
    output <- reg.input(y=y, x=x, data=data, config=config)
  } else {
    output <- lapply(split(data, droplevels(data[by])), 
                     function(i) reg.input(y=y, x=x, data=i, config=config))
    
  }
  if (export)  {
    write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
  }
  class(output) <- c("intsvy.reg", class(output))
  return(output)
  
}

