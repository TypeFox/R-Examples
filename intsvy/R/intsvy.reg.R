intsvy.reg <- 
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
      
      # Replicate weights coefficients for sampling error
      reg.rp <- lapply(1:config$parameters$BRRreps, function(i) summary(lm(formula=as.formula(regform), data=data, 
                      weights=data[[paste(config$variables$weightBRR, i , sep="")]])))
      # Combining coefficients and R-squared replicates
      coef.rp <- sapply(1:config$parameters$BRRreps, function(i) 
                      c(reg.rp[[i]]$coefficients[,1], "R-squared" = reg.rp[[i]]$r.squared))
      
      resid <- sapply(1:config$parameters$BRRreps, function(rep) reg.rp[[rep]]$residuals)
      
      # Total weighted coefficients and R-squared
      reg.tot <- summary(lm(formula=as.formula(regform), data=data, 
                            weights=data[[config$variables$weightFinal]]))
      coef.tot <- c(reg.tot$coefficients[,1], "R-squared" = reg.tot$r.squared)
      
      # Sampling error
      coef.se <- (0.05*apply((coef.rp-coef.tot)^2, 1, sum))^(1/2)
      # T-value
      coef.t <- coef.tot/coef.se
      # Reg Table
      reg.tab <- data.frame("Estimate"=coef.tot, "Std. Error"=coef.se, "t value"=coef.t, check.names=FALSE)
      
      results <- list("replicates"=t(coef.rp), "residuals"= resid, "reg"=reg.tab)
      return(results)
      
    } 
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS
      
      # Replicate weights coefficients for sampling error
      reg.rp <- lapply(1:max(data[[config$variables$jackknifeZone]]), function(i) summary(lm(formula= as.formula(regform), data=data, 
                     weights=ifelse(data[[config$variables$jackknifeZone]] == i, 
                                    2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], data[[config$variables$weight]]))))
      
      # Combine coefficients and R-squared
      coef.rp <- sapply(1:max(data[[config$variables$jackknifeZone]]), function(i) 
        c(reg.rp[[i]]$coefficients[,1], "R-squared" = reg.rp[[i]]$r.squared))
      
      resid <- sapply(1:length(reg.rp), function(rep) reg.rp[[rep]]$residuals)
      
      
      # Total weighted coefficients
      reg.tot <- summary(lm(formula=as.formula(regform), data=data, weights=data[[config$variables$weight]]))
      coef.tot <- c(reg.tot$coefficients[ ,1] , "R-squared" = reg.tot$r.squared)
      
      # Sampling error
      coef.se <- apply((coef.rp-coef.tot)^2, 1, sum)^(1/2)
      # T-value
      coef.t <- coef.tot/coef.se
      # Reg Table
      reg.tab <- data.frame("Estimate"=coef.tot, "Std. Error"=coef.se, "t value"=coef.t, check.names=FALSE)
      
      results <- list("replicates"=t(coef.rp), "residuals"= resid, "reg"=reg.tab)
      return(results)
    }
    if (config$parameters$weights == "mixed_piaac") {
      # mixed design, different for different coutnries
      # PIAAC
      
      # Replicate weights coefficients for sampling error
      Coefrp <- lapply(1:config$parameters$BRRreps, function(i) summary(lm(formula=as.formula(regform), data=data, 
                                                    weights=data[[paste(config$variables$weightBRR, i , sep="")]])))
      # Combining coefficients and R-squared replicates
      Statrp <- sapply(1:config$parameters$BRRreps, function(i) 
                c(Coefrp[[i]]$coefficients[,1], Coefrp[[i]]$r.squared))
      
      # Total weighted coefficients and R-squared
      Reg <- summary(lm(formula=as.formula(regform), data=data, weights=data[[config$variables$weightFinal]]))
      Stattot <- c(Reg$coefficients[,1], Reg$r.squared)
      names(Stattot)[length(Stattot)] <- "R-squared"
      
      cntName <- as.character(unique(data$CNTRYID))[1]
      cc <- piaacReplicationScheme[cntName,"c"]
      if (length(unique(piaacReplicationScheme[as.character(unique(data$CNTRYID)),"c"])) > 1) {
        warning(paste("In PIAAC study different replications schemes were applied in different countries. \n In the selected set of countries more than one scheme was used. \n Further estimation is performed with coefficient c =", cc))
      }
      if (is.na(cc)) cc <- 1
      # Sampling error
      StatSE <- (cc*apply((Statrp-Stattot)^2, 1, sum))^(1/2)
      # T-value
      StatT <- Stattot/StatSE
      # Reg Table
      RegTab <- round(data.frame("Estimate"=Stattot, "Std. Error"=StatSE, "t value"=StatT, check.names=F),2)

      results <- list("replicates"=t(Coefrp), "reg"=RegTab)
      return(results)
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

