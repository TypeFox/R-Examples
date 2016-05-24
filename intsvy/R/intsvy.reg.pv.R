intsvy.reg.pv <-
  function(x, pvlabel, by, data, std=FALSE, export=FALSE, name= "output", folder=getwd(), config) {

  reg.pv.input <- function(x, pvlabel, data, std, config) {
    if (any(sapply(data[x], function(i) all(duplicated(i))))) {
    results <- list("replicates"=NA, "residuals"= NA, "var.w"=NA, "var.b"=NA, "reg"=NA)
    return(results)
    }

    # BRR / JK
    if (config$parameters$weights == "BRR") {
      # balanced repeated replication
      # Replicate weighted %s (sampling error)
      # in PISA

      # PV labels
      pvnames <- paste("PV", 1:config$parameters$PVreps, pvlabel, sep="")
      # List of formulas for each PV
      regform <- lapply(pvnames, function(i) paste(i, "~", paste(x, collapse="+")))

      # Standardise IV and DV variables
      if(std) {
        data <-  cbind(scale(data[c(pvnames, x)]), data[!names(data) %in% c(pvnames, x)])
      }

      # Replicate weighted coefficients for sampling error (5 PVs)
      reg.rep <- lapply(regform, function(pv) lapply(1:config$parameters$BRRreps, function(rep)
        summary(lm(formula=as.formula(pv), data=data, weights=data[[paste0(config$variables$weightBRR, rep)]]))))


      # Combining coefficients and R-squared replicates
      coe.rep <- lapply(1:config$parameters$PVreps, function(pv) sapply(1:config$parameters$BRRreps, function(rep)
        c(reg.rep[[pv]][[rep]]$coefficients[,1], "R-squared"= reg.rep[[pv]][[rep]]$r.squared)))

      resid <- lapply(1:config$parameters$PVreps, function(pv)
                    sapply(1:config$parameters$BRRreps,
                          function(rep) reg.rep[[pv]][[rep]]$residuals))

      # Total weighted coefficient for each PV for imputation (between) error
      reg.pv <- lapply(regform, function(pv)
                    summary(lm(formula=as.formula(pv), data=data, weights=data[[config$variables$weightFinal]])))

      coe.tot <- sapply(1:config$parameters$PVreps, function(pv)
                    c(reg.pv[[pv]]$coefficients[, 1], "R-squared" = reg.pv[[pv]]$r.squared))


      # Mean total coefficients (across PVs)
      stat.tot <- apply(coe.tot, 1, mean)

      # Sampling error (variance within)
      var.w <- apply(0.05*sapply(lapply(1:config$parameters$PVreps, function(pv)
                    (coe.rep[[pv]]-coe.tot[,pv])^2), function(e) apply(e, 1, sum)), 1, mean)

      # Imputation error (variance between)
      var.b <- (1/4)*apply(sapply(1:config$parameters$PVreps, function(pv)
                    (coe.tot[, pv] - stat.tot)^2), 1, sum)

      stat.se <- (var.w +(1+1/config$parameters$PVreps)*var.b)^(1/2)
      stat.t <- stat.tot/stat.se

      # Reg Table
      reg.tab <- data.frame("Estimate"=stat.tot, "Std. Error"=stat.se, "t value"=stat.t, check.names=F)
      results <- list("replicates"=lapply(coe.rep, t), "residuals"= resid, "var.w"=var.w, "var.b"=var.b, "reg"=reg.tab)
      return(results)

    }
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS

      pvnames <- paste(pvlabel, "0", 1:5, sep="")
      # List of formulas for each PV
      regform <- lapply(pvnames, function(i) paste(i, "~", paste(x, collapse="+")))

      # Standardise IV and DV variables
      if(std) {
        data <-  cbind(scale(data[c(pvnames, x)]), data[!names(data) %in% c(pvnames, x)])
      }

      # Replicate weighted coefficients for sampling error (PV1 only)
      reg.rep <- lapply(1:max(data[[config$variables$jackknifeZone]]),
                        function(i) summary(lm(formula=as.formula(regform[[1]]), data=data,
                                    weights=ifelse(data[[config$variables$jackknifeZone]] == i,
                                                   2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]],
                                                   data[[config$variables$weight]]))))

      # Combining coefficients and R-squared replicates
      coe.rep <- sapply(1:max(data[[config$variables$jackknifeZone]]), function(i)
                    c(reg.rep[[i]]$coefficients[,1], "R-squared"= reg.rep[[i]]$r.squared))

      resid <- sapply(1:length(reg.rep), function(rep) reg.rep[[rep]]$residuals)

      # Total weighted coefficient for each PV for imputation (between) error
      reg.pv <- lapply(regform, function(i)
              summary(lm(formula=as.formula(i), data=data, weights=data[[config$variables$weight]])))
      coe.tot <- sapply(1:config$parameters$PVreps, function(pv)
              c(reg.pv[[pv]]$coefficients[, 1], "R-squared" = reg.pv[[pv]]$r.squared))

      # Mean total coefficients (across PVs)
      stat.tot <- apply(coe.tot, 1, mean)

      # Sampling error for PV1 (variance within)
      var.w <- apply((coe.rep-coe.tot[,1])^2, 1, sum)

      # Imputation error (variance between)
      var.b <- (1+1/config$parameters$PVreps)*apply(coe.tot, 1, var)
      stat.se <- (var.w + var.b)^(1/2)
      stat.t <- stat.tot/stat.se

      # Reg Table
      reg.tab <- data.frame("Estimate"=stat.tot, "Std. Error"=stat.se, "t value"=stat.t, check.names=F)

      results <- list("replicates"=coe.rep, "residuals"= resid, "var.w"=var.w, "var.b"=var.b, "reg"=reg.tab)
      return(results)

    }
    if (config$parameters$weights == "mixed_piaac") {
      # mixed design, different for different coutnries
      # PIAAC

      # PV labels
      pvnames <- paste("PV", pvlabel, 1:config$parameters$PVreps, sep="")
      # List of formulas for each PV
      regform <- lapply(pvnames, function(i) paste(i, "~", paste(x, collapse="+")))

      # Replicate weighted coefficients for sampling error (5 PVs)
      Coefrpv <- lapply(regform, function(k) lapply(1:config$parameters$BRRreps, function(i)
        summary(lm(formula=as.formula(k), data=data,
                   weights=data[[paste(config$variables$weightBRR, i , sep="")]]))))

      # Combining coefficients and R-squared replicates
      Statrp <- lapply(1:config$parameters$PVreps, function(pv) sapply(1:config$parameters$BRRreps, function(i)
        c(Coefrpv[[pv]][[i]]$coefficients[,1], Coefrpv[[pv]][[i]]$r.squared)))

      # Total weighted coefficient for each PV for imputation (between) error
      Regpv <- lapply(regform, function(i)
              lm(formula=as.formula(i), data=data, weights=data[[config$variables$weightFinal]]))

      Stattot <- sapply(1:config$parameters$PVreps, function(pv)
              c(summary(Regpv[[pv]])$coefficients[, 1], summary(Regpv[[pv]])$r.squared))
      rownames(Stattot)[nrow(Stattot)] <- "R-squared"

      # Mean total coefficients (across PVs)
      Stattotm <- apply(Stattot, 1, mean)

      cntName <- as.character(unique(data$CNTRYID))[1]
      cc <- piaacReplicationScheme[cntName,"c"]
      if (is.na(cc)) cc <- 1
      if (length(unique(piaacReplicationScheme[as.character(unique(data$CNTRYID)),"c"])) > 1) {
        warning(paste("In PIAAC study different replications schemes were applied in different countries. \n In the selected set of countries more than one scheme was used. \n Further estimation is performed with coefficient c =", cc))
      }
      # Sampling error (variance within)
      Varw <- apply(cc*sapply(lapply(1:config$parameters$PVreps, function(pv)
              (Statrp[[pv]]-Stattot[,pv])^2), function(e) apply(e, 1, sum)), 1, mean)

      # Imputation error (variance between)
      Varb <- (1/(config$parameters$PVreps-1))*apply(sapply(1:config$parameters$PVreps, function(i)
              (Stattot[, i] - Stattotm)^2), 1, sum)

      StatSE <- (Varw+(1+1/config$parameters$PVreps)*Varb)^(1/2)
      StatT <- Stattotm/StatSE

      # Reg Table
      RegTab <- round(data.frame("Estimate"=Stattotm, "Std. Error"=StatSE, "t value"=StatT, check.names=FALSE),2)

      results <- list("replicates"=t(Statrp), "reg"=RegTab)
      return(results)
    }
  }

  # If by no supplied, calculate for the complete sample
  if (missing(by)) {
    output <- reg.pv.input(x=x, pvlabel=pvlabel, data=data, std=std, config=config)
  } else {
    output <- lapply(split(data, droplevels(data[by])), function(i)
      reg.pv.input(x=x, pvlabel=pvlabel, data=i, std=std, config=config))
  }

  if (export)  {
    write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
  }
  class(output) <- "intsvy.reg"
  return(output)

}

