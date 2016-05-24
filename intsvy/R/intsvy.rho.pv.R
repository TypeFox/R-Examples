intsvy.rho.pv <- 
function(variable, pvlabels, by, data, export=FALSE, name= "output", folder=getwd(), config) {
  rho.pv.input <- function(variable, pvlabels, data, config) {
    # BRR / JK
    if (config$parameters$weights == "BRR") {
      stop("Not implemented yet")
    } 
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS
      
      if (length(pvlabels)==2 & missing(variable)) {
        # PV names
        pvnames <- lapply(pvlabels, function(x) paste(x, "0", 1:5, sep=""))
        # Complete dataset (listwise deletion)
        data <- na.omit(data[c(unlist(pvnames), config$variables$weight, config$variables$jackknifeRep, config$variables$jackknifeZone)])
        # Replicate weighted correlations for PV1 (sampling error)
        rhopvrp <- lapply(1:max(data[[config$variables$jackknifeZone]]), function(i) cov.wt(x=data[c(pvnames[[1]][1], pvnames[[2]][1])], cor=T, 
                  wt=ifelse(data[[config$variables$jackknifeZone]] == i, 2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], data[[config$variables$weight]]))[[5]])
        # Total weighted correlation for imputation variance
        rhopvtot <- lapply(1:5, function(i) cov.wt(x=data[c(pvnames[[1]][i],pvnames[[2]][i])], cor=T, 
                 wt= data[[config$variables$weight]])[[5]])
        # Sampling variance, imputation variance, and SEs
        varw <- Reduce("+", lapply(rhopvrp, function(x) (x - rhopvtot[[1]])^2))
        varb <- (1+1/5)* apply(simplify2array(rhopvtot), c(1, 2), var) # slower, but Reduce(var) fails
        rhose <- (varw+varb)^(1/2)
        # Mean total weighted correlation
        rhotot <- Reduce("+", rhopvtot)/length(rhopvtot)
        # Combined rhos and SEs
        rhomat <- round(do.call(cbind, lapply(1:ncol(rhotot), function(x) t(rbind(rhotot[,x], rhose[, x])))), 5)
        # Format and print
        colnames(rhomat) <- unlist(lapply(1:2, function(i) 
          c(paste(pvlabels, "Rho", sep=" ")[i], paste(pvlabels, "s.e.", sep=" ")[i])))
        return(round(rhomat, 3))
      } else {
        
        # Correlation of no PV with PV
        # PV names
        pvnames <- paste(pvlabels, "0", 1:5, sep="")
        # Complete dataset (listwise deletion)
        data <- na.omit(data[c(variable, pvnames, config$variables$weight, config$variables$jackknifeRep, config$variables$jackknifeZone)])
        # Replicate weighted correlations for PV1 (sampling error)
        rhopvrp <- lapply(1:max(data[[config$variables$jackknifeZone]]), function(i) cov.wt(x=data[c(variable, pvnames[1])], cor=T, 
                        wt=ifelse(data[[config$variables$jackknifeZone]] == i, 2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], data[[config$variables$weight]]))[[5]])
        # Total weighted correlation for imputation variance
        rhopvtot <- lapply(pvnames, function(i) cov.wt(x=data[c(variable,i)], cor=TRUE, wt= data[[config$variables$weight]])[[5]])
        # Sampling variance, imputation variance, and SEs
        varw <- Reduce("+", lapply(rhopvrp, function(x) (x - rhopvtot[[1]])^2))
        varb <- (1+1/5)* apply(simplify2array(rhopvtot), c(1, 2), var) # slower, but Reduce(var) fails
        rhose <- (varw+varb)^(1/2)
        # Mean total weighted correlation
        rhotot <- Reduce("+", rhopvtot)/length(rhopvtot)
        # Combined rhos and SEs
        rhomat <- round(do.call(cbind, lapply(1:ncol(rhotot), function(x) t(rbind(rhotot[,x], rhose[, x])))), 5)
        # Format and print
        colnames(rhomat) <- unlist(lapply(1:2, function(i) 
          c(paste(c(variable, pvlabels), "Rho", sep=" ")[i], paste(c(variable, pvlabels), "SE", sep=" ")[i])))
        return(round(rhomat, 3))
      }
      
    }
    if (config$parameters$weights == "mixed_piaac") {
      stop("Not implemented yet")
    } 
  }
  
  # If by no supplied, calculate for the complete sample    
  if (missing(by)) { 
    output <- rho.pv.input(variable=variable, pvlabels=pvlabels,  data=data, config=config)
  } else {
    if (length(pvlabels)==2 & missing(variable)) {
      output <- lapply(split(data, factor(data[[by]])), function(x) rho.pv.input(pvlabels=pvlabels,  data=x, config=config))
    } else {
      output <- lapply(split(data, factor(data[[by]])), function(x) rho.pv.input(variable=variable, pvlabels=pvlabels, data=x, config=config))
    }
  }
  
  if (export)  {
    write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
  }

  class(output) <- c("intsvy.rho", class(output))
  return(output)
  
}

