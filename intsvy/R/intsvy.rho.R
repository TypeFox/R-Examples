intsvy.rho <- 
function(variables, by, data, export=FALSE, name= "output", folder=getwd(), config) {
  rho.input <- function(variables, data, config) {
    # BRR / JK
    if (config$parameters$weights == "BRR") {
      # balanced repeated replication
      # Replicate weighted %s (sampling error)
      # in PISA / PIAAC
      
      data <- na.omit(data[c(variables, config$variables$weightFinal, grep(config$variables$weightBRR, names(data), value=TRUE))]) 
      # Fifth element is correlation matrix
      rhorp <-  lapply(1:config$parameters$BRRreps, function(i) cov.wt(data[variables], wt= data[[paste(config$variables$weightBRR, i , sep="")]], cor = TRUE)[[5]])
      rhotot <- cov.wt(data[variables], wt=data[[config$variables$weightFinal]], cor=TRUE)[[5]]
      
      # SE formula
      
      # Standard error (sampling eror) 
      rhose <- (0.05*Reduce("+", lapply(rhorp, function(x) (x-rhotot)^2)))^(1/2)
      
      # Combined rhos and SEs
      rhomat <- round(do.call(cbind, lapply(1:ncol(rhotot), function(x) t(rbind(rhotot[,x], rhose[, x])))), 3)
      colnames(rhomat) <- unlist(lapply(1:length(variables), function(x) 
        c(paste(variables, "Rho", sep=" ")[x], paste(variables, "s.e.", sep=" ")[x])))
      return(round(rhomat, 6))
      
    } 
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS
      data <- na.omit(data[c(variables, config$variables$weight, config$variables$jackknifeRep, config$variables$jackknifeZone)]) 
      # Fifth element is correlation matrix
      rhorp <- lapply(1:max(data[[config$variables$jackknifeZone]]), function(i) cov.wt(data[variables],
                wt=ifelse(data[[config$variables$jackknifeZone]] == i, 2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], 
                          data[[config$variables$weight]]), cor=TRUE)[[5]])
      rhotot <- cov.wt(data[variables], wt=data[[config$variables$weight]], cor=TRUE)[[5]]
      # SE formula
      rhose <- Reduce("+", lapply(rhorp, function(x) (x-rhotot)^2))^(1/2) 
      # Combined rhos and SEs
      rhomat <- round(do.call(cbind, lapply(1:ncol(rhotot), function(x) t(rbind(rhotot[,x], rhose[, x])))), 3)
      colnames(rhomat) <- unlist(lapply(1:length(variables), function(x) 
        c(paste(variables, "Rho", sep=" ")[x], paste(variables, "s.e.", sep=" ")[x])))
      return(round(rhomat, 6))
      
    }
    if (config$parameters$weights == "mixed_piaac") {
      # mixed design, different for different coutnries
      # PIAAC
      
      stop("Not implemented yet")
    } 
  }
  
  # If by no supplied, calculate for the complete sample    
  if (missing(by)) { 
    output <- rho.input(variables=variables, data=data, config=config)
  } else {
    for (i in by) {
      data[[c(i)]] <- as.factor(data[[c(i)]])
    }
    output <- ddply(data, by, function(x) rho.input(data=x, variables=variables, config=config))
  }
  if (export)  {
    write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
  }
  class(output) <- c("intsvy.rho", class(output))
  return(output)
  
}

