intsvy.mean <- 
function(variable, by, data, export=FALSE, name= "output", folder=getwd(), config) {
  mean.input <- function(variable, data, config) {
    # BRR / JK
    if (config$parameters$weights == "BRR") {
      # balanced repeated replication
      # Replicate weighted %s (sampling error)
      # in PISA / PIAAC
      
      meanrp <- sapply(1:config$parameters$BRRreps, function(i) 
                      weighted.mean(as.numeric(data[[variable]]), 
                               data[[paste(config$variables$weightBRR, i , sep="")]], na.rm = TRUE))
      
      # Replicate weights for SDs (sampling error)
      sdrp <- sapply(1:config$parameters$BRRreps, function(i)  
                      (sum(data[[paste0(config$variables$weightBRR, i)]]*(data[[variable]]-meanrp[i])^2, na.rm = TRUE)/
                               sum(data[[paste0(config$variables$weightBRR, i)]], na.rm = TRUE))^(1/2))
      
      # Total weighted mean                                                                      
      meantot <- weighted.mean(as.numeric(data[[variable]]), data[[config$variables$weightFinal]], na.rm = TRUE)
      
      # Total weighted SD
      sdtot <-  (sum(data[[config$variables$weightFinal]]*(data[[variable]]-meantot)^2, na.rm=TRUE)/sum(data[[config$variables$weightFinal]], na.rm = TRUE))^(1/2)
      
      # Standard error (sampling eror) 
      meanse <- (0.05*sum((meanrp-meantot)^2))^(1/2)
      sdse <- (0.05*sum((sdrp-sdtot)^2))^(1/2)
      
      result <- data.frame("Freq"=sum(!is.na(data[[variable]])), "Mean"= meantot, "s.e."= meanse, "SD" = sdtot, "s.e" = sdse)
      return(round(result, 2))
      
    } 
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS

      # Replicate weight means (sampling error)
      meanrp <- sapply(1:max(data[[config$variables$jackknifeZone]]), function(x) 
        weighted.mean(as.numeric(data[[variable]]), ifelse(data[[config$variables$jackknifeZone]] == x, 
                         2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], data[[config$variables$weight]]), na.rm = TRUE))
      # Total weighted mean                                                                      
      meantot <- weighted.mean(as.numeric(data[[variable]]), data[[config$variables$weight]], na.rm = TRUE)
      # Standard error (sampling eror) 
      meanse <- (sum((meanrp-meantot)^2))^(1/2)
      result <- data.frame("Freq"=sum(!is.na(data[[variable]])), "Mean"= meantot, "s.e."= meanse)
      return(round(result, 2))
      
    }
    if (config$parameters$weights == "mixed_piaac") {
      # mixed design, different for different coutnries
      # PIAAC
      
      meanrp <- sapply(1:config$parameters$BRRreps, function(i) 
              weighted.mean(as.numeric(data[[variable]]), 
                           data[[paste(config$variables$weightBRR, i , sep="")]], na.rm = TRUE))
      
      # Total weighted mean                                                                      
      meantot <- weighted.mean(as.numeric(data[[variable]]), data[[config$variables$weightFinal]], na.rm = TRUE)
      # Standard error (sampling eror) 
      cntName <- as.character(unique(data$CNTRYID))[1]
      cc <- piaacReplicationScheme[cntName,"c"]
      if (is.na(cc)) cc <- 1
      
      meanse <- (cc*sum((meanrp-meantot)^2))^(1/2)
      result <- data.frame("Freq"=sum(!is.na(data[[variable]])), "Mean"= meantot, "s.e."= meanse)
      return(round(result, 2))

    } 
  }
  
  # If by no supplied, calculate for the complete sample    
  if (missing(by)) { 
    output <- mean.input(variable=variable, data=data, config=config)
  } else {
    for (i in by) {
      data[[c(i)]] <- as.factor(data[[c(i)]])
    }
    output <- ddply(data, by, function(x) mean.input(data=x, variable=variable, config=config))
  }
  if (export)  {
    write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
  }
  class(output) <- c("intsvy.mean", "data.frame")
  return(output)
  
}

