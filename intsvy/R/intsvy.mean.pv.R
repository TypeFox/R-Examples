
intsvy.mean.pv <- 
function(pvnames, by, data, export=FALSE, name= "output", folder=getwd(), config) {
  pv.input <- function(pvnames, data, config) {
    # If there is only one observation print NA
    if (nrow(data) <= 1)  
      return(data.frame("Freq"= length(data[[config$variables$weightFinal]]), "Mean"= NA, "s.e."= NA, "SD"=NA, "s.e"=NA))

    # BRR / JK
    if (config$parameters$weights %in% c("BRR","mixed_piaac")) {
      # balanced repeated replication
      # Replicate weighted %s (sampling error)
      # in PISA / PIAAC
      
      # Replicate weighted sds and means of 5 PVs (sampling error)
      R.mean <- sapply(pvnames, function(k) 
        sapply(1:config$parameters$BRRreps, function(i) 
          weighted.mean(data[[k]], 
                        data[[paste0(config$variables$weightBRR, i)]], na.rm = TRUE)))
      
      R.sd <- sapply(pvnames, function(x) 
        sapply(1:config$parameters$BRRreps, function(i)
          (sum(data[[paste0(config$variables$weightBRR, i)]]*
                 (data[[x]]-R.mean[i, x])^2, na.rm = TRUE)/
             sum(data[[paste0(config$variables$weightBRR, i)]], na.rm = TRUE))^(1/2)))
      
      # Grand mean of 5 PVs (imputation variance)
      PV.mean <- sapply(pvnames, function(x) 
        weighted.mean(data[[x]], data[[config$variables$weightFinal]], na.rm = TRUE))
      
      PV.sd <- sapply(pvnames, function(x)
        (sum(data[[config$variables$weightFinal]]*(data[[x]]-PV.mean[x])^2, na.rm=TRUE)/
           sum(data[[config$variables$weightFinal]], na.rm = TRUE))^(1/2))
      
      # Mean of means (the one is reported)
      MEAN.m <- mean(PV.mean)
      SD.m <- mean(PV.sd)
      
      cc = 1/20
      if (config$parameters$weights == "mixed_piaac") {
        cntName <- as.character(unique(data[,config$variables$countryID]))[1]
        cc <- piaacReplicationScheme[cntName,"c"]
        if (is.na(cc)) cc <- 1
        if (length(unique(piaacReplicationScheme[as.character(unique(data[,config$variables$countryID])),"c"])) > 1) {
          warning(paste("In PIAAC study different replications schemes were applied in different countries. \n In the selected set of countries more than one scheme was used. \n Further estimation is performed with coefficient c =", cc))
        }
      }
      
      # Sampling variance; imputation variance; and SEs
      var.mean.w <- mean(sapply(seq_along(pvnames), function(i) cc*sum((R.mean[,i]-PV.mean[i])^2)))
      var.mean.b <- (1/(length(pvnames)-1))*sum(sapply(seq_along(pvnames), function(i) (PV.mean[i]-MEAN.m)^2))
      mean.se <-(var.mean.w+(1+1/length(pvnames))*var.mean.b)^(1/2)
      
      var.sd.w <- mean(sapply(seq_along(pvnames), function(i) cc*sum((R.sd[,i]-PV.sd[i])^2)))
      var.sd.b <- (1/(length(pvnames)-1))*sum(sapply(seq_along(pvnames), function(i) (PV.sd[i]-SD.m)^2))
      sd.se <-(var.sd.w+(1+1/length(pvnames))*var.sd.b)^(1/2)
      
      result <- data.frame("Freq"= length(data[[config$variables$weightFinal]]), "Mean"= mean(MEAN.m), "s.e."= mean.se, 
                           "SD"=mean(SD.m), "s.e"=sd.se)

    } 
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS

      # Replicate weights
      R.wt <- sapply(1:max(data[[config$variables$jackknifeZone]]), function(x) 
              ifelse(data[[config$variables$jackknifeZone]] == x, 
                             2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], data[[config$variables$weight]]))
      
      # Estimates of PV1 (for sampling error)
      R.mean1 <- sapply(1:ncol(R.wt), function(x) 
        weighted.mean(data[[pvnames[[1]]]], R.wt[, x], na.rm = TRUE))                                                                  
      
      R.sd1 <- sapply(1:ncol(R.wt), function (x) 
        (sum(R.wt[, x]*(data[[pvnames[[1]]]]-R.mean1[x])^2, na.rm = TRUE)/sum(R.wt[, x], na.rm = TRUE))^(1/2))
      
      # Grand mean of 5 PVs (for imputation variance)
      R.mean <- sapply(pvnames, function(x) 
        weighted.mean(data[[x]], data[[config$variables$weight]], na.rm = TRUE))
      
      R.sd <- sapply(1:5, function(x) 
        (sum(data[[config$variables$weight]]*(data[[pvnames[x]]]-R.mean[x])^2, na.rm=TRUE)/sum(data[[config$variables$weight]], na.rm = TRUE))^(1/2))
      
      # Sampling variance (1st PV); imputation variance; SEs
      v.meanw <- sum((R.mean1-R.mean[1])^2);    v.meanb <- (1+1/5)*var(R.mean)
      v.sdw <- sum((R.sd1 - R.sd[1])^2);  v.sdb <- (1+1/5)*var(R.sd)
      mean.se <-  (v.meanw+v.meanb)^(1/2); sd.se <- (v.sdw+v.sdb)^(1/2)
      
      result <- data.frame("Freq"= length(data[[config$variables$weight]]), "Mean"= mean(R.mean), "s.e."= mean.se, 
                           "SD"=mean(R.sd), "s.e"=sd.se)
      
    }

    return(round(result, 2))
    
  }
  
  # If by no supplied, calculate for the complete sample    
  if (missing(by)) { 
    output <- pv.input(pvnames=pvnames, data=data, config=config)
  }  else {
    for (i in by) {
      data[[c(i)]] <- as.factor(data[[c(i)]])
    }
    output <- ddply(data, by, function(x) pv.input(data=x, pvnames=pvnames, config=config))
  }

  if (export)  {
    write.csv(output, file=file.path(folder, paste(name, ".csv", sep="")))
  }
  class(output) <- c("intsvy.mean", "data.frame")
  return(output)
  
}

