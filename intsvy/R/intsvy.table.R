# weights = c("BRR / JK")
#intsvy.table <- function(variable, by, data, final_weight="W_FSTUWT", brr_weight="W_FSTR") {
  
intsvy.table <- function(variable, by, data, config) {
  table.input <- function(variable, data, config) {
    # data is empty
    if (sum(is.na((data[[variable]])))==length(data[[variable]])) {
      result <- data.frame(NA, "Freq"=0, "Percentage"=NA, "Std.err."= NA)  
      names(result)[1] <- variable 
      return(result)
    }

    # BRR / JK
    if (config$parameters$weights == "JK") {
      # jack knife
      # in PIRLS / TIMSS
      
      # Replicate weighted %s (sampling error)
      tabrp <- as.matrix(sapply(1:max(data[[config$variables$jackknifeZone]]), function(x) 
        percent(as.factor(as.numeric(data[[variable]])), 
                total=FALSE, 
                weights= ifelse(data[[config$variables$jackknifeZone]] == x, 
                                2*data[[config$variables$weight]]*data[[config$variables$jackknifeRep]], 
                                data[[config$variables$weight]]), na.rm = TRUE)))
      
      # Total weighted %                                                                      
      tabtot <- percent(as.factor(as.numeric(data[[variable]])), 
                        weights= data[[config$variables$weight]], na.rm = TRUE, total=FALSE)
      # Standard error
      if (length(tabtot)!=1) {
        tabse <- apply((tabrp-tabtot)^2, 1, sum)^(1/2)
      } else {
        tabse <-0
      }
      
    } else {
      # balanced repeated replication
      # Replicate weighted %s (sampling error)
      # in PISA / PIAAC
      tabrp <- as.matrix(sapply(1:config$parameters$BRRreps, function(i) 
        percent(as.factor(as.numeric(data[[variable]])), total=FALSE, 
                weights=  data[[paste(config$variables$weightBRR, i , sep="")]], na.rm=TRUE)))     
      
      # Total weighted %                                                                      
      tabtot <- percent(as.factor(as.numeric(data[[variable]])), weights= data[[config$variables$weightFinal]], na.rm = TRUE, total=FALSE)
      # Standard error
      if (length(tabtot)!=1) {
        tabse <- sqrt(rowSums((tabrp-tabtot)^2) / 20)
      } else {
        tabse <- 0
      }

    } 
    
    result <- data.frame(table(data[[variable]][drop=TRUE]), 
                         "Percentage" = round(as.numeric(tabtot), 2), 
                         "Std.err."= round(tabse, 2))
    
    names(result)[1] <- variable # var label for table, otherwise prints "Var1"
    return(result)
  }
  
  # If by not supplied, calculate for complete sample    
  if (missing(by)) { 
    output <- table.input(variable=variable, data=data, config=config)
  } else {
    # Convert by variables to characters for ddply application
    for (i in by) {
      data[[c(i)]] <- as.character(data[[c(i)]])
    }
    output <- ddply(data, by, function(x) table.input(data=x, variable=variable, config=config))
  }
  class(output) <- c("intsvy.table", "data.frame")
  output
}

