# function to impute missing SD | Bracken's (1992) method using the 
# overall SD to mean ratio (coefficient of variation)
impute_SD_Bracken1992 <- function(aDataFrame, 
                                  columnSDnames, 
                                  columnXnames) {
                                  
  for(i in 1:length(columnSDnames)) {
    imputedDataColumn <- aDataFrame[columnSDnames[i]]
    missingness <- is.na(imputedDataColumn) 
    SD_to_mean_ratio <- sum(imputedDataColumn, na.rm = TRUE) / 
                          sum(aDataFrame[columnXnames[i]], na.rm = TRUE)
    X_of_missing_SD <- aDataFrame[[columnXnames[i]]][missingness]
    imputedDataColumn[missingness] <- X_of_missing_SD * SD_to_mean_ratio
    aDataFrame[columnSDnames[i]] <- imputedDataColumn
  }
  return(aDataFrame)
}

# function to impute missing SD | full-random hot deck imputation 
# imputation chosen at random, with replacement)
# m is the number of randomly imputed datasets (default 1)
impute_SD_HotDeck_fullRandom <- function(aDataFrame, 
                                         columnNames, 
                                         M = 1) {
                                         
  if(M == 1) {
    aList <- impute_SD_HotDeck_fullRandom_HELPER(aDataFrame, columnNames)
  } else {
    aList <- lapply(1:M, function(x) impute_SD_HotDeck_fullRandom_HELPER(aDataFrame, columnNames))
  }
  
  return(aList)
}

impute_SD_HotDeck_fullRandom_HELPER <- function(aDataFrame, 
                                                columnNames) {
                                                
  for(i in columnNames) {
    imputedDataColumn <- aDataFrame[i]
    missingness <- is.na(imputedDataColumn)
    imputedDataColumn[missingness] <- sample(aDataFrame[[i]][!missingness], 
                                             sum(missingness), 
                                             replace = TRUE)
    aDataFrame[i] <- imputedDataColumn
  }
  return(aDataFrame)
}


# function to impute missing SD | nearest neighbour hot deck imputation 
# (imputation chosen at random from a range neighbouring values, imputations 
# based on a sorted list of the paired means ) m is the number of randomly 
# imputed datasets (default 1)
impute_SD_HotDeck_nearestNeighbour <- function(aDataFrame, 
                                               columnSDnames, 
                                               columnXnames,
                                               M = 1,
                                               range = 3) {
  
  if(M == 1) {
    for(i in 1:length(columnSDnames)) 
      aDataFrame <- impute_SD_HotDeck_nearestNeighbour_HELPER(aDataFrame, 
                                                              columnSDnames[i], 
                                                              columnXnames[i], 
                                                              range);
    return(aDataFrame)
  }
  
  return( lapply(1:M, function(x) {
    for(i in 1:length(columnSDnames)) 
      aDataFrame <- impute_SD_HotDeck_nearestNeighbour_HELPER(aDataFrame, 
                                                              columnSDnames[i], 
                                                              columnXnames[i], 
                                                              range);
    return(aDataFrame);
  })
  )
  
}

impute_SD_HotDeck_nearestNeighbour_HELPER <- function(aDataFrame, 
                                                      columnSD_name, 
                                                      columnX_name, 
                                                      range = 3) {
                                                      
  aDataFrame <- cbind(aDataFrame, oldOrder = c(1:nrow(aDataFrame)))
  ordered <- aDataFrame[order(aDataFrame[columnX_name]), ]
  NA_index <- which(is.na(ordered[columnSD_name]))
  
  for(i in 1:(length(NA_index))) {
    lowerSub <- subset(ordered, ordered[columnX_name] < ordered[NA_index[i],columnX_name]); 
    lowerSub <- lowerSub[!is.na(lowerSub[columnSD_name]), ]
    upperSub <- subset(ordered, ordered[columnX_name] > ordered[NA_index[i],columnX_name]); 
    upperSub <- upperSub[!is.na(upperSub[columnSD_name]), ]
    if((range * 2) > (nrow(lowerSub) + nrow(upperSub))) 
      stop("Nearest neighbour hot deck: 'range' parameter is too large.")
    
    if ((nrow(lowerSub) - range) < 1) {
      neighbours <- rbind(lowerSub, upperSub[1:(range + abs(nrow(lowerSub) - range)),])
    } else if ((nrow(upperSub) - range) < 1) {
      neighbours <- rbind(lowerSub[(nrow(lowerSub) - range - abs(nrow(upperSub) - range) + 1):nrow(lowerSub), ], upperSub)  
    } else neighbours <- rbind(lowerSub[(nrow(lowerSub)  -range + 1):nrow(lowerSub), ], upperSub[1:range, ])
    
    aDataFrame[ordered[NA_index[i],"oldOrder"], columnSD_name] <- sample(neighbours[[columnSD_name]], 1, replace = TRUE)
  }
  
  return(aDataFrame[, -ncol(aDataFrame)])
}