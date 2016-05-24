# -----------------------------------------------
## Determine distribution quantiles
# -----------------------------------------------
ntiles <- function( data, dv, 
                    factors = NaN,
                    bins = 5,
                    res.labels = FALSE
                    ) {
  
  ## weird bin argument
  if ( is.numeric(bins) == FALSE )
  {
    print("Weird bin argument...don't know what to do.")
    stop(call. = TRUE)
  }
  
  ##############################################################################################
  
  tempdata <- data
  tempdata$tempvarforntiles <- NaN                         # add temporary variable
  if (is.character(dv)) {
    which_column <- which(names(tempdata) == dv )
  } else if (is.numeric(dv)) {
    which_column <- dv 
  }
  
  
  # Check factor definition
  if (length(factors)<=1) {
    if (is.nan(factors))    {
      factors = 0*c(1:length(tempdata[[which_column]][]))
    } else {
      a <- factors
      factors <- list(tempdata[,a])[[1]]
    }
  } else if (TRUE) {
    a <- factors
    factors <- list(as.list(tempdata[,a]))[[1]]
  }


  # Rank-order and cut
  tempdata$tempranks <- unsplit(lapply(split(tempdata[[which_column]], factors), 
                                       rank, ties.method = "average"), factors)
 
  tempdata$tempvarforntiles <- unsplit(lapply(split(as.numeric(tempdata$tempranks), factors), 
                                              cut, bins, include.lowest = TRUE, right = FALSE, labels = res.labels), factors)
    
    
  #names(tempdata)[which(names(tempdata)=="tempvarforntiles")] <- res.name
  #tempdata$tempranks <- NULL
  
  # Output: Column of ranks
  tempdata$tempvarforntiles
}
