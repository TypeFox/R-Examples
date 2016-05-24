# -----------------------------------------------
## Outlier Screening 
# -----------------------------------------------
outlier <- function( data, dv, 
                     todo = "na", res.name = "outlier",
                     upper.limit = NaN, lower.limit = NaN,
                     limit.exact = FALSE,
                     upper.z = NaN, lower.z = NaN,
                     z.exact = FALSE, factors = NaN,
                     z.keep = TRUE, z.name = "zscores",
                     vsj = FALSE,
                     print.summary = TRUE) {
  
  ##############################################################################################
  ## Error Messages
  ## lower.limit >= upper.limit
  if ( !is.nan(upper.limit) & !is.nan(lower.limit) & lower.limit >= upper.limit )
  {
    print("Value for lower.limit must be smaller than value for upper.limit!")  
    stop(call. = TRUE)
  }
  
  ## Weird todo argument
  if ( todo != "na" & todo !="elim" & todo != "nothing")
  {
    print("Unknown specification of parameter todo!")
    stop(call. = TRUE)
  }
  
  ##############################################################################################
  
  tempdata <- data
  tempdata$tempvar <- 0                         # add new variable
  if (is.character(dv)) {
    which_column <- which(names(tempdata) == dv )
  } else if (is.numeric(dv)) {
    which_column <- dv 
  }
  
  ##############################################################################################
  ## Work on absolute limits
  ## upper limit
  if (!is.nan(upper.limit))
  {
    if (limit.exact == FALSE) {
      tempdata$tempvar[tempdata[which_column]>upper.limit] <- 1 
      #### dev note: Might be changed to tempdata[,which_clumn] 
	  #### in future releases (never change a running system though)
    } else {
      tempdata$tempvar[tempdata[which_column]>upper.limit] <- 1 
    }
  }

  ## lower Limit
  if (!is.nan(lower.limit))
  {
    if (limit.exact == FALSE) {
      tempdata$tempvar[tempdata[which_column]<lower.limit] <- 1 
    } else {
      tempdata$tempvar[tempdata[which_column]<=lower.limit] <- 1
    }
  }

  ##############################################################################################
  ## Work on z-scores
  ## z scores  <---- dev note: could be a starting point for a revised version of zscores.R.
  if ((!is.nan(upper.z)) | (!is.nan(upper.z)))
  {
    tempdata$zscores <- NaN 
    ### dev note: Could cause problems if the variable already exists (also applies to tempvar).
    ### might be changed to a more random name in future releases
    
    zindex <- as.logical(1-tempdata$tempvar)
    # Check factor definition
    if (length(factors)<=1) {
      if (is.nan(factors))    {
        factors = 0*c(1:length(tempdata[[which_column]][zindex]))
      } else {
        a <- factors
        factors <- list(tempdata[zindex,a])[[1]]
      }
    } else if (TRUE) {
      a <- factors
      factors <- list(as.list(tempdata[zindex,a]))[[1]]
    }
    # Compute z-scores
    tempdata$zscores[zindex] <- unsplit(lapply(split(tempdata[[which_column]][zindex], factors), scale), factors)
    
    ## upper z
    if (!is.nan(upper.z))   { 
      if (z.exact==FALSE) {
		    tempdata$tempvar[tempdata$zscores>upper.z] <- 1
      } else {
		  tempdata$tempvar[tempdata$zscores>=upper.z] <- 1
    }}
    ## lower z
    if (!is.nan(lower.z))   { 
      if (z.exact==FALSE) {
		  tempdata$tempvar[tempdata$zscores<lower.z] <- 1 
	  } else {
	  	tempdata$tempvar[tempdata$zscores<=lower.z] <- 1 
    }}
    
    ## Correct for cells with n = 1 (z = NA)
    tempdata$tempvar[is.na(tempdata$zscores)] <- 1
  }
  
  ##############################################################################################
  ## Brief summary of screening
  if (print.summary == TRUE)
  {
    summary_table <- table(tempdata$tempvar)
    n_trials <- length(tempdata$tempvar)
    if (length(summary_table) == 2) {
      n_outlier <- summary_table[[2]]
    } else {
      if (names(summary_table) == 0) {
        n_outlier <- 0
      } else if (names(summary_table) == 1) {
        n_outlier <- length(tempdata$tempvar) 
      }
    }
    percent <- round((n_outlier / n_trials) * 100,2)
    table_text <- paste("Function outlier() has screened ",n_trials," trial(s) in total. Of these were ",
                      n_outlier," trial(s) identified as outliers (",percent,"%).",sep="")
    print(table_text)
  }
  
  ##############################################################################################
  ## What to do with outliers? -> Create new dataframe
  if (todo == "na")               # replace with na
  {
     tempdata[[which_column]][tempdata$tempvar == 1] <- NA
  } else if (todo == "elim")      # eliminate cases at all
  {
    tempdata <- tempdata[tempdata$tempvar == 0 ,] 
  }
    
  ##############################################################################################
  ## give variable new name
  names(tempdata)[which(names(tempdata)=="tempvar")] <- res.name
  if ((!is.nan(upper.z)) | (!is.nan(upper.z))) {
    names(tempdata)[which(names(tempdata)=="zscores")] <- z.name
  }
  
  ##############################################################################################
  ## end
  tempdata
}
