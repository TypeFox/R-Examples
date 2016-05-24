stpRule <-
function(method = c("fixed", "precision", "info", "class"),
         it_min, it_left,
         it_crit, se_crit, info_crit,
         it_obs,  se_obs,  info_obs,                        
         categ_est )
{
  
  stp  <- 0   # default to continue (stp <- 1 means stop)
  term <- NA  # default to no termination

# - If the number of items is less than the minimum: continue.
# - If all of the items are used up: stop.
  if( (it_obs < it_min) | (it_left == 0) ){
    if( it_left == 0 ){
      stp <- 1
      term <- "total"
    }
  } else{

# 1) If we have reached the maximum number of items, stop:
    if( any(method == "fixed") ){
      if( !is.null( it_crit ) ){
        if( it_obs >= it_crit ){
          stp <- 1
          term <- "fixed"
        } # END if STATEMENT
      }
    } # END if 1) STATEMENT
  
# 2) If we are using precision termination, check threshold or change:
    if( any(method == "precision") ){
      if( !is.na(se_obs[it_obs]) & !is.null(se_crit$crit) ){
        if( (se_crit$method == "threshold") & (se_obs[it_obs] <= se_crit$crit) ){
          stp <- 1
          term <- "se_thresh"
        } else if( (se_crit$method == "change") & (it_obs >= 2) & (!is.na(se_obs[it_obs - 1])) ){
          if( abs(diff(se_obs)[it_obs - 1]) <= se_crit$crit ){
        	    stp <- 1
        	    term <- "se_change"
        	  } # END if STATEMENT
        	} # END ifelse STATEMENTS
      }
    } # END if 2) STATEMENT    
# If threshold --> the current SEM needs to be less than the crit.
# If change --> the abs dif between the last two SEMs need to be less than the crit.

    
# 3) If we are using info termination, check threshold or change:
    if( any(method == "info") ){
      if( !is.na(info_obs[it_obs]) & !is.null(info_crit$crit) ){
        if( (info_crit$method == "threshold") & (info_obs[it_obs] >= info_crit$crit) ){
          stp <- 1
          term <- "fi_thresh"
        } else if( (info_crit$method == "change") & (it_obs >= 2) & (!is.na(info_obs[it_obs - 1])) ){
        	  if( abs(diff(info_obs)[it_obs - 1]) <= info_crit$crit ){
            stp <- 1
            term <- "fi_change"
          } # END if STATEMENT
        } # END ifelse STATEMENTS
      }
    } # END if 2) STATEMENT
# If threshold --> the current info needs to be greater than the crit.
# If change --> the abs dif between the last two infos need to be less than the crit.

  
# 4) If we have decided on a class, stop:
    if( any( method == "class" ) ){
      if( !is.na(categ_est) ){
        stp <- 1
        term <- "class"
      }
    } # END if 3) STATEMENT

  } # END ifelse STATEMENTS
  
  return( list(stp = stp, term = term) )
  
} # END stpRule FUNCTION
