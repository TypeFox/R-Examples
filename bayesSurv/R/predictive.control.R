####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       predictive.control             ####
####                                            ####
#### FUNCTIONS:  predictive.control             ####
####################################################

### ======================================
### predictive.control
### ======================================
predictive.control <- function(predict, store, only.aver, quantile)
{
  
  if(length(predict) == 0) inpredict <- "arnost"
  else                     inpredict <- names(predict)
  tmp <- match("Et", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$Et <- FALSE
  if (!is.logical(predict$Et)) predict$Et <- FALSE
  tmp <- match("t", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$t <- FALSE
  if (!is.logical(predict$t)) predict$t <- FALSE
  tmp <- match("Surv", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$Surv <- FALSE
  if (!is.logical(predict$Surv)) predict$Surv <- FALSE
  tmp <- match("hazard", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$hazard <- FALSE
  if (!is.logical(predict$hazard)) predict$hazard <- FALSE
  tmp <- match("cum.hazard", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$cum.hazard <- FALSE
  if (!is.logical(predict$cum.hazard)) predict$cum.hazard <- FALSE
  if (!(predict$Et || predict$t || predict$Surv || predict$hazard || predict$cum.hazard))
    stop("Nothing to be predicted.")

  if(length(store) == 0) instore <- "arnost"
  else                   instore <- names(store)
  tmp <- match("Et", instore, nomatch=NA)
  if(is.na(tmp)) store$Et <- FALSE
  if (!is.logical(store$Et)) store$Et <- FALSE
  tmp <- match("t", instore, nomatch=NA)
  if(is.na(tmp)) store$t <- FALSE
  if (!is.logical(store$t)) store$t <- FALSE
  tmp <- match("Surv", instore, nomatch=NA)
  if(is.na(tmp)) store$Surv <- FALSE
  if (!is.logical(store$Surv)) store$Surv <- FALSE
  tmp <- match("hazard", instore, nomatch=NA)
  if(is.na(tmp)) store$hazard <- FALSE
  if (!is.logical(store$hazard)) store$hazard <- FALSE
  tmp <- match("cum.hazard", instore, nomatch=NA)
  if(is.na(tmp)) store$cum.hazard <- FALSE
  if (!is.logical(store$cum.hazard)) store$cum.hazard <- FALSE

  if (!predict$Et) store$Et <- FALSE
  if (!predict$t) store$t <- FALSE
  if (!predict$Surv) store$Surv <- FALSE
  if (!predict$hazard) store$hazard <- FALSE
  if (!predict$cum.hazard) store$cum.hazard <- FALSE

  if (length(quantile)){ if (sum(quantile < 0 | quantile > 1)) stop("Quantiles must lie between 0 and 1.") }
  else                 {only.aver <- TRUE }  
  if (missing(only.aver)) only.aver <- TRUE
  
  back <- list(predict = predict, store = store, only.aver=only.aver, quantile=quantile)
  return(back)
}

