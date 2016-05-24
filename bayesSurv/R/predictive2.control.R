#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2004)                              ####
####                                                 ####
#### FILE:       predictive2.control.R               ####
####                                                 ####
#### FUNCTIONS:  predictive2.control                 ####
#########################################################

### ======================================
### predictive2.control
### ======================================
predictive2.control <- function(predict, only.aver, quantile, obs.dim, time0, Gspline, n)
{
  
  ## predict
  ## ========
  if(length(predict) == 0) inpredict <- "arnost"
  else                     inpredict <- names(predict)

  tmp <- match("density", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$density <- FALSE
  if (!is.logical(predict$density)) predict$density <- FALSE
  
  tmp <- match("Surv", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$Surv <- FALSE
  if (!is.logical(predict$Surv)) predict$Surv <- FALSE
  
  tmp <- match("hazard", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$hazard <- FALSE
  if (!is.logical(predict$hazard)) predict$hazard <- FALSE
  
  tmp <- match("cum.hazard", inpredict, nomatch=NA)
  if(is.na(tmp)) predict$cum.hazard <- FALSE
  if (!is.logical(predict$cum.hazard)) predict$cum.hazard <- FALSE
  
  if (!(predict$density || predict$Surv || predict$hazard || predict$cum.hazard))
    stop("Nothing to be predicted.")

  ## quantile
  ## =========
  if (length(quantile)){ if (sum(quantile < 0 | quantile > 1)) stop("Quantiles must lie between 0 and 1.") }
  else                 {only.aver <- TRUE }

  ## only.aver
  ## =========
  if (missing(only.aver)) only.aver <- TRUE

  ## Gspline
  ## ========
  if(length(Gspline) == 0) inGspline <- "arnost"
  else                     inGspline <- names(Gspline)
  tmp <- match("dim", inGspline, nomatch=NA)
  if(is.na(tmp)) stop("Gspline$dim must be given")
  Gspline$dim <- Gspline$dim[1]
  if (is.na(Gspline$dim)) stop("Gspline$dim must be given")
  if (Gspline$dim < 1 | Gspline$dim > 2) stop("Gspline$dim must be either 1 or 2")

  tmp <- match("K", inGspline, nomatch=NA)
  if(is.na(tmp)) stop("Gspline$K must be given")
  if (length(Gspline$K) < Gspline$dim) stop("Incorrect Gspline$K supplied")
  Gspline$K <- Gspline$K[1:Gspline$dim]
  if (sum(is.na(Gspline$K))) stop("Incorrect Gspline$K supplied")
  if (sum(Gspline$K < 0)) stop("Incorrect Gspline$K supplied")

  Gspline$total.length <- prod(2*Gspline$K + 1)

  ## obs.dim
  ## ========
  if (Gspline$dim >= 2){
     if (!length(obs.dim)) stop("obs.dim must be given")
     if (length(obs.dim) != n) stop("obs.dim has a different length than the rest of data")
     if (sum(is.na(obs.dim))) stop("obs.dim must not be contain NA's")
     if (sum(obs.dim <= 0 | obs.dim > Gspline$dim)) stop("obs.dim must contain only values from 1,...,Gspline$dim")
   }
   else
     obs.dim <- rep(1, n)
  
  ## time0
  ## ======
  if (!length(time0)) time0 <- rep(0, Gspline$dim)
  if (length(time0) == 1) time0 <- rep(time0, Gspline$dim)
  if (length(time0) < Gspline$dim) stop("Incorrect time0 parameter supplied")
  time0 <- time0[1:Gspline$dim]
  if (sum(is.na(time0))) stop("Incorrect time0 parameter supplied")
  if (sum(time0 < 0)) stop("time0 must be non-negative")

  back <- list(predict = predict, only.aver=only.aver, quantile = quantile, obs.dim=obs.dim, time0=time0, Gspline=Gspline)
  return(back)
}
