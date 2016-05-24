ensembleBMAnormal <-
function(ensembleData, trainingDays, dates = NULL, 
         control = controlBMAnormal(), exchangeable = NULL, minCRPS = FALSE)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 if (missing(trainingDays)) stop("trainingDays must be specified")

 call <- match.call()

 warmStart <- FALSE

 if (is.list(trainingDays)) trainingsDays <- trainingsDays[[1]]

 ensMemNames <- ensembleMembers(ensembleData)
 nForecasts <- length(ensMemNames)

 exchangeable <- getExchangeable( exchangeable, ensembleGroups(ensembleData),
                                  nForecasts)

# remove instances missing all forecasts, obs or dates

 M <- !dataNA(ensembleData)
 if (!all(M)) ensembleData <- ensembleData[M,]
 
 nObs <- nrow(ensembleData)
 if (!nObs) stop("no data")

 Dates <- as.character(ensembleValidDates(ensembleData))
 DATES <- sort(unique(Dates))

 julianDATES <- ymdhTOjul(DATES)
 incr <- min(1,min(diff(julianDATES))) ## incr may be fractional for hours

 forecastHour <- ensembleFhour(ensembleData)
 lag <- ceiling( forecastHour / 24 )

## dates that can be modeled by the training data (ignoring gaps)

 dates <- getDates( DATES, julianDATES, dates, trainingDays, lag, incr)
 juliandates <- ymdhTOjul(dates)
 nDates <- length(dates)

 biasCoefs <- array( NA, c(2, nForecasts, nDates),
                     dimnames = list(NULL, ensMemNames, dates))

 if (control$equalVariance) {
   sd <- rep(NA, nDates) 
   names(sd) <- dates
 }
 else {
   sd <- array(NA,c(nForecasts,nDates), 
               dimnames=list(ensMemNames, dates))
 }

 weights <- array( NA, c(nForecasts, nDates))
 dimnames(weights) <- list(ensMemNames, dates)

 trainTable <- rep(0, nDates)
 names(trainTable) <- dates

 nIter <- loglikelihood <- rep(0, nDates)
 names(nIter) <- names(loglikelihood) <- dates

 L <- length(juliandates)
 twin <- 1:trainingDays

## temp <- data.frame(julian = julianDATES,date = DATES)
## print(temp)

 K <- 1:nForecasts

 cat("\n")

 l <- 0
 for(i in seq(along = juliandates)) {

    I <- (juliandates[i]-lag*incr) >= julianDATES
    if (!any(I)) stop("insufficient training data")

    j <- which(I)[sum(I)]

    if (j != l) {

      twin <- (j+1) - (1:trainingDays)
      D <- as.logical(match(Dates, DATES[twin], nomatch=0))
      if (!any(D)) stop("this should not happen")
      d <- ensembleValidDates(ensembleData[D,])
      if (length(unique(d)) != trainingDays) stop("wrong # of training days")
      cat("modeling for date", dates[i], "...")

      kNA <- apply(ensembleForecasts(ensembleData[D,]), 2, 
                   function(x) all(is.na(x)))

      if (any(kNA)) {
        if (!is.null(x <- exchangeable)) x <- exchangeable[-K[kNA]]
        fit <- fitBMAnormal(ensembleData[D,-K[kNA]], control = control,
                            exchangeable = x)
      }
      else {
        fit <- fitBMAnormal(ensembleData[D,], control = control,
                            exchangeable = exchangeable)
      }
  
      if (minCRPS) {

        CRPSobjective <- function(SD) {
            if (any(kNA)) {
              crpsNormal(SD, fit$weights, fit$biasCoefs, ensembleData[D,-K[kNA]])
            }
            else {
              crpsNormal(SD, fit$weights, fit$biasCoefs, ensembleData[D,])
            }
         }

        if (control$equalVariance) {
          fit$sd <- optimize(CRPSobjective, interval = c(0, 6*fit$sd))$minimum
        }
        else {
          opt <- optim(fit$sd, CRPSobjective, method = "BFGS")
          if (!opt$convergence) {
            fit$sd <- opt$par
          }
          else {
            warning("CRPS could not be minimized")
          }
        }
      }
       l <- j ## last model fit
       trainTable[i] <- length(unique(Dates[D]))
       nIter[i] <- fit$nIter
       loglikelihood[i] <- fit$loglikelihood
       if (warmStart) control$start$weights <- weights[,i]
       cat("\n")
     }
   else {
     trainTable[i] <- -abs(trainTable[i-1])
     nIter[i] <- -abs(nIter[i-1])
     loglikelihood[i] <- loglikelihood[i-1]
   }

   biasCoefs[,K[!kNA],i] <- fit$biasCoefs
   weights[K[!kNA],i] <- fit$weights

   if (control$equalVariance) {
     sd[i] <- fit$sd 
   }
   else {
     sd[K[!kNA],i] <- as.vector(fit$sd)
   }

 }

 structure(list(training = list(days=trainingDays,lag= lag,table=trainTable),
                biasCoefs = biasCoefs, sd = sd, weights = weights,
                nIter = nIter, exchangeable = exchangeable),
                forecastHour = forecastHour, 
                initializationTime = ensembleItime(ensembleData),
                call = match.call(), 
                class = c("ensembleBMAnormal", "ensembleBMA"))
}

