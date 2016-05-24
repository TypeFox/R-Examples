ensembleBMAgamma <-
function(ensembleData, trainingDays, dates = NULL, 
         control = controlBMAgamma(), exchangeable = NULL)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 if (missing(trainingDays)) stop("trainingDays must be specified")

 call <- match.call()

 warmStart <- FALSE

 if (is.list(trainingDays)) trainingsDays <- trainingDays[[1]]

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

 biasCoefs <- matrix( NA, nrow = 2, ncol = nDates,
                      dimnames = list(NULL, dates))
 varCoefs <- array( NA, c(2, nDates), dimnames = list(NULL, dates))
 weights <- array( NA, c(nForecasts, nDates),
                      dimnames = list(ensMemNames, dates))

 trainTable <- rep(0, nDates)
 names(trainTable) <- dates

 nIter <- loglikelihood <- rep(0, nDates)
 names(nIter) <- names(loglikelihood) <- dates

 K <- 1:nForecasts

 L <- length(juliandates)
 twin <- 1:trainingDays

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
        fit <- fitBMAgamma(ensembleData[D,-K[kNA]], control = control,
                            exchangeable = x)
      }
      else {
        fit <- fitBMAgamma(ensembleData[D,], control = control,
                            exchangeable = exchangeable)
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

   biasCoefs[,i] <- fit$biasCoefs
   varCoefs[,i] <- fit$varCoefs
   weights[K[!kNA],i] <- fit$weights

 }

 structure(list(training = list(days=trainingDays,lag=lag,table=trainTable),
                biasCoefs = biasCoefs, 
                varCoefs = varCoefs, weights = weights, nIter = nIter,
                exchangeable = exchangeable, power = fit$power),
                forecastHour = forecastHour, 
                initializationTime = ensembleItime(ensembleData),
                call = match.call(), class = c("ensembleBMAgamma", "ensembleBMA"))
}

