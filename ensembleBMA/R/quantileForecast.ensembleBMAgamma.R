quantileForecast.ensembleBMAgamma <-
function(fit, ensembleData, quantiles=0.5, dates=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)
 
 weps <- 1.e-4

 matchITandFH(fit,ensembleData)

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,observations=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 dates <- modelDates(fit)

 Dates <- ensembleValidDates(ensembleData)

 nObs <- nrow(ensembleData)

 nForecasts <- ensembleSize(ensembleData)

 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(dataObsLabels(ensembleData),as.character(quantiles))

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (d in dates) {
    
    l <- l + 1

    WEIGHTS <- fit$weights[,d]
    if (all(Wmiss <- is.na(WEIGHTS))) next

    I <- which(as.logical(match(Dates, d, nomatch = 0)))
    
    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (fit$varCoefs[1,d] + fit$varCoefs[2,d]*f)^2

       fTrans <- sapply( f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans)*fit$biasCoefs[,d], 2, sum)
       MEAN[MEAN < 0] <- 0

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       Q[i,] <- sapply(quantiles, quantBMAgamma, WEIGHTS=W,
                        MEAN=MEAN[!M], VAR=VAR[!M])
    }
 }

  apply(Q, 2, powinv, power = fit$power)
}
