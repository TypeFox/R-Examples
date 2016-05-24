quantileForecast.fitBMAgamma <-
function(fit, ensembleData, quantiles = 0.5, dates=NULL, ...) 
{ 
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)
 
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,observations=FALSE,dates=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]
 nObs <- nrow(ensembleData)

 Q <- matrix(NA, nObs, length(quantiles))
 dimnames(Q) <- list(dataObsLabels(ensembleData),as.character(quantiles))

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights
 if (!all(Wmiss <- is.na(WEIGHTS))) {

   for (i in 1:nObs) {
    
       f <- ensembleData[i,]       

       M <- is.na(f) | Wmiss

       VAR <- (fit$varCoefs[1] + fit$varCoefs[2]*f)^2

       fTrans <- sapply( f, power, fit$power)

       MEAN <- apply(rbind(1, fTrans)*fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       Q[i,] <- sapply( quantiles, quantBMAgamma, WEIGHTS=W,
                        MEAN=MEAN[!M], VAR=VAR[!M])
   }
 }

  apply(Q, 2, powinv, power = fit$power)
}

