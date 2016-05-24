quantileForecast.ensembleBMAnormal <-
function(fit, ensembleData, quantiles = 0.5, dates = NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 weps <- 1.e-4

 matchITandFH(fit,ensembleData)

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData, observations = FALSE)
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

    SD <- if (!is.null(dim(fit$sd))) fit$sd[,d] else rep(fit$sd[d], nForecasts)

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       MEAN <- apply(rbind(1, f)*fit$biasCoefs[,,d], 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W[!M] <- W[!M]/sum(W[!M])
       }

       Q[i,] <- sapply(quantiles, quantBMAnormal,
                       WEIGHTS=W[!M], MEAN=MEAN[!M], SD=SD[!M])
    }
 }

 Q
}

