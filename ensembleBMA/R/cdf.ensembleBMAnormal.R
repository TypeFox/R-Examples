cdf.ensembleBMAnormal <-
function(fit, ensembleData, values, dates = NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 weps <- 1.e-4

 matchITandFH(fit,ensembleData)

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,observations=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=dates)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 if (is.null(dates)) dates <- modelDates(fit)

 Dates <- ensembleValidDates(ensembleData)

 nObs <- nrow(ensembleData)
 if (!nObs) stop("no data")

 Dates <- ensembleValidDates(ensembleData)

 nForecasts <- ensembleSize(ensembleData)

 CDF <- matrix( NA, nrow = nObs, ncol = length(values))
 dimnames(CDF) <- list(dataObsLabels(ensembleData), as.character(values)) 

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0 
 for (d in dates) {

    l <- l + 1

    WEIGHTS <- fit$weights[,d]
     
    if (all(Wmiss <- is.na(WEIGHTS))) next
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd[,d] 
          }
          else rep(fit$sd[d], nForecasts)

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,d], 2, sum)

       M <- is.na(f) | Wmiss

       W <- WEIGHTS

       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       CDF[i,] <- sapply( values, cdfBMAnormal,
                          WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M]) 

    }

 }

 CDF
}

