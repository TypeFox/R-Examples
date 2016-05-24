cdf.fitBMAnormal <-
function(fit, ensembleData, values, dates=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,observations=FALSE,dates=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=NULL)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 nObs <- nrow(ensembleData)
 if (!nObs) stop("no data")

 CDF <- matrix( NA, nrow = nObs, ncol = length(values))
 dimnames(CDF) <- list(dataObsLabels(ensembleData), as.character(values)) 

 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights
     
 if (!all(Wmiss <- is.na(WEIGHTS))) {
     
    SD <- if (!is.null(dim(fit$sd))) {
            fit$sd
          }
          else rep(fit$sd, nForecasts)

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs, 2, sum)

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

