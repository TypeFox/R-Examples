plot.fitBMAnormal <-
function(x, ensembleData, dates=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 exchangeable <- x$exchangeable

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 ensembleData <- ensembleData[,matchEnsembleMembers(x,ensembleData)]

 M <- !dataNA(ensembleData,dates=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(x)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=NULL)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) x <- x[fitDates[M$fit]]

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 nForecasts <- ensembleSize(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- x$weights
     
 if (!all(Wmiss <- is.na(WEIGHTS))) {
     
    SD <- if (!is.null(dim(x$sd))) {
            x$sd
          }
          else rep(x$sd, nForecasts)

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       MEAN <- apply(rbind(1, f) * x$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       plotBMAnormal( WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M],
                      obs = obs[i], exchangeable = exchangeable)

    }

 }

 invisible()
}

