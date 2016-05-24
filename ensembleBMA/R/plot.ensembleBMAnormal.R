plot.ensembleBMAnormal <-
function(x, ensembleData, dates = NULL, ask = TRUE, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 par(ask = ask)

 weps <- 1.e-4

 matchITandFH(x,ensembleData)

 exchangeable <- x$exchangeable

 ensembleData <- ensembleData[,matchEnsembleMembers(x,ensembleData)]

 M <- !dataNA(ensembleData)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(x)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) x <- x[fitDates[M$fit]]

 dates <- modelDates(x)

 Dates <- ensembleValidDates(ensembleData)

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 if (nObs == 0) obs <- rep( NA, nrow(ensembleData))

 nForecasts <- ensembleSize(ensembleData)

 obs <- dataVerifObs(ensembleData)
 lat <- ensembleData$latitude
 lon <- ensembleData$longitude

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0 
 for (d in dates) {

    l <- l + 1

    WEIGHTS <- x$weights[,d]
     
    if (all(Wmiss <- is.na(WEIGHTS))) next
     
    SD <- if (!is.null(dim(x$sd))) {
            x$sd[,d] 
          }
          else rep(x$sd[d], nForecasts)

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]
     
       MEAN <- apply(rbind(1, f) * x$biasCoefs[,,d], 2, sum)

       M <- is.na(f) | Wmiss

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

