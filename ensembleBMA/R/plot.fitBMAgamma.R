plot.fitBMAgamma <-
function(x, ensembleData, dates=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 exchangeable <- x$excchangeable
 
 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

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

 if (nObs == 0) obs <- rep( NA, nrow(ensembleData))

 nForecasts <- ensembleSize(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 obs <- powfun( obs, power = x$power)

 WEIGHTS <- x$weights

 if (!all(Wmiss <- is.na(WEIGHTS)))  {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (x$varCoefs[1] + x$varCoefs[2]*f)^2
        
       fTrans <- sapply(f, powfun, power = x$power)

       MEAN <- apply(rbind(1, fTrans) * x$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       plotBMAgamma( WEIGHTS = W, MEAN = MEAN[!M], VAR = VAR[!M],
                     obs = obs[i], exchangeable = exchangeable)

    }

 }

 invisible()
}

