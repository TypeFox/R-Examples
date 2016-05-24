pit.fitBMAgamma <-
function(fit, ensembleData, dates=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 powfun <- function(x,power) x^power

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,dates=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates = NULL)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 PIT <- numeric(nObs)
 names(PIT) <- dataObsLabels(ensembleData)

 startup <- dataStartupSpeed(ensembleData)
 if (is.null(startup) & !is.na(fit$startup)) {
   if (length(fit$startup) != 1) stop("problem with startup specification")
   startup <- rep(fit$startup, nrow(ensembleData))
 }
 if (is.null(startup)) startup <- rep(controlBMAgamma()$startupSpeed,
                                      nrow(ensembleData))

 if (any(is.na(startup))) {
   if (is.null(controlBMAgamma()$startupSpeed)) 
     stop("default anemometer startup speed not specified")
   startup[is.na(startup)] <- controlBMAgamma()$startupSpeed
 }

 obs <- sapply( dataVerifObs(ensembleData), powfun, power = fit$power)
 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS)))  {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (fit$varCoefs[1] + fit$varCoefs[2]*f)^2
        
       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       PIT[i] <- cdfBMAgamma( obs[i], WEIGHTS = W, 
                          MEAN = MEAN[!M], VAR = VAR[!M]) 

       if (startup[i] > 0 && obs[i] < startup[i]) PIT[i] <- runif(0, max = PIT[i])
    }

 }

 PIT
}

