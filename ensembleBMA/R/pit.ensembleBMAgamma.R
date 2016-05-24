pit.ensembleBMAgamma <-
function(fit, ensembleData, dates = NULL, ...) 
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

 M <- !dataNA(ensembleData)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 dates <- modelDates(fit)

 Dates <- ensembleValidDates(ensembleData)

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 nForecasts <- ensembleSize(ensembleData)
 
 PIT <- numeric(nObs)
 names(PIT) <- dataObsLabels(ensembleData)

 startup <- dataStartupSpeed(ensembleData)
 if (is.null(startup) & !is.null(fit$startup)) {
   if (length(fit$startup) != 1) stop("problem with startup specification")
   startup <- rep(fit$startup, nrow(ensembleData))
 }

 if (is.null(startup)) startup <- controlBMAgamma()$startupSpeed

 if (length(startup == 1)) startup <- rep(startup, length = nObs)

 if (any(is.na(startup))) {
   if (is.null(controlBMAgamma()$startupSpeed)) 
     stop("default anemometer startup speed not specified")
   startup[is.na(startup)] <- controlBMAgamma()$startupSpeed
 }

 obs <- sapply( dataVerifObs(ensembleData), powfun, power = fit$power)
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
        
       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs[,d], 2, sum)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       PIT[i] <- cdfBMAgamma( obs[i], WEIGHTS = W, 
                              MEAN = MEAN[!M], VAR = VAR[!M]) 

       if (startup[i] > 0 && obs[i] <= startup[i]) PIT[i] <- runif(0, max = PIT[i])  
    }

 }

 PIT
}

