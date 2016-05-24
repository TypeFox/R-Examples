pit.ensembleBMAgamma0 <-
function(fit, ensembleData, dates = NULL, randomizeATzero = FALSE, ...) 
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
     
       VAR <- fit$varCoefs[1,d] + fit$varCoefs[2,d]*f
        
       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs[,,d], 2, sum)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs[,,d],
                              2,sum), inverseLogit)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

      if (obs[i] == 0 && randomizeATzero) {
         PIT[i] <- runif(1,min=0,max=sum(W*PROB0[!M]))
       }
       else {
         PIT[i] <- cdfBMAgamma0( obs[i], WEIGHTS = W,
                          MEAN = MEAN[!M], VAR = VAR[!M], PROB0 = PROB0[!M])
       }

    }

 }

 PIT
}

