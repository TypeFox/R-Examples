pit.fitBMAgamma0 <-
function(fit, ensembleData, dates=NULL, randomizeATzero = FALSE, ...) 
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

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=NULL)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 PIT <- numeric(nObs)
 names(PIT) <- dataObsLabels(ensembleData)

 obs <- sapply( dataVerifObs(ensembleData), powfun, power = fit$power)
 nForecasts <- ensembleSize(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS)))  {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- fit$varCoefs[1] + fit$varCoefs[2]*f
        
       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)
       
       PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*fit$prob0coefs,
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

