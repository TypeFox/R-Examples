cdf.ensembleBMAgamma0 <-
function(fit, ensembleData, values, dates = NULL, ...) 
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

 M <- !dataNA(ensembleData,observations=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=dates)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 if (is.null(dates)) dates <- modelDates(fit)

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

       CDF[i,] <- sapply( sapply( values, powfun, power = fit$power), 
                          cdfBMAgamma0, WEIGHTS = W, 
                          MEAN = MEAN[!M], VAR = VAR[!M], PROB0 = PROB0[!M]) 
    }

 }

 CDF
}

