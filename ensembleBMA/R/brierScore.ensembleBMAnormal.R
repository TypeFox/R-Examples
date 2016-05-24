brierScore.ensembleBMAnormal <-
function(fit, ensembleData, thresholds, dates = NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
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

 y <- dataVerifObs(ensembleData)
 nObs <- length(y)

 nForecasts <- ensembleSize(ensembleData) 

 ensembleData <- ensembleForecasts(ensembleData)

 MAT <-  t(outer(y, thresholds, "<="))

 bsClimatology <- apply(sweep(MAT, MARGIN = 1, FUN = "-", 
                        STATS = apply(MAT,1,mean))^2, 1, mean)
 
 bsVotingEns <- apply(ensembleData, 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean, na.rm = TRUE),
                 thresholds = thresholds)

 bsVoting <- apply((bsVotingEns - MAT)^2, 1, mean, na.rm = TRUE)

 MAT <- matrix( NA, nObs, length(thresholds))
 dimnames(MAT) <- list(NULL, as.character(thresholds))

 l <- 0
 for (d in dates) {
# BMA Brier Scores

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

       M <- is.na(f) | Wmiss
     
       MEAN <- apply(rbind(1, f) * fit$biasCoefs[,,d], 2, sum)

       W <- WEIGHTS
       if (any(M)) { 
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

       MAT[i,] <- sapply( thresholds, cdfBMAnormal,
                         WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M]) -
                         (y[i] <= thresholds)

    }

 }

# locations at which forecasts are made (depends on training length and lag)

 bsBMA <- apply(MAT^2, 2, mean, na.rm = TRUE)
 
 safeDiv <- function(x,y) {
              yzero <- !y
              nz <- sum(yzero)
              result <- rep(NA, length(y))
              if (!nz) result <- x/y else result[!yzero] <- x[!yzero]/y[!yzero]
              result
            }  

# data.frame(thresholds = thresholds,
#            ensemble = 1 - safeDiv(bsVoting,bsClimatology), 
#            logistic = 1 - safeDiv(bsLogistic,bsClimatology),  
#            bma = 1 - safeDiv(bsBMA,bsClimatology))

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting,
            bma = bsBMA)
}

