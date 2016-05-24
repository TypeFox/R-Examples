brierScore.fitBMAnormal <-
function(fit, ensembleData, thresholds, dates=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,dates=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=NULL)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 y <- dataVerifObs(ensembleData)
 nObs <- length(y)

 nForecasts <- ensembleSize(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)
 x <- apply( ensembleData, 1, mean, na.rm = TRUE)
 
 MAT <-  outer(y, thresholds, "<=")

 bsClimatology <- apply(sweep(MAT, MARGIN = 2, FUN = "-", 
            STATS = apply(MAT,2,mean))^2, 2, mean)

 bsVoting <- apply((t(apply(ensembleData, 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean, na.rm = TRUE),
                 thresholds = thresholds)) - MAT)^2, 2, mean, na.rm = TRUE)

# fit doesn't have a training period so logistic fit to all data
 logisticFit <- sapply( thresholds, 
            function(thresh, x, y) 
             glm((y <= thresh) ~ x,family=binomial(logit))$coef,
             x = x, y = y)

 logisticFit[2,][is.na(logisticFit[2,])] <- 0

 MAT <- apply(logisticFit, 2, function(coefs,x) 
                      sapply(coefs[1] + coefs[2]*x, inverseLogit),
                      x = x) - outer(y, thresholds, "<=")

 bsLogistic <- apply(MAT^2, 2, mean)

 MAT <- matrix( NA, nObs, length(thresholds))
 dimnames(MAT) <- list(NULL, as.character(thresholds))

# BMA Brier Scores

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

       MAT[i,] <- sapply( thresholds, cdfBMAnormal,
                         WEIGHTS = W, MEAN = MEAN[!M], SD = SD[!M]) -
                         (y[i] <= thresholds)

    }

 }


# locations at which forecasts are made (depends on training length and lag)

 bsBMA <- apply(MAT^2, 2, mean, na.rm = TRUE)
 
# safeDiv <- function(x,y) {
#              yzero <- !y
#              nz <- sum(yzero)
#              result <- rep(NA, length(y))
#             if (!nz) result <- x/y else result[!yzero] <- x[!yzero]/y[!yzero]
#             result
#           }  

# data.frame(thresholds = thresholds,
#            ensemble = 1 - safeDiv(bsVoting,bsClimatology), 
#            logistic = 1 - safeDiv(bsLogistic,bsClimatology),  
#            bma = 1 - safeDiv(bsBMA,bsClimatology))

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting,
            bma = bsBMA)
}

