brierScore.ensembleBMAgamma <-
function(fit, ensembleData, thresholds, dates = NULL, ...) 
{
#
# copyright 2006-present University of Washington. All rights reserved.
# for terms of use, see the LICENSE file.
#
 powfun <- function(x, power) x^power

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

 x <- sapply(apply( ensembleData, 1, mean, na.rm = TRUE), powfun,
             power = fit$power)

 MAT <-  t(outer(y, thresholds, "<="))

 bsClimatology <- apply(sweep(MAT, MARGIN = 1, FUN = "-", 
                        STATS = apply(MAT,1,mean))^2, 1, mean)
 
 bsVotingEns <- apply(ensembleData, 1, function(z, thresholds) 
                 apply(outer(z, thresholds, "<="), 2, mean, na.rm = TRUE),
                 thresholds = thresholds)

 bsVoting <- apply((bsVotingEns - MAT)^2, 1, mean, na.rm = TRUE)

# avoid training data and apply to all data

logisticFit <- sapply( thresholds, 
            function(thresh, x, y) 
             glm((y <= thresh) ~ x,family=binomial(logit))$coef,
             x = x, y = y)

 logisticFit[2,][is.na(logisticFit[2,])] <- 0

 MAT <- apply(logisticFit, 2, function(coefs,x) 
                      sapply(coefs[1] + coefs[2]*x, inverseLogit),
                      x = x) - outer(y, thresholds, "<=")

 bsLogistic <- apply(MAT^2, 2, mean)

 MAT <- matrix( NA, nrow = nObs, ncol = length(thresholds))
 dimnames(MAT) <- list(NULL, as.character(thresholds))

 l <- 0

 for (d in dates) {
# BMA Brier Scores

    l <- l + 1

    WEIGHTS <- fit$weights[,d]
    if (any(Wmiss <- is.na(WEIGHTS)))  next

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

       MAT[i,] <- sapply( sapply( thresholds, powfun, power = fit$power), 
                          cdfBMAgamma, 
          WEIGHTS=W, MEAN=MEAN[!M], VAR=VAR[!M]) - (y[i] <= thresholds)

    }

 }

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
#

 data.frame(thresholds = thresholds,
            climatology = bsClimatology, 
            ensemble = bsVoting, 
            logistic = bsLogistic,
            bma = bsBMA)
}

