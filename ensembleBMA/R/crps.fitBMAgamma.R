crps.fitBMAgamma <-
function(fit, ensembleData, dates=NULL, nSamples=NULL, seed=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 if (!is.null(dates)) warning("dates ignored")

 if (is.null(nSamples)) nSamples <- 10000

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]

 M <- !dataNA(ensembleData,dates=FALSE)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(fit)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates=NULL)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) fit <- fit[fitDates[M$fit]]

 dates <- modelDates(fit)

 Dates <- ensembleValidDates(ensembleData)

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 if (!is.null(seed)) set.seed(seed)

 nForecasts <- ensembleSize(ensembleData) 

 crpsSim <- rep(NA, nObs)
 names(crpsSim) <- dataObsLabels(ensembleData)

 members <- ensembleMembers(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 WEIGHTS <- fit$weights

 if (!all(Wmiss <- is.na(WEIGHTS))) {

    for (i in 1:nObs) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss

       VAR <- (fit$varCoefs[1] + fit$varCoefs[2]*f)^2

       fTrans <- sapply(f, powfun, power = fit$power)

       MEAN <- apply(rbind(1, fTrans) * fit$biasCoefs, 2, sum)

       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M] / sum(W[!M])
       }

       SAMPLES <- sample( (1:nForecasts)[!M], size = nSamples,
                           replace = TRUE, prob = W) 

       tab <- rep(0, nForecasts)
       names(tab) <- members
       for (j in seq(along = tab)) tab[j] <- sum(SAMPLES == j)
       
       SAMPLES[] <- NA

       jj <- 0
       for (j in seq(along = tab)) {
          nsamp <- tab[j]
          if (nsamp == 0) next
          SAMPLES[jj + 1:nsamp] <- rgamma(nsamp,shape=SHAPE[j],rate=RATE[j])
          jj <- jj + nsamp
        }

# crps2 approximates a term that is quadratic in the number of members
       nz <- SAMPLES != 0
       if (any(nz)) SAMPLES[nz] <- sapply(SAMPLES[nz], powinv, power=fit$power)
  
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
    }
 }

 crpsCli <- sapply(obs, function(x,Y) mean(abs(Y-x)), Y = obs)
 crpsCli <- crpsCli - mean(crpsCli)/2

 crpsEns1 <- apply(abs(sweep(ensembleData,MARGIN=1,FUN ="-",STATS=obs)),
                   1, mean, na.rm = TRUE)
 crpsEns2 <- apply(apply(ensembleData, 2, function(z,Z) 
       apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum,na.rm=TRUE),
                  Z = ensembleData), 1, sum, na.rm = TRUE)
 crpsEns <- crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts))

#cbind(climatology = crpsCli, ensemble = crpsEns, BMA = crpsSim)
 cbind(ensemble = crpsEns, BMA = crpsSim)
}

