crps.ensembleBMAgamma <-
function(fit, ensembleData, dates=NULL, nSamples=10000, seed=NULL, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 if (!is.null(seed)) set.seed(seed)

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 if (is.null(nSamples)) nSamples <- 10000

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

 Q <- as.vector(quantileForecast( fit, ensembleData, dates = dates))
 if (any(is.na(Q))) stop("NAs in forecast") # fix like ensembleBMAgamma0

 obs <- dataVerifObs(ensembleData)
 nForecasts <- ensembleSize(ensembleData) 

 crpsSim <- rep(NA, nObs)
 names(crpsSim) <- dataObsLabels(ensembleData)

 members <- ensembleMembers(ensembleData)
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

       RATE <- MEAN/VAR
       SHAPE <- MEAN*RATE

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M] / sum(W[!M])
       }

       if (sum(!M) > 1) {
         SAMPLES <- sample( (1:nForecasts)[!M], size = nSamples,
                             replace = TRUE, prob = W) 
       }
       else {
         SAMPLES <- rep((1:nForecasts)[!M], nSamples)
       }

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

       nz <- SAMPLES != 0
       if (any(nz)) SAMPLES[nz] <- sapply(SAMPLES[nz], powinv, power=fit$power)

# crps2 approximates a term that is quadratic in the number of members
       crps1  <- mean(abs(SAMPLES - obs[i])) 
       crps2 <-  mean(abs(diff(sample(SAMPLES))))
       crpsSim[i]  <- crps1 - crps2/2
    }
 }

##crpsSim <- mean(crpsSim, na.rm = TRUE)

  crpsCli <- sapply(obs, function(x,Y) mean(abs(Y-x)), Y = obs)
##crpsCli <- mean(crpsCli - mean(crpsCli)/2)
  
  crpsCli <- crpsCli - mean(crpsCli)/2

 crpsEns1 <- apply(abs(sweep(ensembleData,MARGIN=1,FUN ="-",STATS=obs))
                   ,1,mean,na.rm=TRUE)

 if (nrow(ensembleData) > 1) {
   crpsEns2 <- apply(apply(ensembleData, 2, function(z,Z) 
     apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum,na.rm=TRUE),
                  Z = ensembleData),1,sum, na.rm = TRUE)
 }
 else {
   crpsEns2 <- sum(sapply(as.vector(ensembleData), 
                   function(z,Z) sum( Z-z, na.rm = TRUE),
                   Z = as.vector(ensembleData)), na.rm = TRUE)
 }

##crpsEns <- mean(crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts)))
  crpsEns <- crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts))

#cbind(climatology = crpsCli, ensemble = crpsEns, BMA = crpsSim)
 cbind(ensemble = crpsEns, BMA = crpsSim)
}

