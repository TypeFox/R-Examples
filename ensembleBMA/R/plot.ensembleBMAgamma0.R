plot.ensembleBMAgamma0 <-
function(x, ensembleData, dates = NULL, ask = TRUE, ...) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 par(ask = ask)

 powfun <- function(x,power) x^power
 powinv <- function(x,power) x^(1/power)

 weps <- 1.e-4

 matchITandFH(x,ensembleData)

 exchangeable <- x$exchangeable

 ensembleData <- ensembleData[,matchEnsembleMembers(x,ensembleData)]

 M <- !dataNA(ensembleData)
 if (!all(M)) ensembleData <- ensembleData[M,]

 fitDates <- modelDates(x)

 M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates)

 if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
 if (!all(M$fit)) x <- x[fitDates[M$fit]]

 dates <- modelDates(x)

 Dates <- ensembleValidDates(ensembleData)

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 if (nObs == 0) obs <- rep(NA,nrow(ensembleData))

 nForecasts <- ensembleSize(ensembleData)
 
 ensembleData <- ensembleForecasts(ensembleData)

 obs <- powfun( obs, power = x$power)
 
 l <- 0
 for (d in dates) {

    l <- l + 1

    WEIGHTS <- x$weights[,d]
     
    if (all(Wmiss <- is.na(WEIGHTS))) next

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
    
       f <- ensembleData[i,]

       M <- is.na(f) | Wmiss
     
       VAR <- x$varCoefs[1,d] + x$varCoefs[2,d]*f
        
       fTrans <- sapply(f, powfun, power = x$power)

       MEAN <- apply(rbind(1, fTrans) * x$biasCoefs[,,d], 2, sum)

       PROB0 <- sapply(apply(rbind( 1, fTrans, f == 0)*x$prob0coefs[,,d],
                              2,sum), inverseLogit)

       W <- WEIGHTS
       if (any(M)) {
         W <- W + weps
         W <- W[!M]/sum(W[!M])
       }

      plotBMAgamma0(WEIGHTS = W, MEAN = MEAN[!M], VAR = VAR[!M], 
                      PROB0 = PROB0[!M], obs = obs[i],
                      exchangeable = exchangeable, power = x$power)

    }

 }

 invisible()
}

