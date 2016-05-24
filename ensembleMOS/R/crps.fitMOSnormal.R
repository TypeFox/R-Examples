crps.fitMOSnormal <-
function(fit, ensembleData, dates=NULL, nSamples = NULL, seed=NULL,  ...)
{

 if(!is.null(dates)) warning("dates ignored")

 #erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

 crpsFunc <- function(mu, sig, y)
  {
   z <- (y - mu)/sig
   crps <- sig * (z*(2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi))
   crps
  }

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts or obs

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]

 if (is.null(obs <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

#nObs <- length(obs)
 nObs <- ensembleNobs(ensembleData)

 if (!is.null(seed)) set.seed(seed)

 nForecasts <- ensembleSize(ensembleData)

 CRPS <- crpsSim <- sampleMedian <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 members <- ensembleMemberLabels(ensembleData)

 obs <- ensembleVerifObs(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)

 B <- fit$B

 if (!all(Bmiss <- is.na(B))) {

    A <- fit$a
    C <- fit$c
    D <- fit$d

    for (i in 1:nObs) {

       f <- ensembleData[i,]
       S.sq <- var(f)
       f <- c(1,f)
       Mu <- c(A,B)%*%f
       Sig <- sqrt(C + D*S.sq)

       CRPS[i] <- crpsFunc(Mu, Sig, obs[i])

   }

}
CRPS
}

