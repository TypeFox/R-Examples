crps.ensembleMOSnormal <-
function(fit, ensembleData, dates=NULL, nSamples=10000, seed=NULL, ...)
{
 crpsFunc <- function(mu, sig, y)
  {
   z <- (y - mu)/sig
   crps <- sig * (z*(2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi))
   crps
  }
 matchITandFH(fit,ensembleData)

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]

## match specified dates with dateTable in fit

 dateTable <- dimnames(fit$B)[[2]]

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable))
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K))
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

  }

 ensDates <- ensembleValidDates(ensembleData)

## match dates in data with dateTable
 if (is.null(ensDates) || all(is.na(ensDates))) {
   if (length(dates) > 1) stop("date ambiguity")
   nObs <- nrow(ensembleData)
   Dates <- rep( dates, nObs)
 }
 else {
## remove instances missing dates
   if (any(M <- is.na(ensDates))) {
     ensembleData <- ensembleData[!M,]
     ensDates <- ensembleValidDates(ensembleData)
   }
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L))
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }

 Q <- as.vector(quantileForecast(fit, ensembleData, dates = dates))
 if (any(is.na(Q))) stop("NAs in forecast") # fix like ensembleBMAgamma0

 obs <- ensembleVerifObs(ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 members <- ensembleMemberLabels(ensembleData)

 CRPS <- crpsSim <- sampleMedian <- rep(NA, nObs)
 names(crpsSim) <- names(sampleMedian) <- ensembleObsLabels(ensembleData)

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (d in dates) {

    l <- l + 1
    k <- K[l]

    B <- fit$B[,k]
    if (all(Bmiss <- is.na(B))) next

    A <- fit$a[,k]
    C <- fit$c[,k]
    D <- fit$d[,k]

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {

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

