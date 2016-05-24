MAE.default <-
function(fit, ensembleData, dates=NULL, ...) 
{
 weps <- 1.e-4

 matchITandFH(fit,ensembleData)

 ensembleData <- ensembleData[,matchEnsembleMembers(fit,ensembleData)]
 
## remove instances missing data

 ignoreDates  <- substring(class(fit)[1], 1, 6) == "fitBMA"
#ignoreDates  <- inherits( fit, "fitBMA")

 if (ignoreDates) {

   M <- !dataNA(ensembleData, dates = FALSE)
   if (!all(M)) ensembleData <- ensembleData[M,]

 }
 else {
  
   M <- !dataNA(ensembleData)
   if (!all(M)) ensembleData <- ensembleData[M,]

   fitDates <- modelDates(fit)

   M <- matchDates( fitDates, ensembleValidDates(ensembleData), dates)

   if (!all(M$ens)) ensembleData <- ensembleData[M$ens,]
   if (!all(M$fit)) fit <- fit[fitDates[M$fit]]
}

 obs <- dataVerifObs(ensembleData)
 nObs <- length(obs)

 nForecasts <- ensembleSize(ensembleData) 

# NA removal and date matching has already been taken care of
 Q <- as.vector(quantileForecast( fit, ensembleData))
 if (any(is.na(Q))) stop("NAs in forecast")
 
 ensembleData <- ensembleForecasts(ensembleData)

## maeCli <- mean(abs(obs - median(obs)))

## maeEns <- mean(abs(obs - apply(ensembleData, 1, median)))

 maeCli <- mean(abs(obs - mean(obs)))

 maeEnsMean <- mean(abs(obs - apply(ensembleData, 1, mean, na.rm = TRUE)))
 maeEnsMedian <- mean(abs(obs - apply(ensembleData, 1, median, na.rm = TRUE)))
 mseEnsMean <- sum((obs - apply(ensembleData, 1, mean, na.rm = TRUE))^2)/nObs
 mseEnsMedian <- sum((obs - apply(ensembleData, 1, median, na.rm = TRUE))^2)/nObs

 maeBMAmedian <- mean(abs(obs - Q), na.rm = TRUE) 
 mseBMAmedian <- sum((obs - Q)^2)/length(obs)

##c(climatology = maeCli, ensemble = maeEns, BMA = maeBMA)

 c(ensemble = maeEnsMedian, BMA = maeBMAmedian)
}

