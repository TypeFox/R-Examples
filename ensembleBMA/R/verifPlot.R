verifPlot <-
function( fit, ensembleData, dates=NULL) 
{

 forc <- quantileForecast( fit, ensembleData, dates = dates,
                              quantiles = c( .1, .5, .9))

 if (is.null(dates)) dates <- names(fit$nIter)
 
 use <- as.logical(match(ensembleValidDates(ensembleData), dates, nomatch = 0))

 lon <- ensembleData$longitude[use]; lat <- ensembleData$latitude[use]
 verif <- dataVerifObs(ensembleData)[use]

 nObs <- length(verif)
 if (nObs == 0) stop("no observations")

 ord <- order(forc[,"0.9"])
 ylim <- c(0, max(forc,verif))
 plot(1:nObs, forc[,"0.9"][ord], ylim = ylim, type = "l", col = "black",
      xlab = "Observations in order of increasing 90th percentile forecast",
      ylab = "Precipitation (hundreths of an inch)", xaxt = "n")
 lines(1:nObs, forc[,"0.1"][ord], type = "l", col = "black")
 lines(1:nObs, forc[,"0.5"][ord], type = "l", col = "red")
 points(1:nObs, verif, pch = 16, col = "black", cex = 0.5)

 invisible(forc)
}

