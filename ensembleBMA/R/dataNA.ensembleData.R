dataNA.ensembleData <-
function(ensembleData, forecasts = TRUE, observations = TRUE, dates = TRUE)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

# remove instances missing all forecasts, obs or dates

 M <- rep( TRUE, nrow(ensembleData))
 if (forecasts) M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 if (observations & !is.null(dataVerifObs(ensembleData))) {
   nObs <- dataNobs(ensembleData)
   if (nObs > 0) M <- M | is.na(dataVerifObs(ensembleData))
 }
 if (dates & !is.null(ensembleData$dates)) M <- M | is.na(ensembleValidDates(ensembleData))

 M
}

