trainingData <-
function(ensembleData, trainingDays, date = NULL)
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 if (is.list(trainingDays)) trainingsDays <- trainingDays[[1]]

 if (length(trainingDays) > 1 || trainingDays <= 0 
       || (trainingDays - trunc(trainingDays)) != 0) 
   stop("trainingDays improperly specified")
 
 lag <- ceiling( ensembleFhour(ensembleData) / 24 )

 ensMemNames <- ensembleMembers(ensembleData)
 nForecasts <- length(ensMemNames)

 nObs <- dataNobs(ensembleData)
 if (!nObs) stop("no observations")

# remove instances missing all forecasts, obs or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(dataVerifObs(ensembleData))
 M <- M | is.na(ensembleValidDates(ensembleData))
 ensembleData <- ensembleData[!M,]
 
 nObs <- dataNobs(ensembleData)
 if (!nObs) stop("no observations")

 ensDates <- ensembleValidDates(ensembleData)

 Dates <- as.character(ensDates)
 DATES <- sort(unique(Dates))

 if (trainingDays > length(DATES)) 
   stop("insufficient training data")

 julianDATES <- ymdhTOjul(DATES)
 origin <- attr( julianDATES, "origin")
 incr <- min(1,min(diff(julianDATES))) ## incr may be fractional for hours

## dates that can be modeled by the training data (ignoring gaps)

 Jdates <- seq(from = julianDATES[trainingDays]+lag*incr,
               to = max(julianDATES)+lag*incr, by = incr)

## determine the modeling dates

 DATEShh <- getHH(DATES)

 if (length(DATEShh) != 1) 
    warning("valid dates do not have a unique forcast hour")

 if (!(lD <- unique(sapply(DATES,nchar)))) 
   stop("all dates in data should have same character length")

   dates <- date

   if (!all(dateCheck(dates))) 
     stop("improperly specified date(s) in dates argument")

    datesHH <- getHH(dates)

    if (length(datesHH) != 1) 
       warning("dates do not have a unique forcast hour")

    if (any(datesHH != DATEShh)) stop("specified dates incompatible with data")

    if (!(ld <- unique(sapply(dates,nchar)))) 
      stop("all specified dates should have same character length")

    if (ld < lD) {
      dates <- sapply( dates, function(s) paste(s, "00", sep =""))
    }
    else if (ld < lD) {
      dates <- sapply( dates, function(s) substring(s, 1, 8))
    }

    if (any(dates < julTOymdh(min(Jdates),origin=origin,dropHour=(lD == 8)))) {
     stop("some dates precede the first training period")
   }

   if (any(dates > julTOymdh(max(Jdates),origin=origin,dropHour=(lD == 8)))) {
     warning("there are dates beyond the last training period")
   }

 juliandates <- ymdhTOjul( dates, origin = origin)

 K <- 1:nForecasts

 L <- length(juliandates)
 twin <- 1:trainingDays

 cat("\n")
 l <- 0
 for(i in seq(along = juliandates)) {

    I <- (juliandates[i]-lag*incr) >= julianDATES
    if (!any(I)) stop("insufficient training data")

    j <- which(I)[sum(I)]

    if (j != l) {

      twin <- (j+1) - (1:trainingDays)

      D <- as.logical(match(Dates, DATES[twin], nomatch=0))
      if (!any(D)) stop("this should not happen")

      trainData <- ensembleData[D,]
   }

 }

  trainData

}

