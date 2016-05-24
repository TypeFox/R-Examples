trainingData <-
function(ensembleData, trainingDays, consecutive = FALSE, date = NULL)
{
 if (!is.logical(consecutive)) stop("consecutive mispecified")

 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 if (is.list(trainingDays)) trainingsDays <- trainingDays[[1]]

 if (length(trainingDays) > 1 || trainingDays <= 0
       || (trainingDays - trunc(trainingDays)) != 0)
   stop("trainingDays improperly specified")

 forecastHour <- ensembleFhour(ensembleData)
 lag <- ceiling(forecastHour / 24)

 ensMemNames <- ensembleMemberLabels(ensembleData)
 nForecasts <- length(ensMemNames)

# remove instances missing all forecasts, obs or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 M <- M | is.na(ensembleValidDates(ensembleData))
 ensembleData <- ensembleData[!M,]

 nObs <- ensembleNobs(ensembleData)
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

   dates <- as.character(date)

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

    if(!consecutive){
      I <- (juliandates[i]-lag*incr) >= julianDATES
      if (!any(I)) stop("insufficient training data")

      j <- which(I)[sum(I)]

      if (j != l) {

        twin <- (j+1) - (1:trainingDays)

        D <- as.logical(match(Dates, DATES[twin], nomatch=0))
        if (!any(D)) stop("this should not happen")

        trainData <- ensembleData[D,]
      }
   }else{

      dMax <- DATES[which(juliandates[i]  == julianDATES)- lag]
      dMin <- DATES[which(juliandates[i]  == julianDATES) - lag -
                          trainingDays + 1]
      tDATES <- DATES[as.logical((DATES <= dMax) * (DATES >= dMin))]
      D <- as.logical(match(Dates, tDATES, nomatch=0))
              if (!any(D)) stop("this should not happen")

        trainData <- ensembleData[D,]
   }


 }

  trainData

}

