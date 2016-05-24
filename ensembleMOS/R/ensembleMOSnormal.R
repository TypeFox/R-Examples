ensembleMOSnormal <-
function(ensembleData, trainingDays, consecutive = FALSE, dates = NULL,
         control = controlMOSnormal(), warmStart = FALSE,
         exchangeable = NULL)
{
  if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")

 call <- match.call()

 if(!is.logical(warmStart)) stop("warmStart improperly specified")

 if (!is.logical(consecutive)) stop("consecutive improperly specified")

 if (is.list(trainingDays)) trainingsDays <- trainingsDays[[1]]

 if (length(trainingDays) > 1 || trainingDays <= 0
       || (trainingDays - trunc(trainingDays)) != 0)
   stop("trainingDays improperly specified")

 forecastHour <- ensembleFhour(ensembleData)
 lag <- ceiling( forecastHour / 24 )

 ensMemNames <- ensembleMemberLabels(ensembleData)
 nForecasts <- length(ensMemNames)

 exchangeable <- getExchangeable( exchangeable,
                                 ensembleGroups(ensembleData),nForecasts)

# remove instances missing all forecasts, obs or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 M <- M | is.na(ensembleValidDates(ensembleData))
 ensembleData <- ensembleData[!M,]

 nObs <- ensembleNobs(ensembleData)
 if (!nObs) stop("no observations")

 ensDates <- ensembleValidDates(ensembleData)
 if (is.null(ensDates)) stop("dates unavailable")

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
   warning("valid dates do not have a unique forecast hour")

 lD <- nchar(DATES[1])

 if (!(lD <- unique(sapply(DATES,nchar))))
    stop("all dates in data should have same character length")

 if (nullDates <- is.null(dates)) {

   dates <- julTOymdh(Jdates, origin = origin, dropHour = (lD == 8))

 }
 else {

   dates <- sort(unique(as.character(dates)))

   if (!all(dateCheck(dates)))
     stop("improperly specified date(s) in dates argument")

   datesHH <- getHH(dates)

   if (length(datesHH) != 1)
     warning("dates do not have a unique forecast hour")

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

 }

 juliandates <- ymdhTOjul( dates, origin = origin)

 nDates <- length(dates)

 a <- array(NA, c(1,nDates))
 dimnames(a) <- list("a", dates)

 B <- array( NA, c(nForecasts, nDates))
 dimnames(B) <- list(ensMemNames, dates)

 c <- array( NA, c(1, nDates))
 dimnames(c) <- list(c("c"), dates)

 d <- array( NA, c(1, nDates))
 dimnames(d) <- list(c("d"), dates)

 trainTable <- rep(0, nDates)
 names(trainTable) <- dates

 nIter <- rep(0, nDates)
 names(nIter) <- dates

 L <- length(juliandates)
 twin <- 1:trainingDays

## temp <- data.frame(julian = julianDATES,date = DATES)
## print(temp)

 cat("\n")

 l <- 0
 for(i in seq(along = juliandates)) {
    cat("modeling for date", dates[i], "...")

    if(!consecutive){
      I <- (juliandates[i]-lag*incr) >= julianDATES
      if (!any(I)) stop("insufficient training data")

      j <- which(I)[sum(I)]

      if (j != l) {

        twin <- (j+1) - (1:trainingDays)

        D <- as.logical(match(Dates, DATES[twin], nomatch=0))
             if (!any(D)) stop("this should not happen")

        fit <- fitMOSnormal(ensembleData[D,], control = control,
                          exchangeable = exchangeable)

      }
   }else{

      dMax <- DATES[which(juliandates[i]  == julianDATES)- lag]
      dMin <- DATES[which(juliandates[i]  == julianDATES) - lag -
                          trainingDays + 1]
      tDATES <- DATES[as.logical((DATES <= dMax) * (DATES >= dMin))]
      D <- as.logical(match(Dates, tDATES, nomatch=0))
              if (!any(D)) stop("this should not happen")

      fit <- fitMOSnormal(ensembleData[D,], control = control,
                          exchangeable = exchangeable)
   }
    trainTable[i] = sum(D)
    c[,i] <- fit$c
    d[,i] <- fit$d
    a[i] <- fit$a
    B[,i] <- fit$B
    if (warmStart) {
      control$start <- list(a = fit$a, B = fit$B, c = fit$c, d = fit$d)
    }
    cat("\n")
    print(round(c(fit$a, fit$B),2))
    print(round(c(fit$c, fit$d),2))
    cat("\n")
  }


 structure(list(training = c(days = trainingDays, lag = lag,
                table = trainTable),
		a = a, B = B, c = c, d = d,
                exchangeable = attr(ensembleData, "exchangeable")),
                call = match.call(), class = "ensembleMOSnormal")
}

