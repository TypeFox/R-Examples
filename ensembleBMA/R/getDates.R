getDates <-
function( DATES, julianDATES, dates, trainingDays, lag, incr) {
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

  if (length(trainingDays) > 1 || trainingDays <= 0 
       || (trainingDays - trunc(trainingDays)) != 0) 
   stop("trainingDays improperly specified")

 if (trainingDays > length(julianDATES)) 
   stop("insufficient training data")

 Jdates <- seq(from = julianDATES[trainingDays]+lag*incr,
               to = max(julianDATES)+lag*incr, by = incr)

 origin <- attr( julianDATES, "origin")
 lD <- unique(sapply(DATES,nchar))

 modelDates <- julTOymdh(Jdates, origin = origin, dropHour = (lD == 8))

 if (nullDates <- is.null(dates)) {

   dates <- modelDates

 }
 else {

   dates <- sort(unique(as.character(dates)))

   M <- match( dates, modelDates, nomatch = 0)
   if (any(!M)) stop("specified dates not present in data")

   if (any(dates < julTOymdh(min(Jdates),origin=origin,dropHour=(lD == 8)))) {
     stop("dates precede the first training period")
   }

   if (any(dates > julTOymdh(max(Jdates),origin=origin,dropHour=(lD == 8)))) {
     warning("there are dates beyond the last training period")
   }

 }

 dates
}

