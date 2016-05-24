matchDates <-
function(fitDates, ensDates, dates=NULL) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 fsub <- rep( TRUE, length(fitDates))
 esub <- rep(TRUE,length(ensDates))

 if (!is.null(dates)) {

## restrict dates to specified dateses
   dates <- sort(unique(as.character(dates)))
   M <- match2side(ensDates, dates)
   if (any(M$yx == 0)) stop("specified dates not present in data")
   if (any(M$xy == 0)) esub <- as.logical(M$xy)

# restrict dates and fitDates to specified dateses

   M <- match2side( dates, fitDates)
   if (any(M$xy == 0)) stop("specified dates not present in fitDates")
   if (any(M$yx == 0)) {
     if (all(M$yx == 0)) stop("fit dates dates not represented in data")
     warning("some fit dates dates not represented in data")
     fsub <- as.logical(M$yx)
   }

 }
 else {

# restrict data and fitDates to specified dates

   M <- match2side( ensDates, fitDates)
   esub <- if (any(M$xy == 0)) as.logical(M$xy) else rep(TRUE,length(M$xy))
   if (any(M$yx == 0)) {
     if (all(M$yx == 0)) stop("fit dates dates not represented in data")
     warning("some fit dates dates not represented in data")
     fsub <- as.logical(M$yx)
   }

 }

 list( fit = fsub, ens = esub)
}

