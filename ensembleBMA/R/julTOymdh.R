`julTOymdh` <-
function (julianDates, origin = NULL, dropHour = NULL) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
# require("chron")
 
 if (!is.null(orig <- attr( julianDates, "origin"))) {
   if (!is.null(origin)) {
     bad <- origin["month"] != orig["month"]
     bad <- origin["day"] != orig["day"] || bad
     bad <- origin["year"] != orig["year"] || bad
     if (bad) stop("origin is not uniquely specified")
   }
   origin <- orig
 }
 else if (is.null(origin)) stop("origin is not specified")


 eps <- abs(round(julianDates) - julianDates)

 if (any(eps > 0)) {
   julianDates <- round(julianDates*24)/24
   hour <- round(24*as.vector(julianDates - floor(julianDates)))
 }
 else hour <- 0

 x <- month.day.year( as.vector(floor(julianDates)), origin. = origin)

 if (is.null(dropHour)) {
   l <- attr(julianDates, "nchar")
   dropHour <- is.null(l) || l == 8
 }
 
 if (any(hour != 0) || !dropHour) {
   x <- lapply(c(x[c("year", "month", "day")], list(hour = hour)), 
                 as.character)
 }
 else {
   x <- lapply(x[c("year", "month", "day")], as.character)
 }

 pad0mdh <- function(x) {
 pad0 <- function(x) if (nchar(x) == 2) x else paste("0", x, sep = "")
         as.vector(sapply(x, pad0))
         }

 x[-1] <- lapply(x[-1], pad0mdh)
 as.vector(apply(data.frame(x), 1, paste, collapse = ""))
}

