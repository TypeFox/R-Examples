`matchITandFH.default` <-
function( fit, ensembleData) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 naNULL <- function(x) {
            if (is.na(x)) x <- -2^20
            if (is.null(x)) x <- -2^25
            if (is.infinite(x)) x <- -2^30
            x
           }

 fitFH <- naNULL(attr(fit, "forecastHour"))
 fitIT <- naNULL(attr(fit, "initializationTime"))
 datFH <- naNULL(attr(fit, "forecastHour"))
 datIT <- naNULL(attr(fit, "initializationTime"))

 out <- TRUE

 if (fitFH != datFH & fitIT != datIT) {
   warning("forecast hour and initialization time inconsistent in data and fit\n")
   out <- FALSE
 }
 else if (fitFH != datFH) {
   warning("forecast hour inconsistent in data and fit\n")
   out <- FALSE
 }
 else if (fitIT != datIT) {
   warning("initialization time inconsistent in data and fit\n")
   out <- FALSE
 }

 out
}

