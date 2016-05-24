`getExchangeable` <-
function (argument, attribute, nForecasts) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (is.null(argument)) return(attribute)
 if (length(argument) != nForecasts) {
   stop("exchangeable specification not consistent with forecasts")
 }
 namAtr <- names(attribute)
 if (is.null(namAtr))
   stop("exchangeable data attribute has no member names")
 namArg <- names(argument)
 if (is.null(namArg)) {
   warning("no member names associated with exchangeable specification")
   names(argument) <- names(attribute)
   return(argument)
 }
 m <- match( namArg, namAtr, nomatch = 0)
 if (any(!m) || (length(unique(m)) != length(m))) {
   stop("exchangeable argument names inconsistent with data")
 }
 argument[order(m)]
}

