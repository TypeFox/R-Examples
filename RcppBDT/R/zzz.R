
## zzz.R: Loading Rcpp and Boost Date_Time glue
##
## Copyright (C) 2010 - 2013  Dirk Eddelbuettel and Romain Francois
##
## This file is part of RcppBDT.
##
## RcppBDT is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppBDT is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppBDT.  If not, see <http://www.gnu.org/licenses/>.


loadModule("bdtDdMod", TRUE)
loadModule("bdtDtMod", TRUE)
loadModule("bdtDuMod", TRUE)
loadModule("bdtPtMod", TRUE)
loadModule("bdtTzMod", TRUE)

## create a variable 'bdt' from the bdtDt Module (formerly: bdtMod)
## this variable is used as a package-global instance in some R access function
delayedAssign("bdt", local({
    x <- new(bdtDt)
    x$setFromUTC()
    x
}))


## define an onLoad expression to set some methods
evalqOnLoad({
    setMethod("show", "Rcpp_bdtDt", function(object) print(object$getDate()))
    setMethod("show", "Rcpp_bdtTz", function(object) print(object$getRegion()))
    setMethod("show", "Rcpp_bdtPt", function(object) print(object$getDatetime()))
    setMethod("show", "Rcpp_bdtDu", function(object) print(as.difftime(object$getTotalSeconds() + object$getFractionalSeconds() / 1.0e9, units="secs")))
    setMethod("show", "Rcpp_bdtDd", function(object) print(as.difftime(object$getDays(), units="days")))

    setGeneric("format", function(x,...) standardGeneric("format"))
    setMethod("format", "Rcpp_bdtDt", function(x, ...) format(x$getDate(), ...))
    setMethod("format", "Rcpp_bdtTz", function(x, ...) format(x$getRegion(), ...))
    setMethod("format", "Rcpp_bdtPt", function(x, ...) format(x$getDatetime(), ...))
    setMethod("format", "Rcpp_bdtDu", function(x, ...) format(as.difftime(x$getTotalSeconds() + x$getFractionalSeconds() / 1.0e9, units="secs")))
    setMethod("format", "Rcpp_bdtDd", function(x, ...) format(as.difftime( x$getDays(), units="days")))

    setMethod("Arith",   signature(e1 = "Rcpp_bdtDu", e2 = "Rcpp_bdtDu" ), function(e1, e2) arith_bdtDu_bdtDu(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDu", e2 = "integer" ),    function(e1, e2) arith_bdtDu_int(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDu", e2 = "numeric" ),    function(e1, e2) arith_bdtDu_int(e1, as.integer(e2), .Generic))
    setMethod("Arith",   signature(e1 = "integer",    e2 = "Rcpp_bdtDu"),  function(e1, e2) arith_int_bdtDu(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "numeric",    e2 = "Rcpp_bdtDu"),  function(e1, e2) arith_int_bdtDu(as.integer(e1), e2, .Generic))
    setMethod("Compare", signature(e1 = "Rcpp_bdtDu", e2 = "Rcpp_bdtDu"),  function(e1, e2) compare_bdtDu_bdtDu(e1, e2, .Generic))

    setMethod("Arith",   signature(e1 = "Rcpp_bdtPt", e2 = "Rcpp_bdtDu" ), function(e1, e2) arith_bdtPt_bdtDu(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDu", e2 = "Rcpp_bdtPt" ), function(e1, e2) arith_bdtDu_bdtPt(e1, e2, .Generic))
    setMethod("Compare", signature(e1 = "Rcpp_bdtPt", e2 = "Rcpp_bdtPt"),  function(e1, e2) compare_bdtPt_bdtPt(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtPt", e2 = "numeric" ),    function(e1, e2) arith_bdtPt_double(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "numeric",    e2 = "Rcpp_bdtPt"),  function(e1, e2) arith_double_bdtPt(e1, e2, .Generic))

    setMethod("Compare", signature(e1 = "Rcpp_bdtDt", e2 = "Rcpp_bdtDt"),  function(e1, e2) compare_bdtDt_bdtDt(e1, e2, .Generic))
    setMethod("Compare", signature(e1 = "Rcpp_bdtDt", e2 = "Date"),        function(e1, e2) compare_bdtDt_bdtDt(e1, new(bdtDt, e2), .Generic))
    setMethod("Compare", signature(e1 = "Date",       e2 = "Rcpp_bdtDt"),  function(e1, e2) compare_bdtDt_bdtDt(new(bdtDt, e1), e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDt", e2 = "integer"),     function(e1, e2) arith_bdtDt_int(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDt", e2 = "numeric"),     function(e1, e2) arith_bdtDt_int(e1, as.integer(e2), .Generic))
    setMethod("Arith",   signature(e1 = "integer",    e2 = "Rcpp_bdtDt"),  function(e1, e2) arith_int_bdtDt(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "numeric",    e2 = "Rcpp_bdtDt"),  function(e1, e2) arith_int_bdtDt(as.integer(e1), e2, .Generic))

    setMethod("Arith",   signature(e1 = "Rcpp_bdtDd", e2 = "Rcpp_bdtDd" ), function(e1, e2) arith_bdtDd_bdtDd(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDd", e2 = "integer" ),    function(e1, e2) arith_bdtDd_int(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDd", e2 = "numeric" ),    function(e1, e2) arith_bdtDd_int(e1, as.integer(e2), .Generic))
    setMethod("Arith",   signature(e1 = "integer",    e2 = "Rcpp_bdtDd"),  function(e1, e2) arith_int_bdtDd(e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "numeric",    e2 = "Rcpp_bdtDd"),  function(e1, e2) arith_int_bdtDd(as.integer(e1), e2, .Generic))
    setMethod("Compare", signature(e1 = "Rcpp_bdtDd", e2 = "Rcpp_bdtDd"),  function(e1, e2) compare_bdtDd_bdtDd( e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDd", e2 = "Rcpp_bdtDt" ), function(e1, e2) arith_bdtDd_bdtDt( e1, e2, .Generic))
    setMethod("Arith",   signature(e1 = "Rcpp_bdtDt", e2 = "Rcpp_bdtDd" ), function(e1, e2) arith_bdtDt_bdtDd( e1, e2, .Generic))
  })






