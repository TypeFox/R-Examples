##
## bdt.R: Some accessor functions for Boost Date_Time functionality
##
## Copyright (C) 2010 - 2012  Dirk Eddelbuettel and Romain Francois
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

getEndOfBizWeek <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getEndOfBizWeek()
}

getEndOfMonth <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getEndOfMonth()
}

getYear <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getYear()
}

getMonth <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getMonth()
}

getDay <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getDay()
}

getDayOfWeek <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getDayOfWeek()
}

getDayOfYear <- function(date = Sys.Date()) {
    stopifnot(inherits(date, "Date"))
    bdt$setDate(date)
    bdt$getDayOfYear()
}

getIMMDate <- function(mon, year) {     # defined as third Wednesday
    bdt$getIMMDate(mon, year)
}

getNthDayOfWeek <- function(nthday, dow, mon, year) {
    bdt$getNthDayOfWeek(nthday, dow, mon, year)
}

getLastDayOfWeekInMonth <- function(nthday, mon, year) {
    bdt$getLastDayOfWeekInMonth(nthday, mon, year)
}

getFirstDayOfWeekInMonth <- function(nthday, mon, year) {
    bdt$getFirstDayOfWeekInMonth(nthday, mon, year)
}

getFirstDayOfWeekAfter <- function(dow, date) {
    bdt$getFirstDayOfWeekAfter(dow, date)
}

getLastDayOfWeekBefore <- function(dow, date) {
    bdt$getLastDayOfWeekBefore(dow, date)
}
