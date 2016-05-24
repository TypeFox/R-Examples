###############################################################################
# R (http://r-project.org/) Instrument Class Model
#
# Copyright (c) 2009-2012
# Peter Carl, Dirk Eddelbuettel, Jeffrey Ryan, 
# Joshua Ulrich, Brian G. Peterson, and Garrett See
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: MonthCodes.R 899 2012-01-01 19:00:09Z gsee $
#
###############################################################################

#' Month-to-Code and Code-to-Month
#'
#' Convert month code (used for futures contracts) 
#' to abbreviated month name, or convert abbreviated month name to month code
#' @aliases C2M M2C
#' @param code Month code: F, G, H, J, K, M, N , Q, U, V, X, or Z
#' @param month Abbreviated month: jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, or dec
#' @return corresponding code or month.
#' @author Garrett See
#' @seealso \code{\link{MC2N}}
#' @examples
#' C2M()
#' C2M("M")
#' C2M()[6]
#' M2C()
#' M2C("Sep")
#' M2C()[9]
#' @export
C2M <- function(code) {
    if (missing(code)) c(F='Jan',G='Feb',H='Mar',
                    J='Apr',K='May',M='Jun',
                    N='Jul',Q='Aug',U='Sep',
                    V='Oct',X='Nov',Z='Dec')    
    else switch(toupper(code), F='Jan', G='Feb',H='Mar',
                    J='Apr',K='May',M='Jun',
                    N='Jul',Q='Aug',U='Sep',
                    V='Oct',X='Nov',Z='Dec')
}

#' @export
#' @rdname C2M
M2C <- function(month) {
    if (missing(month)) c(jan='F',feb='G',mar='H',
                    apr='J',may='K',jun='M',
                    jul='N',aug='Q',sep='U',
                    oct='V',nov='X',dec='Z')
    else switch(toupper(month), JAN=, JANUARY='F',
                FEB=, FEBRUARY='G',MAR=, MARCH='H', 
                APR=, APRIL='J', MAY='K', JUN=, JUNE='M',
                JUL=, JULY='N', AUG=, AUGUST='Q',
                SEP=, SEPTEMBER='U', OCT=, OCTOBER='V',
                NOV=, NOVEMBER='X', DEC=, DECEMBER='Z')
}


