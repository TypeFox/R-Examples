##
## Copyright (c) 2010, Brandon Whitcher
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: $
##

#' Convert DICOM Time/Date Entry
#' 
#' The DICOM time entry (TM) is converted into two alternative formats: a text
#' version of the original format and a number in seconds.  The DICOM date
#' entry (DA) is converted into a simple alternative format.
#' 
#' DICOM \dQuote{TM} format consists of a string of characters of the format
#' \code{hhmmss.frac}; where \code{hh} contains hours (range \dQuote{00} -
#' \dQuote{23}), \code{mm} contains minutes (range \dQuote{00} - \dQuote{59}),
#' \code{ss} contains seconds (range \dQuote{00} - \dQuote{59}), and
#' \code{frac} contains a fractional part of a second as small as 1 millionth
#' of a second (range 000000 - 999999).  A 24 hour clock is assumed.  Midnight
#' can be represented by only 0000 since 2400 would violate the hour range.
#' The string may be padded with trailing spaces.  Leading and embedded spaces
#' are not allowed.  One or more of the components \code{mm}, \code{ss}, or
#' \code{frac} may be unspecified as long as every component to the right of an
#' unspecified component is also unspecified.  If \code{frac} is unspecified
#' the preceding \dQuote{.} may not be included.  \code{Frac} shall be held to
#' six decimal places or less to ensure its format conforms to the ANSI HISPP
#' MSDS Time common data type.  Examples: \describe{ \item{1.}{070907.0705
#' represents a time of 7 hours, 9 minutes and 7.0705 seconds.} \item{2.}{1010
#' represents a time of 10 hours, and 10 minutes.} \item{3.}{021 is an invalid
#' value.} } Notes: For reasons of backward compatibility with versions of this
#' standard prior to V3.0, it is recommended that implementations also support
#' a string of characters of the format \code{hh:mm:ss.frac} for this VR.
#' 
#' DICOM \dQuote{DA} format A string of characters of the format yyyymmdd;
#' where yyyy shall contain year, mm shall contain the month, and dd shall
#' contain the day. This conforms to the ANSI HISPP MSDS Date common data type.
#' Example: \describe{ \item{1.}{19930822 would represent August 22, 1993.} }
#' Notes: For reasons of backward compatibility with versions of this standard
#' prior to V3.0, it is recommended that implementations also support a string
#' of characters of the format yyyy.mm.dd for this VR.
#' 
#' @aliases str2time str2date
#' @param tt TM field from a DICOM header.
#' @param dd DA field from a DICOM header.
#' @param format.in,format.out Appropriate formatting of input or output.
#' @return For \dQuote{TM}, a list structure containing two fields \item{txt}{A
#' text version of the time where colons have been inserted for readability.}
#' \item{time}{Time in seconds from midnight.} for \dQuote{DA}, a simple
#' character string.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{readDICOM}}
#' @references Digital Imaging and Communications in Medicine (DICOM)\cr
#' \url{http://medical.nema.org}\cr
#' \url{http://en.wikipedia.org/wiki/Digital_Imaging_and_Communications_in_Medicine}
#' 
#' @examples 
#' 
#' str2date("19930822")
#' str2time("112308")
#' 
#' @keywords misc
#' @export str2time
str2time <- function(tt, format.out="%02i:%02i:%08.5f") {
  tt <- as.numeric(tt)
  hh <- as.integer(trunc(tt / 10000))
  tt <- tt %% 10000
  mm <- as.integer(trunc(tt / 100))
  ss <- tt %% 100
  list(txt = sprintf(format.out, hh, mm, ss), time = 3600*hh + 60*mm + ss)
}
#' @rdname str2time
#' @export str2date
str2date <- function(dd, format.in="%Y%m%d", format.out="%d %b %Y") {
  format(as.Date(dd, format.in), format.out)
}
