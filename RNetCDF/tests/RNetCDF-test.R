#===============================================================================#
#                                                                               #
#  Name:       RNetCDF-test.R                                                   #
#                                                                               #
#  Version:    1.8-2                                                          #
#                                                                               #
#  Purpose:    Test functions to the NetCDF interface for R.                    #
#                                                                               #
#  Author:     Pavel Michna (michna@giub.unibe.ch)                              #
#              Milton Woods (m.woods@bom.gov.au)                                #
#                                                                               #
#  Copyright:  (C) 2010-2014 Pavel Michna                                       #
#                                                                               #
#===============================================================================#
#                                                                               #
#  This program is free software; you can redistribute it and/or modify         #
#  it under the terms of the GNU General Public License as published by         #
#  the Free Software Foundation; either version 2 of the License, or            #
#  (at your option) any later version.                                          #
#                                                                               #
#  This program is distributed in the hope that it will be useful,              #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of               #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
#  GNU General Public License for more details.                                 #
#                                                                               #
#  You should have received a copy of the GNU General Public License            #
#  along with this program; if not, write to the Free Software                  #
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    #
#                                                                               #
#===============================================================================#
#  Implementation and Revisions                                                 #
#-------------------------------------------------------------------------------#
#  Author   Date       Description                                              #
#  ------   ----       -----------                                              #
#  pm       29/12/10   First implementation                                     #
#  mw       18/07/12   Test packed variables                                    #
#  mw       02/09/14   Test 1D character arrays and character scalars           #
#  mw       05/09/14   Test reading/writing NC_CHAR as raw bytes                #
#  mw       26/01/16   Test utcal.nc and utinvcal.nc with POSIXct type          #
#                                                                               #
#===============================================================================#


#===============================================================================#
#  Load library                                                                 #
#===============================================================================#

library(RNetCDF)


#===============================================================================#
#  Run tests                                                                    #
#===============================================================================#

#-------------------------------------------------------------------------------#
#  NetCDF library functions                                                     #
#-------------------------------------------------------------------------------#

#--Initialize ------------------------------------------------------------------#
cat("Starting NetCDF tests...\n")

testfun <- function(x,y,tally=NULL) {
  if (is.null(tally)) {
    tally <- c(pass=0,fail=0)
  }
  # Compare numeric values with single precision tolerance:
  if (isTRUE(all.equal(x,y,tolerance=2^(-23)))) {
    cat("OK\n")
    return(tally+c(1,0))
  } else {
    cat("failed\n")
    return(tally+c(0,1))
  }
}

##  Create a new NetCDF dataset and define dimensions
nc <- create.nc("foo.nc")

nstation <- 5
ntime <- 2
nstring <- 32
nempty <- 0

dim.def.nc(nc, "station", nstation)
dim.def.nc(nc, "time", ntime)
dim.def.nc(nc, "max_string_length", nstring)
dim.def.nc(nc, "empty", unlim=TRUE)

##  Define variables
var.def.nc(nc, "time", "NC_INT", "time")
var.def.nc(nc, "temperature", "NC_DOUBLE", c(0,1))
var.def.nc(nc, "packvar", "NC_BYTE", c("station"))
var.def.nc(nc, "name", "NC_CHAR", c("max_string_length", "station"))
var.def.nc(nc, "qcflag", "NC_CHAR", c("station"))
var.def.nc(nc, "int0", "NC_INT", NA)
var.def.nc(nc, "char0", "NC_CHAR", NA)
var.def.nc(nc, "numempty", "NC_FLOAT", c("station","empty"))

##  Put some missing_value attribute for temperature
att.put.nc(nc, "temperature", "missing_value", "NC_DOUBLE", -99999.9)

## Define the packing used by packvar
att.put.nc(nc, "packvar", "scale_factor", "NC_DOUBLE", 10)
att.put.nc(nc, "packvar", "add_offset", "NC_DOUBLE", -5)

##  Define variable values
mytime        <- c(1:2)
mytemperature <- matrix(c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, NA, NA, 9.9),ncol=ntime)
mypackvar     <- seq_len(5)*10-5
myname        <- c("alfa", "bravo", "charlie", "delta", "echo")
myqcflag      <- "ABCDE"
myint0        <- 12345
mychar0       <- "?"

##  Put the data
var.put.nc(nc, "time", mytime, 1, length(mytime))
var.put.nc(nc, "temperature", mytemperature, c(1,1), c(nstation,ntime))
var.put.nc(nc, "packvar", mypackvar, pack=TRUE)
var.put.nc(nc, "name", myname, c(1,1), c(nstring,nstation))
var.put.nc(nc, "qcflag", charToRaw(myqcflag))
var.put.nc(nc, "int0", myint0)
var.put.nc(nc, "char0", mychar0)

sync.nc(nc)

## Read tests
cat("Read numeric vector ... ")
x <- mytime
dim(x) <- length(x)
y <- var.get.nc(nc, 0)
tally <- testfun(x,y)

cat("Read numeric matrix ... ")
x <- mytemperature
y <- var.get.nc(nc, "temperature")
tally <- testfun(x,y,tally)

cat("Read numeric matrix slice ... ")
x <- mytemperature[,2]
dim(x) <- length(x)
y <- var.get.nc(nc, "temperature", c(NA,2), c(NA,1))
tally <- testfun(x,y,tally)

cat("Read numeric matrix empty slice ... ")
x <- numeric(0)
dim(x) <- c(0,1)
y <- var.get.nc(nc, "temperature", c(NA,2), c(0,1),collapse=FALSE)
tally <- testfun(x,y,tally)

cat("Read numeric scalar ... ")
x <- myint0
dim(x) <- 1
y <- var.get.nc(nc, "int0")
tally <- testfun(x,y,tally)

cat("Read numeric empty array ... ")
x <- numeric(0)
dim(x) <- c(nstation,nempty)
y <- var.get.nc(nc, "numempty")
tally <- testfun(x,y,tally)

cat("Read 2D char array ... ")
x <- myname
dim(x) <- length(x)
y <- var.get.nc(nc, "name")
tally <- testfun(x,y,tally)

cat("Read 2D char slice ... ")
x <- substring(myname[2:3],1,4)
dim(x) <- length(x)
y <- var.get.nc(nc, "name", c(1,2), c(4,2))
tally <- testfun(x,y,tally)

cat("Read 2D char slice as raw bytes ... ")
x <- substring(myname[2:3],1,4)
dim(x) <- length(x)
x <- apply(x,MARGIN=1,FUN=charToRaw)
y <- var.get.nc(nc, "name", c(1,2), c(4,2), rawchar=TRUE)
tally <- testfun(x,y,tally)

cat("Read 2D char slice as characters ... ")
x <- myname[2:3]
dim(x) <- length(x)
y <- var.get.nc(nc, "name", c(1,2), c(NA,2))
tally <- testfun(x,y,tally)

cat("Read empty 2D char array ... ")
x <- character(0)
dim(x) <- 0
y <- var.get.nc(nc, "name", NA, c(0,0),collapse=FALSE)
tally <- testfun(x,y,tally)

cat("Read 1D char slice ... ")
x <- substring(myqcflag,2,3)
dim(x) <- 1
y <- var.get.nc(nc, "qcflag", c(2), c(2))
tally <- testfun(x,y,tally)

cat("Read scalar char ... ")
x <- mychar0
dim(x) <- 1
y <- var.get.nc(nc, "char0")
tally <- testfun(x,y,tally)

cat("Read and unpack numeric array ... ")
x <- mypackvar
dim(x) <- length(x)
y <- var.get.nc(nc, "packvar", unpack=TRUE)
tally <- testfun(x,y,tally)

#-- Close file -----------------------------------------------------------------#
close.nc(nc)


#-------------------------------------------------------------------------------#
#  UDUNITS calendar functions                                                   #
#-------------------------------------------------------------------------------#

cat("utcal.nc - numeric values ...")
x <- matrix(data=c(1899, 1900, 1900, 1900, 1900, 1900,
                     12,    1,    1,    1,    1,    1,
		     31,    1,    1,    1,    1,    1,
		     23,    0,    1,    2,    3,    4,
		      0,    0,    0,    0,    0,    0,
		      0,    0,    0,    0,    0,    0),
	    ncol=6)
colnames(x) <- c("year","month","day","hour","minute","second")
y <- utcal.nc("hours since 1900-01-01 00:00:00 +01:00", c(0:5))
tally <- testfun(x,y,tally)

cat("utcal.nc - string values ...")
x <- c("1899-12-31 23:00:00", "1900-01-01 00:00:00", "1900-01-01 01:00:00",
       "1900-01-01 02:00:00", "1900-01-01 03:00:00", "1900-01-01 04:00:00")
y <- utcal.nc("hours since 1900-01-01 00:00:00 +01:00", c(0:5), type="s")
tally <- testfun(x,y,tally)

cat("utcal.nc - POSIXct values ...")
x <- ISOdatetime(c(1899,1900,1900,1900,1900,1900),
                 c(  12,   1,   1,   1,   1,   1),
                 c(  31,   1,   1,   1,   1,   1),
                 c(  23,   0,   1,   2,   3,   4),
                 c(   0,   0,   0,   0,   0,   0),
                 c(   0,   0,   0,   0,   0,   0), tz="UTC")
y <- utcal.nc("hours since 1900-01-01 00:00:00 +01:00", c(0:5), type="c")
tally <- testfun(x,y,tally)

cat("utinvcal.nc - numeric values ...")
x <- 6.416667
y <- utinvcal.nc("hours since 1900-01-01 00:00:00 +01:00", c(1900,1,1,5,25,0))
tally <- testfun(x,y,tally)

cat("utinvcal.nc - string values ...")
x <- 6.416667
y <- utinvcal.nc("hours since 1900-01-01 00:00:00 +01:00", "1900-01-01 05:25:00")
tally <- testfun(x,y,tally)

cat("utinvcal.nc - POSIXct values ...")
x <- 6.416667
y <- utinvcal.nc("hours since 1900-01-01 00:00:00 +01:00",
         ISOdatetime(1900,1,1,5,25,0,tz="UTC"))
tally <- testfun(x,y,tally)

#-------------------------------------------------------------------------------#
#  Overall summary                                                              #
#-------------------------------------------------------------------------------#
cat("Summary:", tally["pass"], "pass /", tally["fail"], "fail. ")

if (tally["fail"]==0) {
  cat("Package seems to work properly.\n")
} else {
  cat("Some problems were detected.\n")
}

#===============================================================================#

#===============================================================================#
#  SCRATCH									#
#===============================================================================#

