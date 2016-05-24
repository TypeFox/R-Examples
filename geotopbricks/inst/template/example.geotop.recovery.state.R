# TODO: Add comment
# 
# Author: ecor
###############################################################################

# example.geotop.recovery.state.R
#
# This file contains a script which plots the measurement sites on a GoogleMap support
#
#
# author: Emanuele Cordano on 15-05-2013

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################


rm(list=ls())
library(geotopbricks)
wpath <- system.file('template/friuli',package="geotopbricks")


rec <- get.geotop.inpts.keyword.value("SubfolderRecoveryFiles",wpath=wpath,add_wpath=TRUE)
soilfile <- paste(get.geotop.inpts.keyword.value("SoilParFile",wpath=wpath,add_wpath=TRUE),"0001.txt",sep="") # path to soil type file 
soildata <- read.table(soilfile,sep=",",header=TRUE)
nsoillayers <- nrow(soildata)
recstate <- get.geotop.recovery.state(recFolder=rec,nsoillayers=nsoillayers)
newRecFolder <- "./newRecFolder_justatest"

 # Not Run !!! Please uncomment the following lines to run them. 
 if (!file.exists(newRecFolder)) dir.create(newRecFolder)
 set.geotop.recovery.state(rec=recstate,newRecFolder=newRecFolder)
 
 

 # Not Run !!! Please uncomment the following lines to run them. 
newtextfile <- paste(newRecFolder,"newtextfile.txt",sep="/")
write.vectorized.geotop.recovery(rec=recstate,file=newtextfile,overwrite=TRUE)

newrecstate <- read.vectorized.geotop.recovery(file=newtextfile)

 # Check a generic recovery map!!

it <- newrecstate$names[3]
newrecstate[[it]]==recstate[[it]]


