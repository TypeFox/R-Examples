# file generations.R
#
# This file contains a script example with an application of software RMAWGEN on temperature daily data of two thermometric stations of Trentino.
#
#
# author: Emanuele Cordano, Annalisa Di Piazza on 21-02-2013

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
library(RMAWGEN)
library(RClimMAWGEN)
# trentino_1958_2010 contains max and min temperature and precipitation data for some site of Trentino
data(trentino_1958_2010)

# cosidered stations
vstation <- c("T0129","T0139")

generation_p1_calc <- ComprehensiveTemperatureGenerator(station=vstation,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,n_GPCA_iteration=0,year_min=1981,year_max=2010,p=1,nscenario=5,seed=456,yearly=FALSE)

# Check generation output with generation_p1

data(generation_p1)

# write the 'str' 

str(generation_p1_calc$output$Tx_gen00004)
## 
##
str(generation_p1$output$Tx_gen00004)
## 
str(generation_p1_calc$output$Tx_gen00004==generation_p1$output$Tx_gen00004)

min((generation_p1_calc$output$Tx_gen00004==generation_p1$output$Tx_gen00004))
## It shoul be returned 1 if the two date frames are equal!!