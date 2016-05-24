# file temperature-generator.R
#
# This file contains a script example with a precipitation stochastic generation 
#
#
# author: Emanuele Cordano on 02-03-2012

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
data(trentino)
set.seed(1233)
vstation <- c("B2440","B6130","B8570","B9100","LAVIO","POLSA","SMICH","T0001",
		"T0010","T0014","T0018","T0032","T0064","T0083","T0090","T0092","T0094","T0099",
		"T0102","T0110","T0129","T0139","T0147","T0149","T0152","T0157","T0168","T0179","T0189","T0193","T0204","T0210","T0211","T0327","T0367","T0373")
station <- c("T0090","T0083")
generation0 <- ComprehensiveTemperatureGenerator(station=station,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=1961,year_max=1990)
VAR_model <- generation0$var
CX <- generation0$input$monthly_mean_Tx[,station]
CN <- generation0$input$monthly_mean_Tn[,station]
generation1 <- ComprehensiveTemperatureGenerator(station=station,varmodel=generation0$var,
		onlygeneration=TRUE,year_min=1961,year_max=1990,mean_climate_Tn_sim=CN,mean_climate_Tx_sim=CX,year_min_sim=1961,year_max_sim=1990,
		Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,nscenario=3,)

str(generation0)
str(generation1)
