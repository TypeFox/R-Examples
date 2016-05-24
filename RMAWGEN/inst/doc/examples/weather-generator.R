# file weather-generator.R
# 
# This file contains a script example with two coupled temperature and precipitation stochastic generations 
#
#
# author: Emanuele Cordano on 12-01-2012
#
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
set.seed(1222)
library(RMAWGEN)
data(trentino)

year_max <- 1990
year_min <- 1961
origin <- "1961-1-1"

n_GPCA_iter <- 10
n_GPCA_iteration_residuals <- 10
n_GPCA_iter_prec <- 20
n_GPCA_iteration_residuals_prec <- 20
station <- c("T0090","T0083","T0099","T0001") 

#vstation <- c("B2440","B6130","B8570","B9100","LAVIO","POLSA","SMICH","T0001",
#		"T0010","T0014","T0018","T0032","T0064","T0083","T0090","T0092","T0094","T0099",
#		"T0102","T0110","T0129","T0139","T0147","T0149","T0152","T0157","T0168","T0179","T0189","T0193","T0204","T0210","T0211","T0327","T0367","T0373")		

# generation of temperature max and min 
generation00_temp <- ComprehensiveTemperatureGenerator(station=station,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,p=1,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,sample="monthly")

# Use of measured and observed temperature as exogenous variables

exogen_sim <- cbind(generation00_temp$output$Tx_gen,generation00_temp$output$Tn_gen)
names(exogen_sim) <- cbind(paste(names(generation00_temp$output$Tx_gen),"_Tx",sep=""),paste(names(generation00_temp$output$Tn_gen),"_Tn",sep=""))
exogen <- cbind(generation00_temp$input$Tx_mes,generation00_temp$input$Tn_mes)
names(exogen) <- cbind(paste(names(generation00_temp$input$Tx_mes),"_Tx",sep=""),paste(names(generation00_temp$input$Tn_mes),"_Tn",sep=""))

# Precipitation Generator (temperture enters as exogenous variable)

valmin <- 1.0
generation00_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,exogen=exogen,exogen_sim=exogen_sim,p=1,n_GPCA_iteration=n_GPCA_iter_prec,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals_prec,sample="monthly",valmin=1.0,extremes=TRUE)

# Post-processing calculation on precipitation 

prec_mes <- generation00_prec$prec_mes
prec_gen <- generation00_prec$prec_gen

vprec_mes <- prec_mes[,1]
vprec_gen <- prec_gen[,1]

qqplot(vprec_mes[!is.na(vprec_mes) & vprec_mes>valmin],vprec_gen[vprec_gen>valmin & !is.na(vprec_gen)],xlab="measured",ylab="generated",main=paste("Q-Qplot precipitation at ",names(vprec_gen),sep=""))
qqplot(vprec_mes,vprec_gen,xlab="measured",ylab="generated",main=paste("Q-Qplot precipitation at ",names(vprec_gen),sep=""))
mes <- length(vprec_mes[!is.na(vprec_mes) & vprec_mes>0])/length(vprec_mes[!is.na(vprec_mes)])
gen <- length(vprec_mes[!is.na(vprec_gen) & vprec_gen>0])/length(vprec_gen[!is.na(vprec_gen)])

data_gen <- extractmonths(data=generation00_prec$prec_gen,when=c("Jun","Jul","Aug"),origin="1961-1-1")
data_mes <- extractmonths(data=generation00_prec$prec_mes,when=c("Jun","Jul","Aug"),origin="1961-1-1")

c_mes <- continuity_ratio(data_gen,valmin=1.0)
c_gen <- continuity_ratio(data_mes,valmin=1.0)

print(generation00_temp$var)
plot_sample(vprec_mes,vprec_gen,sample="monthly",origin=origin,sort=TRUE)
