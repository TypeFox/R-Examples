# file example.generate.R 
#
# This file script contains makes a multi-site genaration of daily precipitation in some sites of the Trentino Dataset
# The data and te results contained in this script are for educational use only and may not be realistic.
#
#
# author: Emanuele Cordano on 06-01-2015
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

library(RGENERATEPREC)

data(trentino)

year_min <- 1961
year_max <- 1990

origin <- paste(year_min,1,1,sep="-")
end <- paste(year_max,12,31,sep="-")

period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
period_temp <- TEMPERATURE_MAX$year>=year_min & TEMPERATURE_MAX$year<=year_max

prec_mes <- PRECIPITATION[period,]
Tx_mes <- TEMPERATURE_MAX[period_temp,]
Tn_mes <- TEMPERATURE_MIN[period_temp,]
accepted <- array(TRUE,length(names(prec_mes)))
names(accepted) <- names(prec_mes)
for (it in names(prec_mes)) {
	acc <- TRUE
	acc <- (length(which(!is.na(Tx_mes[,it])))==length(Tx_mes[,it]))
	acc <- (length(which(!is.na(Tn_mes[,it])))==length(Tn_mes[,it])) & acc
	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it])) & acc
	
}

valmin <- 1.0
prec_mes <- prec_mes[,accepted]



Tx_mes <- Tx_mes[,accepted]
Tn_mes <- Tn_mes[,accepted]
prec_occurence_mes <- prec_mes>=valmin

station <- names(prec_mes)[!(names(prec_mes) %in% c("day","month","year"))]
it <- station[2]
vect <- Tx_mes[,it]-Tn_mes[,it]
months <- factor(prec_mes$month)

#
### Not Run!!!
###  Please uncomment the following lines to run them


model <-
		PrecipitationOccurenceModel(x=prec_mes[,it],exogen=vect,
				monthly.factor=months,valmin=valmin)
#
obs <- prec_mes[,it]>=valmin
#
gen <- generate(model,exogen=vect,monthly.factor=months,n=length(months))


### MultiSite Generation


station <- station[1:2]
exogen <- Tx_mes[,station]-Tn_mes[,station]

months <- factor(prec_mes$month)

#
### Not Run!!!
###  Please uncomment the following lines to run them

model_multisite <-
		PrecipitationOccurenceMultiSiteModel(x=prec_mes[,station],
				exogen=exogen,origin=origin,multisite_type="wilks")

#
## LOGIT-type Model
model_multisite_logit <-
		PrecipitationOccurenceMultiSiteModel(x=prec_mes,exogen=exogen,
				origin=origin,multisite_type="logit",station=station)
#
#
obs_multisite <- prec_mes[,station]>=valmin
#
gen_multisite <- generate(model_multisite,exogen=exogen,origin=origin,end=end)
#
gen_multisite_logit <- generate(model_multisite_logit,exogen=exogen,origin=origin,end=end)
