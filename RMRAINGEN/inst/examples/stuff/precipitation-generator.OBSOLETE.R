# file PRECIPITATION_19582010-generator.R
#
# This file contains a script example with a PRECIPITATION_19582010 stochastic generation 
#
#
# author: Emanuele Cordano on 17-12-2011

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
library(RMRAINGEN)
library(mvtnorm)

#load(trentino_19582010)
load('/Users/ecor/Dropbox/iasma/RMAWGENdev/RMRAINGEN/inst/doc/examples/data/trentino_19582010.rda')

year_min <- 1971
origin <- "1-1-1971"

year_max <- 2000
months <- c("Jul")
station <- c("T0001","T0083","T0090","T0099","T0010","T0064")

vstation <- c("B2440","B6130","B8570","B9100","LAVIO","POLSA","SMICH","T0001","T0010","T0014","T0018","T0032","T0064","T0083","T0090","T0092","T0094","T0099","T0102","T0110","T0129","T0139","T0147","T0149","T0152","T0157","T0168","T0179","T0189","T0193","T0204","T0210","T0211","T0327","T0367","T0373")
station <- vstation #[-c(20,33)] #vstation[c(1:19,21:32,34:36)] #[c(1:19,21:25)] stazione 20 con problemi 

NREALIZATION <- 10000
prec_mes <- as.data.frame(extractyears(PRECIPITATION_19582010,year_min=year_min,year_max=year_max,station=station))
prec_mes2 <- extractmonths(prec_mes,when=months,origin=origin)

continuity_ratio <- continuity_ratio(prec_mes2,valmin=0.5)

R <- omega_inv(continuity_ratio$nooccurence,tolerance=0.005)

x <-rmvnorm(NREALIZATION,mean=array(0,ncol(R)),sigma=R,method="eigen")

O <- omega(R,p0_v1=diag(continuity_ratio$nooccurence))

CCGamma <- CCGamma(data=prec_mes,lag=0:2,tolerance=0.001,only.matrix=FALSE)






# FARE OMEGA DI X DANDO ANCHE COME INPUT LA DIAGONLE DI continuity_ratio$nooccurence
#occurence <- prec_mes>0.5
#occurence2 <- diff(occurence)
#pwd <- length(occurence[i,!is.na(occurence) & !is.na(occurences) & occurence!=occurence2])/length(occurence[i,])
# CoSTRUIRE MARKOV CHAIN RAIN O NORAIN
