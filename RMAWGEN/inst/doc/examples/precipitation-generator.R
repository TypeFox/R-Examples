# file precipitation-generator-enviro.R
#
# This file contains a script example with a precipitation stochastic generation 
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
data(trentino)

year_max <- 1990
year_min <- 1961
origin <- "1961-1-1"

n_GPCA_iter <- 10
nscenario=20
station <- c("T0090","T0083") #,"T0099","T0001") 

generation00 <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=3,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=0,sample="monthly",nscenario=nscenario,no_spline=FALSE)#VARselect(generation00$var@GPCA_data$final_results)

		

DJF <- extractmonths(data=1:nrow(generation00$prec_gen),when=c("Dec","Jan","Feb"),origin="1961-1-1")
MAM <- extractmonths(data=1:nrow(generation00$prec_gen),when=c("Mar","Apr","May"),origin="1961-1-1")
JJA <- extractmonths(data=1:nrow(generation00$prec_gen),when=c("Jun","Jul","Aug"),origin="1961-1-1")
SON <- extractmonths(data=1:nrow(generation00$prec_gen),when=c("Sep","Oct","Nov"),origin="1961-1-1")

istation <- "T0090"
iseason <- DJF
title <- "DJF"
cex <- 0.5
pch <- 1
lag <- 2

prec <- data.frame(prec_mes=generation00$prec_mes[,istation],prec_gen=generation00$prec_gen[,istation])

for (i in 1:nscenario) {
	
	str <- sprintf("prec_gen%05d",i)
	prec[,str] <- generation00[[str]][,istation]
	print(str)
	
}


# check observed vs generated with Wilcoxon test
wilcox.test(generation00$prec_mes[iseason,istation],generation00$prec_gen[iseason,istation])

# check observed vs generated with Kolgomorov-Smirnov
ks.test(generation00$prec_mes[iseason,istation],generation00$prec_gen[iseason,istation],exact=FALSE)

# CHECK lagged-averaged 

lag <- 1 
titlelag <- paste(title,lag,"day lag",sep=" ")


qqplot.lagged(x=prec,when=iseason,lag=lag,main=titlelag,xlab="observed",ylab="generated",cex=cex,pch=pch)
abline(0,1)






#prec_gen <- extractmonths(data=generation00$prec_gen,when=c("Jun","Jul","Aug"),origin="1961-1-1")[,i]
#prec_mes <- extractmonths(data=generation00$prec_mes,when=c("Jun","Jul","Aug"),origin="1961-1-1")[,i]


# GAUSSIANIZED 

qqplot(generation00$data_prec[iseason,istation],generation00$data_prec_gen[iseason,istation],main=title,xlab="observed",ylab="generated",cex=cex,pch=pch)
abline(0,1)
# FARE ISTOGRAMMI!!!
#c_mes <- continuity_ratio(data_gen,valmin=1.0)
#c_gen <- continuity_ratio(data_mes,valmin=1.0)

#print(generation00$var,rmax=2)

#prec_gen <- generation00$prec_gen[date_gen,2] 
#prec_mes <- generation00$prec_mes[date_mes,2]
#qqplot(prec_mes,prec_gen)
#abline(0,1)



