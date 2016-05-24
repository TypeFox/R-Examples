# file marginal_gaussianization1.R
#
# This file contains a trivial example of One-dimensional margial Gaussianization
#
#
# author: Emanuele Cordano on 05-06-2012

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





# This script plots the measurement sites on a GoogleMap support 

rm(list=ls())
library(RMAWGEN)
z <- rexp(10000,rate=0.5) 
u <- z
x <- normalizeGaussian(x=z,data=z)
z1 <- normalizeGaussian(x=x,data=z,inverse=TRUE) 
#-- Plot probability histogram of z, x and z1 
zhist <- hist(z,probability=TRUE,breaks=50)
xhist <- hist(x,probability=TRUE,breaks=50)
z1hist <- hist(z1,probability=TRUE,breaks=50)
zhist$counts <- zhist$density
xhist$counts <- xhist$density
z1hist$counts <- z1hist$density
def.par <- par(no.readonly = TRUE) # save default, for resetting...
par(mfrow=c(1,3),oma=c(15.0,0.0,15.0,0.0)) 
plot(zhist)
plot(xhist)
plot(z1hist)
par(def.par)
# end of the script

library("RMAWGEN") 
data("trentino") 
col <- rainbow(n=12,start=0.1,end=0.9) 
col[6:1] <- col[1:6]
col[7:12] <- col[12:7] 

pdf <- "/Users/ecor/R-packages/RMAWGENCodeCorner/script_for_plotting/plot/temperature_gaussianization.pdf"

pdf(pdf)
plot_sample(x=TEMPERATURE_MIN$T0090,sample="monthly",
origin="1958-1-1",axes=FALSE,
xlab="Tn [degC]",ylab="x",abline=NULL,col=col)
###plot_sample(x=TEMPERATURE_MIN\$T0090,sample="monthly",\\ %origin="1958-1-1",axes=FALSE,xlab="Tn [degC]",ylab="x")\\
dev_off()

pdf(pdf)

pdf <- "/Users/ecor/R-packages/RMAWGENCodeCorner/script_for_plotting/plot/precipitation_gaussianization.pdf"

plot_sample(x=PRECIPITATION$T0090,sample="monthly",
origin="1958-1-1",axes=FALSE,xlab="Prec [mm]",
ylab="x",abline=NULL,valmin=0.5,step=0,col=col)

dev.off()
 
