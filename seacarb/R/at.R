# Copyright (C) 2011 Jean-Pierre Gattuso and Heloise Lavigne and Steeve Comeau
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#

"at" <-
function(S=35, T=25, C=0.1, d=1, pHTris=NULL, ETris=NULL, weight, E, volume){

# S salinity    (constant)
# T temperature (°C)  (vector or constant)
# C molarity of acid   (constant)
# d density of acid    (constant)
# pHTris, pH used for the calibration of the electrode  (constant)
# ETris, voltage used for the calibration of the electrode (constant)
# E  voltage recorded during the titration in milliVolt(vector)
# volume : volume of acid (vector)
# weight : mass of the sample  (constant)

R<-8.31447215        #constant
F<-96485.339924

# check that T is given as a vector
test <- length(T) != length(E)
if(test) { T <- rep(T[1], length(E)) }

Tk <- T + 273.15

# creation of a table p
p <- data.frame(E=E, volume=volume, Tk=Tk)
z <- p

# transform mV in pH (total scale)
if(!is.null(pHTris)&!is.null(ETris)){
pH <- pHTris + (ETris/1000-E/1000)/(R*Tk*log(10)/F)     

# creation of a table p
p <- data.frame(p, pH=pH)
#to use only the value between pH 3.5 and 3 (according to Dickson, 2007):
iii <- which((3<= p$pH)&(p$pH<=3.5))
z<- p[iii,]
}

options(digits=9)

m <- z$volume*d	# Mass of acid
m0 <- weight		# Mass of the sample

#linear estimation of the total alkalinity (gran function):
#####################################################
F1<-(m0+m)*exp((z$E/1000)/(R*(z$Tk)/F))
f<-lm(m~F1)
TA<-f$coefficients[1]*C/m0[1]


#non linear estimation:
E0 <- z$E/1000-(R*z$Tk/F)*log((-m0*TA+m*C)/(m0+m))
Hprime <- exp((z$E/1000-E0)/(R*z$Tk/F))

Cl <- S / 1.80655;             # Cl = chlorinity; S = salinity (per mille)
St <- 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
Ks <- Ks(S,T[1], 0)
Z <- 1+ St/Ks
Ft <- 6.7e-5 * Cl/18.9984      # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
Kf <- exp(874/z$Tk-9.68+0.111*S**(1/2))
y <- (m/m0)
regr <- nls(y ~ ((At + (St/(1 + Ks*Z/(f*Hprime)))+(Ft/ (1+Kf/(f*Hprime)))+(f*Hprime/Z))/(C-f*Hprime/Z)),start=list(At=TA, f=1))

ALK <- summary(regr)$parameters[1]
attr(ALK,"unit") <- "mol/kg-soln"
attr(ALK,"name") <- "Total Alkalinity"

## At= Total Alkalinity (mol/kg)

return(ALK)
}
