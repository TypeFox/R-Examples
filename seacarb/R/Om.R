# Copyright (C) 2010 Heloise Lavigne, Andreas Anderson and Jean-Pierre Gattuso
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


"Om" <- function(x, flag, var1, var2, k1k2='x', kf='x', ks="d", pHscale="T", b="u74"){

c <- carb(flag=flag, var1=var1, var2=var2, S=35, T=25, P=0, Pt=0, Sit=0, k1k2= k1k2, kf= kf, ks=ks,  pHscale= pHscale, b=b)

S=35
x <- x*100  #conversion to percentage

#Activity coefficients Millero and Pierrot 1998
gCO3 <- 0.037
gCa <- 0.191
gMg <- 0.2029

Mg <- 0.05282*S/35   #Mg2+ concentration for a seawater at S=35 Dickson 2007
Ca <- 0.01028*S/35    #Dickson 2007
CO3 <- c$CO3

aCa <- gCa*Ca
aMg <- gMg*Mg
aCO3 <- gCO3*CO3

## to find the good IAP (equations derived from Bishoff, 1993)
IAPbiogenic <- 10^(-8.45022726020772 + -0.0363899891218405*x^1 + 0.0107168547236804*x^2 + -0.000307011230719308*x^3)
IAPcleaned <- 10^(-8.34594750694325 + -0.0492826022095371*x^1 + 0.0125595815829352*x^2 + -0.0010724925092009*x^3 + 4.33859823603212e-05*x^4 + -6.69727692992637e-07*x^5)
x <- x/100 # reconversion to proportion between 0 and 1
omegaBiogenic <- (aMg^x)*(aCa^(1-x))*(aCO3)/IAPbiogenic
omegaCleaned <- (aMg^x)*(aCa^(1-x))*(aCO3)/IAPcleaned
RES <- list(OmegaMgCa_biogenic=omegaBiogenic, OmegaMgCa_biogenic_cleaned=omegaCleaned)
return(RES)
}
