# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


"Kn" <- 
function (S=35, T=25, P=0, pHscale="T")

#--------------------------------------------------------------
# Dissociation constant of ammonium 
#--------------------------------------------------------------

{

nK <- max(length(S), length(T), length(P), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1],nK)}

#-------Constantes----------------

tk = 273.15           # [K] (for conversion [deg C] <-> [K])
TC = T + tk           # TC [C]; T[K]

#--------------------------------------------------------------
# Dissociation constant of ammonium on seawater scale - Millero 1995
#--------------------------------------------------------------
lnK <- -6285.33/TC+0.0001635*TC-0.25444+(0.46532-123.7184/TC)* sqrt(S)+
      (-0.01992+3.17556/TC)*S
Kn <- exp(lnK)

# ----------------- Pressure Correction ------------------	
Kn <- Pcorrect(Kvalue=Kn, Ktype="Kn", T=T, S=S, P=P, pHscale="SWS")


###----------------pH scale corrections
factor <- rep(NA,nK)
pHsc <- rep(NA,nK)
for(i in (1:nK)){   
 if(pHscale[i]=="T"){factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total ; pHsc[i] <- "total scale"}
 if(pHscale[i]=="F"){factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHsc[i] <- "free scale"}
 if(pHscale[i]=="SWS"){factor[i] <- 1 ; pHsc[i] <- "seawater scale"}
Kn[i] <- Kn[i]*factor[i]
}

##------------Warnings

for(i in 1:nK){
if((T[i]>45)|(S[i]>45)|(T[i]<0)){warning("S and/or T is outside the range of validity of the formulation available for Kn in seacarb.")}
}

  attr(Kn,"pH scale") = pHsc
  attr(Kn,"unit")     = "mol/kg-soln"
  return(Kn)

}  ## END Kn
