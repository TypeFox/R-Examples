# Copyright (C) 2008 Jean-Pierre Gattuso and Heloise Lavigne
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"pHconv" <- function (flag=1, pH=8.100, S=35, T=25, P=0, ks="d")

{

n <- max(length(pH), length(flag), length(S), length(T), length(S), length(P), length(ks))

if(length(pH)!=n){pH <- rep(pH[1],n)}
if(length(flag)!=n){flag <- rep(flag[1],n)}
if(length(S)!=n){S <- rep(S[1],n)}
if(length(T)!=n){T <- rep(T[1],n)}
if(length(P)!=n){P <- rep(P[1],n)}
if(length(ks)!=n){ks <- rep(ks[1],n)}

pHconv <- rep(NA, n)
conv <- rep(NA, n)

for(i in (1:n)){  

	Ks = Ks(S=S[i], T=T[i], P=P[i], ks=ks[i])[1]                 # on free pH scale
        Cl = S[i] / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
        ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
	Kf = Kf(S=S[i], T=T[i], P=P[i], pHscale="F")[1]  # on the free pH scale
        FT = 6.7e-5 * Cl/18.9984    # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)

# flag 1 : seawater to total
if (flag[i]==1) {pHconv[i] <- pH[i] + log10(1 + ST/Ks + FT/Kf) - log10(1 + ST/Ks); conv[i] <- "seawater to total"} 
# flag 2 : free to total 
if (flag[i]==2) {pHconv[i] <- pH[i] - log10(1 + ST/Ks); conv[i] <- "free to total"} 
# flag 3 :  total to seawater
if (flag[i]==3) {pHconv[i] <- pH[i] - log10(1 + ST/Ks + FT/Kf) + log10(1 + ST/Ks); conv[i] <- "total to seawater"}
# flag 4 : total to free 
if (flag[i]==4) {pHconv[i] <- pH[i] + log10(1 + ST/Ks); conv[i] <- "total to free"} 
# flag 5 : seawater to free
if (flag[i]==5) {pHconv[i] <- pH[i] + log10(1 + ST/Ks + FT/Kf); conv[i] <- "seawater to free"}
# flag 6 : free to seawater
if (flag[i]==6) {pHconv[i] <- pH[i] - log10(1 + ST/Ks + FT/Kf); conv[i] <- "free to seater"}
}

attr(pHconv, "conversion") <- conv
return(pHconv)

}

