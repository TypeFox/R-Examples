# Copyright (C) 2008 Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# sys is 0 if the system is closed and 1 if it is open
# co3 and hco3 are the   amount added in mol/kg


# One adds a certain volume (vol) of HCl (vol  is then negative) or NaOH (vol is then positive) to 1 kg of seawater. 
# vol is in liter
# N is the normality of the HCl or NaOH
# HCl or NaOH are, respectively, a strong acid and a strong base which are therefore fully dissociated. The concentration of H+ or OH- is N mole/l
# HCl 1N: 0.1 ml therefore adds 0.1 e-3 mol of H+, decreasing TA
# NaOH 1N: 0.1 ml therefore adds 0.1 e-3 mol of OH-, increasing TA
# DIC is constant in a closed system
"ppH" <-
function(flag, sys, var1, var2, pCO2a, vol, N, S=35, T=20, P=0, Pt=0, Sit=0, pHscale="T", k1k2='x', kf='x', ks="d"){
  # if the concentrations of total silicate and total phosphate are NA
  # they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  if (sys==0) {
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit,  pHscale=pHscale, k1k2=k1k2, kf=kf, ks=ks)
		alkf <- (ci$ALK + vol*N)/(1+abs(vol)) #final alk - dilution is taken into account
		dicf <- (ci$DIC)/(1+abs(vol))  #final dic - dilution is taken into account
		cf <- carb(flag=15, var1=alkf, var2=dicf, S=S, T=T, P=P,  Pt=Pt, Sit=Sit, pHscale=pHscale, k1k2=k1k2, kf=kf, ks=ks)
		co <- as.data.frame(c("ppH-closed-initial", rep("ppH-closed-final", nrow(cf)))) 
	}
	if (sys==1) {
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit,  pHscale=pHscale, k1k2=k1k2, kf=kf, ks=ks)
		alkf <- (ci$ALK + vol*N)/(1+abs(vol)) # final total alkalinity  - dilution is taken into account
		dicf <- (ci$DIC)/(1+abs(vol))	# final dic  before requilibration - dilution is taken into account
		#pHf <- ci$pH + ci$PhiH * (-vol) *N	# final pH using a buffer factor (see Frankignoulle, 1994)
		cc <- carb(flag=15, var1=alkf, var2=dicf, S=S, T=T, P=P, Pt=Pt, Sit=Sit,  pHscale=pHscale, k1k2=k1k2, kf=kf, ks=ks) #  before requilibration	
		cf <- carb(flag=24,var1=pCO2a, var2=alkf, S=S, T=T, P=P, Pt=Pt, Sit=Sit,  pHscale=pHscale, k1k2=k1k2, kf=kf, ks=ks)		
		co <- as.data.frame(c("ppH-open-initial", rep("ppH-open-final", nrow(cf))))
	}
	out <- rbind(ci, cf)
	out <- cbind(co, out)
	names(out)[1] <- "comment"
	return(out)
}
