# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
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
"pHinsi" <- function(pH=8.2, ALK=2.4e-3, Tinsi=20, Tlab=25, S=35, Pt = 0, Sit = 0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74") {
  # if the concentrations of total silicate and total phosphate are NA
  # they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  
  # according to Hunter (1998)
	#si TA=0 alors calculer TA: TA=660+47.6 * S
	dat1 <- carb(flag = 8, var1 = pH, var2 = ALK, S = S, T = Tlab, P = 0, Pt = Pt, Sit = Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
	DIC <- dat1$DIC
	dat2 <- carb(flag = 15, var1 = ALK, var2 = DIC, S = S, T = Tinsi, P = 0, Pt = Pt, Sit = Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
	ph_insi <- dat2$pH
	#utiliser DIC et TA pour calculer pH in situ (flag=15)
	#attr(bor,"unit") <- "mol/kg"
	return(ph_insi)
}
