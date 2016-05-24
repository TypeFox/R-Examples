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
"pTA" <-
function(flag, sys=0, var1, var2, pCO2a, co3, hco3, S=35, T=20, P=0, Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74"){
  # if the concentrations of total silicate and total phosphate are NA
  # they are set at 0
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  
  if (sys==0) {
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
		alkf <- ci$ALK + 2*co3 + hco3 # final alkalinity
		dicf <- ci$DIC + co3 + hco3	# final dic
		cf <- carb(flag=15, var1=alkf, var2=dicf, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
		co <- as.data.frame(c("pTA-closed-initial", rep("pTA-closed-final", nrow(cf))))
	}
	if (sys==1) {
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks,  pHscale=pHscale, b=b)
		alkc <- ci$ALK + 2*co3 + hco3 # total alkalinity before requilibration
		dicc <- ci$DIC + co3 + hco3	# dic before requilibration
		cc <- carb(flag=15, var1=alkc, var2=dicc, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
		alkf <- cc$ALK	# final total alkalinity
		cf <- carb(flag=24, var1=pCO2a, var2=alkf, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
		co <- as.data.frame(c("pTA-open-initial", rep("pTA-open-final", nrow(cf))))
	}
out <- rbind(ci, cf)
out <- cbind(co, out)
names(out)[1] <- "comment"
return(out)
}
