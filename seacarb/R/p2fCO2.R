# Copyright (C) 2010  Héloïse Lavigne and Jean-Pierre Gattuso
# with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr>
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr>
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
#
p2fCO2 <- function(T=25, Patm=1, P=0, pCO2){
tk <- 273.15;           # [K] (for conversion [deg C] <-> [K])
TK <- T + tk;           # TK [K]; T[C]
Phydro_atm = P / 1.01325  # convert hydrostatic pressure from bar to atm (1.01325 bar / atm)
Ptot = Patm + Phydro_atm  # total pressure (in atm) = atmospheric pressure + hydrostatic pressure

# Original "seacarb" f2pCO2 calculation:
# B <- (-1636.75+12.0408*TK-0.0327957*(TK*TK)+0.0000316528*(TK*TK*TK))*1e-6
# fCO2 <-  pCO2*(1/exp((1*100000)*(B+2*(57.7-0.118*TK)*1e-6)/(8.314*TK)))^(-1)
# Above calculation:
# - uses incorrect R (wrong units, incompatible with pressure in atm)
# - neglects a term "x2" (see below)
# - assumes pressure is always 1 atm (wrong for subsurface)

# To compute fugcoeff, we need 3 other terms (B, Del, xc2) in addition to 3 others above (TK, Ptot, R)
  B   <- -1636.75 + 12.0408*TK - 0.0327957*TK^2 + 0.0000316528*TK^3
  Del <- 57.7-0.118*TK

# "x2" term often neglected (assumed = 1) in applications of Weiss's (1974) equation 9
# x2 = 1 - x1 = 1 - xCO2 (it is close to 1, but not quite)
# Let's assume that xCO2 = pCO2. Resulting fugcoeff is identical to at least 8th digit after the decimal.
  xCO2approx <- pCO2
  xc2 <- (1 - xCO2approx*1e-6)^2 

  fugcoeff = exp( Ptot*(B + 2*xc2*Del)/(82.057*TK) )
  fCO2 <- pCO2 * fugcoeff


return(fCO2)
}
