# Copyright (C) 2015 James C. Orr
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"theta" <- 
function (S=35, T=25, P=0, Pref=0)

#--------------------------------------------------------------
# Potential temperature
#--------------------------------------------------------------

{
    nK <- max(length(S), length(T), length(P), length(Pref))

    ##-------- If any input argument is a vector, make vectors out of others having only a single entry

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(P)!=nK){Pref <- rep(Pref[1], nK)}

    #The 'oce' package's 'swTheta' routine is used to compute potential temperature

    #Convert pressure in bar (seacarb default) to decibar (oce swTheta default)
    pressure = 10*P
    referencePressure = 10*Pref

#   Compute potential temperature on ITS 90 scale (with routine that requires IPTS 68 scale)
#   - IPTS 68 = International Practical Temperature Scale of 1968
#   - ITS 90  = International Temperature Scale of 1990
#   (for T scale conversion: see Dickson et al., 2007, Chap. 5, p. 7, including footnote)
#   -------------------------------------------------------------------------------------
#   a) Compute "in-situ Temperature" (ITS 90) from "in situ Temperature" (IPTS 68)
    T68 <- (T - 0.0002) / 0.99975

#   b) Compute "potential temperature" (IPTS68) from "in situ temperature" (IPTS 68)
    theta68 <- swTheta(salinity=S, temperature=T68, pressure=pressure, referencePressure=referencePressure, "unesco")

#   c) Compute potential temp (ITS 90) from potential temp (IPTS 68)
    theta90 <- 0.99975*theta68 + 0.0002

#   Note: parts (a) and (c) above are tiny corrections (often neglected)
#         part  (b) is a big correction for deep waters (but zero at surface)

    attr(theta90,"unit")     = "degrees C (ITS 90)"
    return(theta90)
}  ## END theta
