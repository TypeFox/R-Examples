# Copyright (C) 2008 Karline Soetaert (K.Soetaert@nioo.knaw.nl) and Heloise Lavigne
# Revised by James Orr, 2012-01-17
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


#--------------------------------------------------------------
# Conversion factors for converting dissociation constants 
# from total pH scale to free pH scale (ktotal2free) 
# from total pH scale to seawater pH scale (ktotal2SWS)
# and from free pH scale to seawater scale (kfree2SWS)
# kfree = ktotal * ktotal2free
# kSWS  = ktotal * ktotal2SWS
# kSWS  = kfree  * kfree2SWS
# kfree = kSWS * kSWS2free
# ktotal = kfree * kfree2total
# ktotal = kSWS * kSWS2total
#--------------------------------------------------------------


"kconv" <- function (S=35,T=25,P=0,kf='x',Ks=0,Kff=0)

{
    #nK <- max(length(S), length(T), length(P), length(kf), length(Ks), length(Kff))
    #if(length(S)!=nK){S <- rep(S[1], nK)}
    #if(length(T)!=nK){T <- rep(T[1], nK)}
    #if(length(P)!=nK){P <- rep(P[1], nK)}
    #if(length(kf)!=nK){kf <- rep(kf[1], nK)}
    #if(length(Ks)!=nK){Ks <- rep(Ks[1], nK)}
    #if(length(Kff)!=nK){Kff <- rep(Kff[1], nK)}
    
    if ( missing(kf) && exists("kfg",envir=parent.frame()) )
       {
	 kf <- get("kfg", envir = parent.frame())
       }

	#--------------------------------------------------------------
    # CONVERT equilibrium constants to free scale:
	#--------------------------------------------------------------
	#------------------ Ks ----------------------------------------
	#       Dickson and Goyet (1994), Chapter 5, p.13
	#       (required for total2free)
	#       Equilibrium constant for HSO4- = H+ + SO4--
	#
	#       K_S  = [H+]free [SO4--] / [HSO4-]
	#       pH-scale: free scale !!!
	#
	#       the term log(1-0.001005*S) converts from
	#       mol/kg H2O to mol/kg soln
	
    # if Ks not given
    if (missing(Ks))
	{
	    Ks = Ks(S=S, T=T, P=P)                 # on free pH scale
	}
        Cl = S / 1.80655              # Cl = chlorinity; S = salinity (per mille)
        ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
	total2free  = 1/(1+ST/Ks)      # Kfree = Ktotal*total2free
	total2free <- as.numeric(total2free)	

	#---------------------------------------------------------------------
	# --------------------- Kf  ------------------------------------------
	#  Kf = [H+][F-]/[HF]
	#
	#   (Dickson and Riley, 1979 in Dickson and Goyet,
	#   1994, Chapter 5, p. 14)
	#   pH-scale: 'free'

    # if Kf on free pH scale is not given as an argument 
	if (missing(Kff))
	{
	    Kff = Kf(S=S, T=T, P=P, kf=kf, pHscale="F")
	}


	#------- sws2free -----------------------------------------------
	#
	#       convert from pH_sws ('seawater scale`) to pH ('free`):
	#       pH_sws = pH_free - log(1+S_T/K_S(S,T)+F_T/K_F(S,T))
        Cl = S / 1.80655            # Cl = chlorinity; S = salinity (per mille)
        FT = 6.7e-5 * Cl/18.9984    # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
	free2SWS  = 1+ST/Ks+FT/Kff         # Kfree = Ksws*sws2free
	free2SWS <- as.numeric(free2SWS)
	total2SWS = total2free * free2SWS # KSWS = Ktotal*total2SWS
	total2SWS <- as.numeric(total2SWS)

    return (list(ktotal2SWS=total2SWS, ktotal2free=total2free,kfree2SWS=free2SWS, kfree2total=(1/total2free), kSWS2total=(1/total2SWS), kSWS2free=(1/free2SWS)))
}

