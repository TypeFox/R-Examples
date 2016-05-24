# Copyright (C) 2008 Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye
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
"Ks" <-
function(S=35,T=25,P=0, ks="d"){

    nK <- max(length(S), length(T), length(P), length(ks))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(ks)!=nK){ks <- rep(ks[1], nK)}

    #-------Constants----------------

    #---- issues de equic----
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK = T + tk;           # T [C]; TK [K]
    Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
    iom0 = 19.924*S/(1000-1.005*S);
    ST = 0.14 * Cl/96.062  # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)


	#--------------------------------------------------------------
	#------------------ Ks ----------------------------------------
	#       Dickson (1990)
	#     
	#       Equilibrium constant for HSO4- = H+ + SO4--
	#
	#       K_S  = [H+]free [SO4--] / [HSO4-]
	#       pH-scale: free scale !!!
	#
	#       the term log(1-0.001005*S) converts from
	#       mol/kg H2O to mol/kg soln
	
	
	tmp1 = -4276.1 / TK + 141.328 -23.093*log(TK);
	tmp2 = +(-13856 / TK + 324.57 - 47.986 * log(TK))*sqrt(iom0);
	tmp3 = +(35474 / TK - 771.54 + 114.723 * log(TK))*iom0;
	tmp4 = -2698 / TK *(sqrt(iom0)*iom0) + 1776 / TK *iom0 *iom0;
	                                       
	
	lnKs = tmp1 + tmp2 + tmp3 + tmp4 + log(1-0.001005*S);
		
	Ks_d = exp(lnKs);

	#--------------------------------------------------------------
	#------------------ Ks ----------------------------------------
	#       Khoo et al. 1977
	#     
	#       Equilibrium constant for HSO4- = H+ + SO4--
	#
	#       K_S  = [H+]free [SO4--] / [HSO4-]
	#       pH-scale: free scale !!!
	#        
	#      correct for T range : 5 - 40°C
	#      correct for S range : 20 - 45

    I35 <- 0.7227
    I <- I35*(27.570*S)/(1000-1.0016*S)   #formal ionic strength
    logbeta <- 647.59/TK - 6.3451 + 0.019085*TK - 0.5208*I^(0.5)
    beta <- 10^(logbeta)
    Ks_k <- 1/beta


    #-------------------------------------------------------------------
    #--------------- choice between the formulations -------------------
    is_k <- (ks=='k')    # everything that is not "k" is "d" by default
    Ks <- Ks_d
    Ks[is_k] <- Ks_k[is_k]

    method <- rep(NA, nK)
    method[!is_k] <- "Dickson (1990)"
    method[ is_k] <- "Khoo et al. (1977)"
    
    # ------------------- Pressure effect --------------------------------

    if (any(P != 0))
        Ks <- Pcorrect(Kvalue=Ks, Ktype="Ks", T=T, S=S, P=P, pHscale="F")

    ##------------Warnings

    if (any (is_k & (T<5 | T>40 | S<20 | S>45)))
        {warning("S and/or T is outside the range of validity of the formulation chosen for Ks.")}
    if (any (T>45 | S>45 | T<0 | S<5))
        warning("S and/or T is outside the range of validity of the formulations available for Ks in seacarb.")

    ##------------Attributes

    attr(Ks,"unit") = "mol/kg-soln"	
    attr(Ks,"pH scale") = "free scale"
    attr(Ks, "method") = method
    return(Ks)
}

