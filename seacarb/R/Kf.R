# Copyright (C) 2008 Jean-Pierre Gattuso and Jean-Marie Epitalon and Heloise Lavigne and Aurelien Proye
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
#
"Kf" <-
function(S=35,T=25,P=0,kf='x',pHscale="T",Ks_p0=0,Ks_p=0){

    nK <- max(length(S), length(T), length(P), length(kf), length(pHscale), length(Ks_p0), length(Ks_p))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(kf)!=nK){kf <- rep(kf[1], nK)}
    if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}
    
##----------Make kf a global variable to facilitate pH scale changes in multiple routines

    if (missing(kf)) {
        if (exists("kfg",envir=parent.frame()) )
           {kf <- get("kfg", envir = parent.frame()) }
	else    
           {assign("kfg", "x", envir = parent.frame()) }
	}
    else
        {assign("kfg", kf, envir = parent.frame()) } 

    ##----------Check the validity of the method regarding the T/S range
    is_x <- kf=='x'
    is_outrange <- T>33 | T<10 | S<10 | S>40
    kf[is_x] <- "pf"  ## Perez and Fraga by default
    kf [is_x & is_outrange] <- "dg"
    
    #-------Constantes----------------

    #---- issues de equic----
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK = T + tk;           # TC [C]; T[K]
    Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
    iom0 = 19.924*S/(1000-1.005*S);
    ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
    FT = 6.7e-5 * Cl/18.9984    # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)

    #---------------------------------------------------------------------
    #---------------------- Kf Perez and Fraga ---------------------------
    #  Kf = [H+][F-]/[HF]  
    #
    #   Perez and Fraga, 1987 in Guide to the Best Practices for Ocean CO2 Measurements
    #   Dickson, Sabine and Christian, 2007, Chapter 5, p. 14)
    #  
    #   pH-scale: 'total'   

    lnKfpf <- 874/TK - 9.68 + 0.111*S^(1/2)
    Kfpf <- exp(lnKfpf)

    # --------------- Conversion to free scale for pressure corrections
    
    # if Ks at zero pressure is not given
    if (missing(Ks_p0))
        # Ks at zero pressure NOT given --> compute it
    	Ks = Ks(S=S, T=T, P=rep(0,nK))                 # on free pH scale
    else
    {
    	Ks <- Ks_p0
    	if (length(Ks)!=nK) Ks <- rep(Ks[1], nK)
    }
    
    Cl = S / 1.80655;             # Cl = chlorinity; S = salinity (per mille)
    ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
    total2free  = 1/(1+ST/Ks)      # Kfree = Ktotal*total2free
    total2free <- as.numeric(total2free)	

    factor <- total2free

    Kfpf <- Kfpf * factor

    #---------------------------------------------------------------------
    # --------------------- Kf Dickson and Goyet -------------------------
    #  Kf = [H+][F-]/[HF]  
    #
    #   (Dickson and Riley, 1979 in Dickson and Goyet, 
    #   1994, Chapter 5, p. 14)
    #   pH-scale: 'free' (require to convert in total scale after pressure corrections 

    lnKfdg = 1590.2/TK - 12.641 + 1.525*sqrt(iom0) + log(1-0.001005*S);

    Kfdg <- exp(lnKfdg)

    # ---------- Choice between methods (Perez and Dickson) ----------

    is_pf <- (kf!='dg')    # everything that is not "dg" is "pf" by default
    Kf <- Kfpf
    Kf[!is_pf] <- Kfdg[!is_pf]


    # ------------------- Pressure effect --------------------------------
    if (any (P != 0)) Kf <- Pcorrect(Kvalue=Kf, Ktype="Kf", T=T, S=S, P=P, pHscale="F") 
    # All the Kf constants are on free scale, whatever the method or the pressure

    # Which pH scales are required ?
    is_total <- pHscale=="T"
    is_SWS   <- pHscale=="SWS"

    # pH scale correction factor
    factor <- rep(1.0,nK)

    # if any pH scale corrections required (from free scale)
    if (any(is_total) || any(is_SWS))
    {
        # If Ks at given pressure is not given 
        if (missing(Ks_p))
            # compute Ks
            Ks = Ks(S=S, T=T, P=P)                 # on free pH scale
        else
        {
        	Ks <- Ks_p
        	if (length(Ks)!=nK) Ks <- rep(Ks[1], nK)
        }
        # Compute pH scale correction factor
        factor[is_total] <- 1 + ST[is_total]/Ks[is_total]
        factor[is_SWS]   <- 1 + ST[is_SWS]/Ks[is_SWS] + FT[is_SWS]/Kf[is_SWS]
    }
    # Perform pH scale correction
    Kf <- Kf*factor

    # Return full name of pH scale
    pHsc <- rep(NA,nK)
    pHsc[is_total] <- "total scale"
    pHsc[is_SWS]   <- "seawater scale"
    pHsc[!is_total & !is_SWS] <- "free scale"


    ##------------Warnings

    if (any (is_pf & (T>33 | T<9 | S<10 | S>40)) ) {warning("S and/or T is outside the range of validity of the formulation chosen for Kf.")}
    if (any( T>45 | S>45 )) {warning("S and/or T is outside the range of validity of the formulations available for Kf in seacarb.")}

    ##---------------Attributes
    method <- rep(NA,nK)
    method[ is_pf] <- "Perez and Fraga (1987)"
    method[!is_pf] <- "Dickson and Riley (1979 in Dickson and Goyet, 1994)"


    attr(Kf,"unit")     = "mol/kg-soln"
    attr(Kf,"pH scale") = pHsc
    attr(Kf, "method") = method
    return(Kf)
}
