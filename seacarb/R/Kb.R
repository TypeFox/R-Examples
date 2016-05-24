# Copyright (C) 2008 Jean-Pierre Gattuso and Jean-Marie Epitalon and Heloise Lavigne and Aurelien Proye
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

"Kb" <-
function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0){

    nK <- max(length(S), length(T), length(P), length(pHscale), length(kSWS2scale) ,length(ktotal2SWS_P0))

    ##-------- Create vectors for all input (if vectorial)

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}
            
    #-------Constants----------------

    #---- issues of equic----
    tk <- 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK <- T + tk;           # T [C]; TK [K]

            
    #---------------------------------------------------------------------
    # --------------------- Kb  --------------------------------------------
    #  Kbor = [H+][B(OH)4-]/[B(OH)3]
    #
    #   (Dickson, 1990 in Guide to Best Practices in Ocean CO2 Measurements 2007)
    #   pH-scale: 'total'. mol/kg-soln
    
    
    tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S^(3/2)-0.0996*S*S);
    tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
    tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S)*log(TK);
    
    lnKb = tmp1 / TK + tmp2 + tmp3 + 0.053105*sqrt(S)*TK;
    Kb <- exp(lnKb)

    ## ---- Conversion from Total scale to seawater scale before pressure corrections
    
    # if correction factor (from Total scale to seawater at P=0) not given
    if (missing(ktotal2SWS_P0))
    {
        # Compute it
        ktotal2SWS_P0 <- kconv(S=S, T=T, P=0)$ktotal2SWS
    }
    else
        # Check its length
        if(length(ktotal2SWS_P0)!=nK) ktotal2SWS_P0 <- rep(ktotal2SWS_P0[1], nK)
    Kb <- Kb * ktotal2SWS_P0

    # -------------------Correct for Pressure effect, if needed ---------------------------
    if (any (P != 0))
        Kb <- Pcorrect(Kvalue=Kb, Ktype="Kb", T=T, S=S, P=P, pHscale="SWS", 1., 1.)

    ###----------------pH scale corrections

    # Which pH scales are required ?
    is_total <- pHscale=="T"
    is_free   <- pHscale=="F"

    # if any pH scale correction required (from SWS scale)
    if (any(is_total) || any(is_free))
    {
        # if pH scale correction factor not given
        if (missing(kSWS2scale))
        {
            # Compute it
            kSWS2scale <- rep(1.0,nK)
            if (any(is_total)) {kSWS2scale[is_total] <- kconv(S=S[is_total], T=T[is_total], P=P[is_total])$kSWS2total}
            if (any(is_free))  {kSWS2scale[is_free]  <- kconv(S=S[is_free], T=T[is_free], P=P[is_free])$kSWS2free}
        }
        else
            # Check its length
            if(length(kSWS2scale)!=nK){kSWS2scale <- rep(kSWS2scale[1], nK)}
        # Apply pH scale correction
        Kb <- Kb*kSWS2scale
    }

    # Return full name of pH scale
    pHsc <- rep(NA,nK)
    pHsc[is_total] <- "total scale"
    pHsc[is_free]  <- "free scale"
    pHsc[!is_total & !is_free] <- "seawater scale"

    
    ##------------Warnings
    if (any (T>45 | S>45 | T<0 | S<5) ) {warning("S and/or T is outside the range of validity of the formulation available for Kb in seacarb.")}


    attr(Kb,"unit")     = "mol/kg-soln"
    attr(Kb,"pH scale") = pHsc
    return(Kb)
}
