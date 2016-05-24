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
"K1p" <-
function(S=35,T=25,P=0,pHscale="T",kSWS2scale=0){

    nK <- max(length(S), length(T), length(P), length(pHscale), length(kSWS2scale))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}


    #-------Constantes----------------

    #---- issues de equic----
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TC = T + tk;           # TC [C]; T[K]


    # --------------------- Phosphoric acid ---------------------
    #
    #   (DOE, 1994)  (Dickson and Goyet): pH_T, mol/(kg-soln)
    #   Ch.5 p. 16
    #	*** J. Orr (15 Jan 2013): Formulation changed to be on the SWS scale (without later conversion)
    
    #lnK1P = -4576.752 / TC + 115.525 - 18.453*log(TC) + (-106.736 / TC + 0.69171) * sqrt(S) + (-0.65643 / TC - 0.01844) * S;

    # From J. C. Orr on 15 Jan 2013:
    # The formulation above was a modified version of Millero (1995) where Dickson et al. (2007) subtracted 0.015
    # from Millero's original constant (115.54) to give 115.525 (the 2nd term above). BUT Dickson's reason for that 
    # operation was to "convert--approximately--from theSWS pH scale (including HF) used by Millero (1995) to the 'total' 
    # scale ...". 
    # This subtraction of 0.015 to switch from the SWS to Total scale is not good for 2 reasons:
    # (1) The 0.015 value is inexact (not constant), e.g., it is 0.022 at T=25, S=35, P=0;
    # (2) It makes no sense to switch to the Total scale when just below you switch back to the SWS scale.
    # The best solution is to reestablish the original equation (SWS scale) and delete the subsequent scale conversion.

    # now the original formulation: Millero (1995)
    lnK1P = -4576.752 / TC + 115.54  - 18.453*log(TC) + (-106.736 / TC + 0.69171) * sqrt(S) + (-0.65643 / TC - 0.01844) * S;
    
    K1P = exp(lnK1P);


    # ---- Conversion from Total scale to seawater scale before pressure corrections
    #      *** JC Orr: This is no longer necessary: with original formulation (Millero, 1995), K1P is on "seawater scale"!

    #factor <- kconv(S=S, T=T, P=rep(0,nK))$ktotal2SWS
    #K1P <- K1P * factor

    # ----------------- If needed, Pressure Correction ------------------
            
    if (any (P != 0)) 
        K1P <- Pcorrect(Kvalue=K1P, Ktype="K1p", T=T, S=S, P=P, pHscale="SWS", 1., 1.)

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
            if (any(is_total))
                kSWS2scale[is_total] <- kconv(S=S[is_total], T=T[is_total], P=P[is_total])$kSWS2total
            if (any(is_free))
                kSWS2scale[is_free]  <- kconv(S=S[is_free], T=T[is_free], P=P[is_free])$kSWS2free
        }
        else
            # Check its length
            if(length(kSWS2scale)!=nK){kSWS2scale <- rep(kSWS2scale[1], nK)}
        # Apply pH scale correction
        K1P <- K1P*kSWS2scale
    }

    # Return full name of pH scale
    pHsc <- rep(NA,nK)
    pHsc[is_total] <- "total scale"
    pHsc[is_free]  <- "free scale"
    pHsc[!is_total & !is_free] <- "seawater scale"


    ##------------Warnings

    if (any (T>45 | S>45 | T<0) ) {warning("S and/or T is outside the range of validity of the formulation available for K1p in seacarb.")}


    attr(K1P,"unit")     = "mol/kg-soln"
    attr(K1P,"pH scale") = pHsc
    return(K1P)
}
