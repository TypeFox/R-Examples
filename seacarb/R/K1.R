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
"K1" <-
function(S=35,T=25,P=0,k1k2='x',pHscale="T",kSWS2scale=0,ktotal2SWS_P0=0)
{

    nK <- max(length(S), length(T), length(P), length(k1k2), length(pHscale), length(kSWS2scale) ,length(ktotal2SWS_P0))

    ##-------- Create vectors for all input (if vectorial)

    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(k1k2)!=nK){k1k2 <- rep(k1k2[1], nK)}
    if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}

    ##---------- pHsc : this vector is not the actual pHscale because it can change during processing
    pHsc <- rep(NA,nK)
    pHlabel <- rep(NA,nK)
    K1 <- rep(NA, nK)

    ##----------Check the validity of the method regarding the T/S range
    is_x <- k1k2 == 'x'
    is_outrange <- T>35 | T<2 | S<19 | S>43
    k1k2[is_x] <- 'l'  ## luecker by default
    k1k2[is_x & is_outrange] <- "m10"  # Millero 2010 if outrange

    #-------Constants----------------

    #---- issues of equic----
    tk <- 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK <- T + tk;           # TC [C]; T[K]
    
    
    # --------------------- K1 ---------------------------------------
    #   first acidity constant:
    #   [H^+] [HCO_3^-] / [CO2] = K_1
    #
    #     Mehrbach et al (1973) refit by Lueker et al. (2000).
    #
    #(Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
    #   Dickson, Sabine and Christian , 2007, Chapter 5, p. 13)
    #
    #   pH-scale: 'total'. mol/kg-soln
    is_l <- k1k2 == "l"
    logK1lue <- (-3633.86)/TK[is_l] + 61.2172 - 9.67770*log(TK[is_l]) + 0.011555*S[is_l] - 0.0001152*S[is_l]*S[is_l]
    K1[is_l]<- 10^logK1lue
    pHsc[is_l] <- "T"


    # --------------------- K1 ---------------------------------------
    #   first acidity constant:
    #   [H^+] [HCO_3^-] / [CO2] = K_1
    #
    #   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
    #   pH-scale: 'total'. mol/kg-soln
    is_r <- k1k2 == "r"
    tmp1  <- 2.83655 - 2307.1266 / TK[is_r] - 1.5529413 * log(TK[is_r])
    tmp2  <- - (0.20760841 + 4.0484 / TK[is_r]) * sqrt(S[is_r])
    tmp3 <- 0.08468345 * S[is_r] - 0.00654208 * S[is_r] * sqrt(S[is_r])   
    tmp4 <- log(1 - 0.001005 * S[is_r])
    lnK1roy <- tmp1 + tmp2 + tmp3 + tmp4
    K1[is_r] <- exp(lnK1roy)
    pHsc[is_r] <- "T"

        
    # --------------------- K1 ---------------------------------------
    #   first acidity constant:
    #   [H^+] [HCO_3^-] / [CO2] = K_1
    #
    #   Millero et al. 2006 Marine Chemistry
    #   pH-scale: 'SWS scale'. mol/kg-soln
    is_m06 <- k1k2 == "m06"
    pK1o <- 6320.813/TK[is_m06] + 19.568224*log(TK[is_m06]) -126.34048
    A1 <- 13.4191*S[is_m06]^(0.5) + 0.0331*S[is_m06] - (5.33e-5)*S[is_m06]^2
    B1 <- -530.123*S[is_m06]^(0.5) - 6.103*S[is_m06]
    C1 <- -2.06950*S[is_m06]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_m06] + C1*log(TK[is_m06])
    K1[is_m06] <- 10^(-pK1)
    pHsc[is_m06] <- "SWS"


    # --------------------- K1 ---------------------------------------
    #   first acidity constant:
    #   [H^+] [HCO_3^-] / [CO2] = K_1
    #
    #   Millero 2010 Marine and Fresh water research

    is_m10 <- k1k2 == "m10"
    pK1o <- 6320.813/TK[is_m10] + 19.568224*log(TK[is_m10]) -126.34048

    #   pH-scale: 'SWS scale'. mol/kg-soln
    is_m10_SWS <- is_m10 & (pHscale=="SWS" | P>0)
    A1 <- 13.4038*S[is_m10_SWS]^(0.5) + 0.03206*S[is_m10_SWS] - (5.242e-5)*S[is_m10_SWS]^2
    B1 <- -530.659*S[is_m10_SWS]^(0.5) - 5.8210*S[is_m10_SWS]
    C1 <- -2.0664*S[is_m10_SWS]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_m10_SWS] + C1*log(TK[is_m10_SWS])
    K1[is_m10_SWS] <- 10^(-pK1)   # K1 according to Millero et al. 2010 at Seawater scale
    pHsc[is_m10_SWS] <- "SWS" 

    #   pH-scale: 'Total scale'. mol/kg-soln
    is_m10_T <- is_m10 & (pHscale=="T" & P==0)
    A1 <- 13.4051*S[is_m10_T]^(0.5) + 0.03185*S[is_m10_T] - (5.218e-5)*S[is_m10_T]^2
    B1 <- -531.095*S[is_m10_T]^(0.5) - 5.7789*S[is_m10_T]
    C1 <- -2.0663*S[is_m10_T]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_m10_T] + C1*log(TK[is_m10_T])
    K1[is_m10_T] <- 10^(-pK1)   # K1 according to Millero et al. 2010 at Total scale 
    pHsc[is_m10_T] <- "T"

    #   pH-scale: 'Free scale'. mol/kg-soln
    is_m10_F <- is_m10 & (pHscale=="F" & P==0)
    A1 <- 5.09247*S[is_m10_F]^(0.5) + 0.05574*S[is_m10_F] - (9.279e-5)*S[is_m10_F]^2
    B1 <- -189.879*S[is_m10_F]^(0.5) - 11.3108*S[is_m10_F]
    C1 <- -0.8080*S[is_m10_F]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_m10_F] + C1*log(TK[is_m10_F])
    K1[is_m10_F] <- 10^(-pK1)  # K1 according to Millero et al. 2010 at Total scale
    pHsc[is_m10_F] <- "F"

    # --------------------- K1 ---------------------------------------
    #   first acidity constant:
    #   [H^+] [HCO_3^-] / [CO2] = K_1
    #
    #   Waters, Millero, Woosley (Mar. Chem., 165, 66-67, 2014)

    is_w14 <- k1k2 == "w14"
    pK1o <- 6320.813/TK[is_w14] + 19.568224*log(TK[is_w14]) -126.34048

    #   pH-scale: 'SWS scale'. mol/kg-soln
    is_w14_SWS <- is_w14 & (pHscale=="SWS" | P>0)
    A1 <- 13.409160*S[is_w14_SWS]^(0.5) + 0.031646*S[is_w14_SWS] - (5.1895e-5)*S[is_w14_SWS]^2
    B1 <- -531.3642*S[is_w14_SWS]^(0.5) - 5.713*S[is_w14_SWS]
    C1 <- -2.0669166*S[is_w14_SWS]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_w14_SWS] + C1*log(TK[is_w14_SWS])
    K1[is_w14_SWS] <- 10^(-pK1)   # K1 according to Millero et al. 2010 at Seawater scale
    pHsc[is_w14_SWS] <- "SWS" 

    #   pH-scale: 'Total scale'. mol/kg-soln
    is_w14_T <- is_w14 & (pHscale=="T" & P==0)
    A1 <- 13.568513*S[is_w14_T]^(0.5) + 0.031645*S[is_w14_T] - (5.3834e-5)*S[is_w14_T]^2
    B1 <- -539.2304*S[is_w14_T]^(0.5) - 5.635*S[is_w14_T]
    C1 <- -2.0901396*S[is_w14_T]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_w14_T] + C1*log(TK[is_w14_T])
    K1[is_w14_T] <- 10^(-pK1)   # K1 according to Millero et al. 2010 at Total scale 
    pHsc[is_w14_T] <- "T"

    #   pH-scale: 'Free scale'. mol/kg-soln
    is_w14_F <- is_w14 & (pHscale=="F" & P==0)
    A1 <- 5.592953*S[is_w14_F]^(0.5) + 0.028845*S[is_w14_F] - (6.388e-5)*S[is_w14_F]^2
    B1 <- -225.7489*S[is_w14_F]^(0.5) - 4.761*S[is_w14_F]
    C1 <- -0.8715109*S[is_w14_F]^(0.5)
    pK1 <- pK1o + A1 + B1/TK[is_w14_F] + C1*log(TK[is_w14_F])
    K1[is_w14_F] <- 10^(-pK1)  # K1 according to Millero et al. 2010 at Total scale
    pHsc[is_w14_F] <- "F"

    ##----------------- Conversion from total to SWS scale
    ##                  if pressure correction needed
    ##                  or pH scale conversion required anyway
    ##                 (in which case SWS may be an intermediate stage of conversion)
    convert <- (pHsc == "T") & ((P > 0) | (pHscale != pHsc))
    if (any (convert))
    {
        ##------------- Convert from total to SWS scale
        # if correction factor (from Total scale to seawater at P=0) not given
        if (missing(ktotal2SWS_P0))
        {
            # Compute it
            ktotal2SWS_P0  <- rep(1.0,nK)
            ktotal2SWS_P0 <- kconv(S=S[convert], T=T[convert], P=0)$ktotal2SWS
        }
        else
        {
            # Check its length
            if(length(ktotal2SWS_P0)!=nK) ktotal2SWS_P0 <- rep(ktotal2SWS_P0[1], nK)
            # Filter
            ktotal2SWS_P0 <- ktotal2SWS_P0[convert]
        }
        K1[convert] <- K1[convert] * ktotal2SWS_P0
        pHsc[convert] <- "SWS"

        # Just computed constant is on Free scale only if required scale is free scale
        # and if no pressure correction needed  
        # --> No need to convert from free to other scale
        # --> No need to determine conversion factor from free to SWS scale
    }

    # ------------------- Pressure effect --------------------------------
    i_press <- which (P > 0)
    if (length(i_press) > 0)
    {
        # Call Pcorrect() on SWS scale
	# issue (Orr): Why are the last 2 argupments set to one here?
        K1[i_press] <- Pcorrect(Kvalue=K1[i_press], Ktype="K1", T=T[i_press], 
             S=S[i_press], P=P[i_press], pHscale=pHsc[i_press], 1., 1.)
        #K1[i_press] <- Pcorrect(Kvalue=K1[i_press], Ktype="K1", T=T[i_press], 
        #    S=S[i_press], P=P[i_press], pHscale=pHsc[i_press])
    }


    ## --------------- Last conversion in the required pHscale ----------------

    convert <- (pHscale != pHsc)
    # At this stage, if any conversion is needed, it is from SWS scale
    # Which pH scales are required ?
    is_total <- convert & (pHscale=="T")
    is_free  <- convert & (pHscale=="F")

    # if any pH scale correction required (from SWS scale)
    if (any(convert))
    {
        # if pH scale correction factor not given
        if (missing(kSWS2scale))
        {
            # Compute it
            kSWS2scale <- rep(1.0,nK)
	    if (any(is_total)){ kSWS2scale[is_total] <- kconv(S=S[is_total], T=T[is_total], P=P[is_total])$kSWS2total }
            if (any(is_free)){  kSWS2scale[is_free]  <- kconv(S=S[is_free],  T=T[is_free],  P=P[is_free])$kSWS2free }
        }
        else
            # Check its length
            if (length(kSWS2scale)!=nK) kSWS2scale <- rep(kSWS2scale[1], nK)
        # Apply pH scale correction
        K1[convert] <- K1[convert] * kSWS2scale[convert]
    }

    # Return full name of pH scale
    is_total <- (pHscale=="T")
    is_free  <- (pHscale=="F")
    pHlabel[is_total] <- "total scale"
    pHlabel[is_free]  <- "free scale"
    pHlabel[!is_total & !is_free] <- "seawater scale"


    ##------------Warnings

    if (any(is_l & (T>35 | T<2 | S<19 | S>43))  || any(is_r & (T>45 | S<5 | S>45)) )
        warning("S and/or T is outside the range of validity of the formulation chosen for K1.")

    if (any(T>50 | S>50)) 
        warning("S and/or T is outside the range of validity of the formulations available for K1 in seacarb.")

    ##---------------Attributes
    method <- rep(NA, nK)
    method[is_m06] <- "Millero et al. (2006)"
    method[is_m10] <- "Millero (2010)"
    method[is_w14] <- "Waters et al. (2014)"
    method[is_r]   <- "Roy et al. (1993)"
    method[! (is_m06 | is_m10 | is_w14 | is_r) ] <- "Luecker et al. (2000)"

    attr(K1,"unit") <- "mol/kg-soln"
    attr(K1,"pH scale") <- pHlabel
    attr(K1, "method") <- method
    return(K1)
}
