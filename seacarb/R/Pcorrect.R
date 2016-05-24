# Copyright (C) 2010 Jean-Pierre Gattuso and Heloise Lavigne
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"Pcorrect" <-
function(Kvalue, Ktype, T=25, S=35, P=0, pHscale="T", kconv2ScaleP0=0, kconv2Scale=0)
{
  #Pcoeffs <- get("Pcoeffs", envir  = environment()) # added to silence the CRAN note "Pcorrect: no visible binding for global variable ‘Pcoeffs’"
  nK <- max(length(Kvalue), length(Ktype), length(P), length(T), length(pHscale), length(S), length(kconv2ScaleP0), length(kconv2Scale))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

    if(length(Kvalue)!=nK){Kvalue <- rep(Kvalue[1], nK)}
    if(length(Ktype)!=nK){Ktype <- rep(Ktype[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}

    # If conversion factor at zero pressure is given 
    if (! missing(kconv2ScaleP0))
        # Check length of given conversion factor vector
        if (length(kconv2ScaleP0) != nK)
            kconv2ScaleP0 <- rep(kconv2ScaleP0[1], nK)
            
    # If conversion factor at given pressure is given 
    if (! missing(kconv2Scale))
        if (length(kconv2Scale) != nK)
            kconv2Scale <- rep(kconv2Scale[1], nK)
            

    # Constants 
    TK <- T + 273.15   # Temperature in Kelvin
    R = 83.14472;             # mol bar deg-1
    	
    ## loading table with coefficients
   # if(!exists("Pcoeffs", where = .GlobalEnv)) data(Pcoeffs)

    # ------------------- Pressure effect --------------------------------
    
	## 
	##  here, pHscale given in argument is the pH scale of Kvalue also given in argument
	##  if pHscale is different from seawater scale, it is necessary to convert K in seawater scale
	##  before applying the pressure correction, except with Kf, Kspc or Kspa.
	##  whith Kf, Kspc and Kspa pressure correction applied on free scale


    ##
    ## 1)
    ## Process  Ktype %in% c("K1", "K2", "K1p", "K2p", "K3p", "Kb", "Khs", "Kn", "Ksi", "Kw")
    ## 

	# Indices of Kvalue elements where pressure correction to apply on seawater scale
    i_SWscale <- which (P > 0.0 & Ktype %in% c("K1", "K2", "K1p", "K2p", "K3p", "Kb", "Khs", "Kn", "Ksi", "Kw"))     
    # if there is any such element
    if (length(i_SWscale) > 0)
    {
        # If conversion factor at zero pressure is not given 
        if (missing(kconv2ScaleP0))
        {
            ## ---------------- Compute this factor where pressure correction to apply on seawater scale
            
            # Initialise vector of conversion factor
            n_SWscale <- length(i_SWscale)
            conv <- rep(1., n_SWscale)
            
            # Indices in subvector Kvalue[i_SWscale] relevant to pHscale "T" 
            i_from_T <- which (pHscale[i_SWscale] == "T")
            if (length(i_from_T) > 0)
            {
                # Indices in initial vector Kvalue candidate to conversion from T to SWS
                i_from_T_SWS <- i_SWscale[i_from_T]             
                # compute conversion factor from Total pH scale to SW pH scale at zero pressure
                conv[i_from_T] <- kconv(S=S[i_from_T_SWS], T=T[i_from_T_SWS], P=0)$ktotal2SWS        
            }

            # Indices in subvector Kvalue[i_SWscale] relevant to pHscale "F" 
            i_from_F <- which (pHscale[i_SWscale] == "F")
            if (length(i_from_F) > 0)
            {
                # Indices in initial vector Kvalue candidate to conversion from F to SWS
                i_from_F_SWS <- i_SWscale[i_from_F] 
                # compute conversion factor from Free pH scale to SW pH scale at zero pressure
                conv[i_from_F] <- kconv(S=S[i_from_F_SWS], T=T[i_from_F_SWS], P=0)$free2SWS
            }        
        }
        else
        {
            # extract subset of relevant conversion factors
            conv <- kconv2ScaleP0[i_SWscale] 
        }
        
        # Apply conversion to SWS at zero pressure 
        Kvalue[i_SWscale] <- Kvalue[i_SWscale] * conv
    
        ## ------------------ Pressure correction
            
        l <- match(Ktype[i_SWscale], Pcoeffs$K)
        deltav  <-  Pcoeffs$a0[l] + Pcoeffs$a1[l] *T[i_SWscale] + Pcoeffs$a2[l] *T[i_SWscale]*T[i_SWscale]
        deltak  <-  Pcoeffs$b0[l]  + Pcoeffs$b1[l] *T[i_SWscale] + Pcoeffs$b2[l] *T[i_SWscale]*T[i_SWscale]
        lnkpok0 <-  -(deltav /(R*TK[i_SWscale]))*P[i_SWscale] + (0.5*deltak /(R*TK[i_SWscale]))*P[i_SWscale]*P[i_SWscale];
        Kvalue[i_SWscale] = Kvalue[i_SWscale]*exp(lnkpok0);
    
    
        # If conversion factor at given pressure is not given 
        if (missing(kconv2Scale))
        {
            ## ---------------  Compute this factor where pressure correction is applied on seawater scale
            ##
    
            # Initialise vector of conversion factor
            n_SWscale <- length(i_SWscale)
            conv <- rep(1., n_SWscale)
            
            # Indices in subvector Kvalue[i_SWscale] relevant to pHscale "T" 
            i_from_T <- which (pHscale[i_SWscale] == "T")
            if (length(i_from_T) > 0)
            {
                # Indices in initial vector Kvalue candidate to conversion from T to SWS
                i_from_T_SWS <- i_SWscale[i_from_T]
                # compute conversion factor from Total pH scale to SW pH scale at given pressure
                conv[i_from_T] <- kconv(S=S[i_from_T_SWS], T=T[i_from_T_SWS], P=P[i_from_T_SWS])$ktotal2SWS        
            }
    
            # Indices in subvector Kvalue[i_SWscale] relevant to pHscale "F" 
            i_from_F <- which (pHscale[i_SWscale] == "F")
            if (length(i_from_F) > 0)
            {
                # Indices in initial vector Kvalue candidate to conversion from F to SWS
                i_from_F_SWS <- i_SWscale[i_from_F] 
                # compute conversion factor from Free pH scale to SW pH scale at given pressure
                conv[i_from_F] <- kconv(S=S[i_from_F_SWS], T=T[i_from_F_SWS], P=P[i_from_F_SWS])$free2SWS
            }        
        }
        else
        {
            # extract subset of relevant conversion factors
            conv <- kconv2Scale[i_SWscale] 
        }
        
        # Apply conversion from SWS to given pH scale 
        Kvalue[i_SWscale] <- Kvalue[i_SWscale] / conv
    }

    ##
    ## 2) 
    ## Process  Ktype %in% c("Ks", "Kspa", "Kspc", "Kf")
    ## 

	# Indices of Kvalue elements where pressure correction to apply on free pH scale
    i_Fscale <- which (P > 0.0 & Ktype %in% c("Ks", "Kspa", "Kspc", "Kf"))     
    # if there is any such element
    if (length(i_Fscale) > 0)
    {
        # There is no pH scale choice for 
        #    Ks -> "free scale"
        #    Kspa and Kspc do not need pH scale
          
        # For Kf, there is pH scale choice
        #   must convert from given scale to free pH scale
    
        # Indices in subvector Kvalue[i_Fscale] relevant to "Kf" 
        i_Kf <- which (Ktype[i_Fscale] == "Kf")
        # Indices in initial vector Ktype candidate to conversion to free scale
        i_Kf_to_free <- i_Fscale[i_Kf] 
        
        # if there is any such candidate
        if (length(i_Kf_to_free) > 0)
        {                 
            # If conversion factor at zero pressure is not given 
            if (missing(kconv2ScaleP0))
            {
                ## ---------------- Compute this factor where applying to Kf where P > 0.0
                
                # Initialise vector of conversion factor
                n_Kf_to_free <- length(i_Kf_to_free)
                conv <- rep(1., n_Kf_to_free)
                
                # Indices in subvector Kvalue[i_Kf_to_free] relevant to Kf on pHscale "T" 
                i_from_T <- which (pHscale[i_Kf_to_free] == "T")
                # if there is any such element in subvector
                if (length(i_from_T) > 0)
                {
                    # Indices in full vector Kvalue
                    i_from_T_full <- i_Kf_to_free[i_from_T]                
                  
                    # compute conversion factor from Total pH scale to Free pH scale at zero pressure
                    Cl = S[i_from_T_full] / 1.80655;  # Cl = chlorinity; S = salinity (per mille)
                    ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
                    FT = 6.7e-5 * Cl/18.9984    # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
                    conv[i_from_T] <- 1./(1+ST/Ks(S=S[i_from_T_full], T=T[i_from_T_full], P=0))
                }
              
                # Indices in subvector Kvalue[i_Kf_to_free] relevant to Kf on pHscale "SWS" 
                i_from_SWS <- which (pHscale[i_SWscale] == "SWS")
                # if there is any such element in subvector
                if (length(i_from_SWS) > 0)
                {
                    # Indices in full vector Kvalue
                    i_from_SWS_full <- i_Kf_to_free[i_from_SWS]                
                  
                    # compute conversion factor from SW pH scale to Free pH scale at zero pressure
                    Cl = S[i_from_SWS_full] / 1.80655;  # Cl = chlorinity; S = salinity (per mille)
                    ST = 0.14 * Cl/96.062         # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
                    FT = 6.7e-5 * Cl/18.9984    # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
                    conv[i_from_SWS] <- 1. / (1 + ST/Ks(S=S[i_from_SWS_full], T=T[i_from_SWS_full], P=0)
                                      + FT/Kf(T=T[i_from_SWS_full], P=0, S=S[i_from_SWS_full], pHscale="F"))
                }
            }
            else
            {
                # extract subset of relevant conversion factors
                conv <- kconv2ScaleP0[i_Kf_to_free] 
            }
            
            # Apply conversion at zero pressure (for Kf only)
            Kvalue[i_Kf_to_free] <- Kvalue[i_Kf_to_free] * conv
        }
            
        ## ------------------ Pressure correction
        
        l <- match(Ktype[i_Fscale], Pcoeffs$K)
        deltav  <-  Pcoeffs$a0[l] + Pcoeffs$a1[l] *T[i_Fscale] + Pcoeffs$a2[l] *T[i_Fscale]*T[i_Fscale]
        deltak  <-  Pcoeffs$b0[l]  + Pcoeffs$b1[l] *T[i_Fscale] + Pcoeffs$b2[l] *T[i_Fscale]*T[i_Fscale]
        lnkpok0 <-  -(deltav /(R*TK[i_Fscale]))*P[i_Fscale] + (0.5*deltak /(R*TK[i_Fscale]))*P[i_Fscale]*P[i_Fscale];
        
        Kvalue[i_Fscale] = Kvalue[i_Fscale]*exp(lnkpok0);
    }
    
    # Update pH scale information
    pHscale [ (P > 0.0) & (Ktype %in% c("Ks", "Kf")) ] <- "F" 

    phs <- rep(NA, nK)
    phs [pHscale == "SWS"] <- "Seawater scale"
    phs [pHscale == "T"]   <- "Total scale"
    phs [pHscale == "F"]   <- "Free scale"
    phs [Ktype %in% c("Kspa", "Kspc")] <- NA

    attr(Kvalue, "pH scale") <- phs
    attr(Kvalue, "unit") <- "mol/kg-soln"
    return(Kvalue)
}

