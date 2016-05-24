# Copyright (C) 2008 Jean-Marie Epitalon and Jean-Pierre Gattuso and Heloise Lavigne and Aurelien Proye
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

buffer <- 
  function(flag, var1, var2, S=35, T=25, Patm=1, P=0, Pt=0, Sit=0, k1k2='x', kf='x', ks="d", pHscale="T", b="u74"){
    n <- max(length(flag), length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit), length(k1k2), length(kf), length(pHscale), length(ks), length(b))
    if(length(flag)!=n){flag <- rep(flag[1],n)}
    if(length(var1)!=n){var1 <- rep(var1[1],n)}
    if(length(var2)!=n){var2 <- rep(var2[1],n)}
    if(length(S)!=n){S <- rep(S[1],n)}
    if(length(T)!=n){T <- rep(T[1],n)}
    if(length(Patm)!=n){Patm <- rep(Patm[1],n)}
    if(length(P)!=n){P <- rep(P[1],n)}
    if(length(Pt)!=n){Pt <- rep(Pt[1],n)}
    if(length(Sit)!=n){Sit <- rep(Sit[1],n)}
    if(length(k1k2)!=n){k1k2 <- rep(k1k2[1],n)}
    if(length(kf)!=n){kf <- rep(kf[1],n)}
    if(length(ks)!=n){ks <- rep(ks[1],n)}
    if(length(pHscale)!=n){pHscale <- rep(pHscale[1],n)}
    if(length(b)!=n){b <- rep(b[1],n)}

    # if the concentrations of total silicate and total phosphate are NA
    # they are set at 0
    Sit[is.na(Sit)] <- 0
    Pt[is.na(Pt)] <- 0
    
    Carb <- carb(flag=flag, var1=var1, var2=var2, S=S, T=T, Patm=Patm, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, ks=ks, pHscale=pHscale, b=b)
    
	pH   <- Carb$pH
	h    <- 10^(-pH)
	CO2  <- Carb$CO2
	HCO3 <- Carb$HCO3
	CO3  <- Carb$CO3
	DIC  <- Carb$DIC
	ALK  <- Carb$ALK
	Oa   <- Carb$OmegaAragonite
	Oc   <- Carb$OmegaCalcite
    
    #-------Constants----------------
    
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TK = T + tk;           # TK [K]; T[C]
    
    Cl = S / 1.80655;            # Cl = chlorinity; S = salinity (per mille)
    ST = 0.14 * Cl/96.062        # (mol/kg) total sulfate  (Dickson et al., 2007, Table 2)
    FLUO = 6.7e-5 * Cl/18.9984   # (mol/kg) total fluoride (Dickson et al., 2007, Table 2)
    BOR = bor(S=S , b=b);        # (mol/kg) total boron

    #---------------------------------------------------------------------
    #--------------------- compute K's ----------------------------------
    #---------------------------------------------------------------------
    
    # Ks (free pH scale) at zero pressure and given pressure
    Ks_P0 <- Ks(S=S, T=T, P=0, ks=ks)
    Ks    <- Ks(S=S, T=T, P=P, ks=ks)
    
    # Kf on free pH scale
    Kff_P0 <- Kf(S=S, T=T, P=0, pHscale="F", kf=kf, Ks_P0, Ks)
    Kff <- Kf(S=S, T=T, P=P, pHscale="F", kf=kf, Ks_P0, Ks)
    # Kf on given pH scale
    Kf <- Kf(S=S, T=T, P=P, pHscale=pHscale, kf=kf, Ks_P0, Ks)
    
    # Conversion factor from total to SWS pH scale at zero pressure
    ktotal2SWS_P0 <- kconv(S=S,T=T,P=P,kf=kf,Ks=Ks_P0,Kff=Kff_P0)$ktotal2SWS

    # Conversion factor from SWS to chosen pH scale
    conv <- kconv(S=S,T=T,P=P,kf=kf,Ks=Ks,Kff=Kff)
    kSWS2chosen <- rep(1.,n)
    kSWS2chosen [pHscale == "T"] <- conv$kSWS2total [pHscale == "T"]
    kSWS2chosen [pHscale == "F"] <- conv$kSWS2free [pHscale == "F"]  

    K1 <- K1(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0)   
    K2 <- K2(S=S, T=T, P=P, pHscale=pHscale, k1k2=k1k2, kSWS2chosen, ktotal2SWS_P0)
    Kw <- Kw(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen)
    K0 <- K0(S=S, T=T, Patm=Patm, P=P)
    Kb <- Kb(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen, ktotal2SWS_P0)
    K1p <- K1p(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen)
    K2p <- K2p(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen)
    K3p <- K3p(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen)
    Ksi <- Ksi(S=S, T=T, P=P, pHscale=pHscale, kSWS2chosen)
    Kspa <- Kspa(S=S, T=T, P=P)
    Kspc <- Kspc(S=S, T=T, P=P)

    rho <- rho(S=S,T=T,P=P)

    # Compute potential K0 with S, potential temperature, and atmospheric pressure (usually 1 atm)
    K0pot <- K0(S=S, T=theta(S=S, T=T, P=P, Pref=0), Patm=Patm, P=0)
    
    #---------------------------------------------------------------------
    #--------------------    buffer effects    ---------------------------
    #---------------------------------------------------------------------
    
    DD=-((-Kb*BOR)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
    A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
    B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
    C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
    PhiD=-1/(h*log(10) * ( B+A+C ) );
    BetaD=-h*log(10)*DIC/CO2*B*PhiD;
    
    
    Q=(h+2*K2);
    V=(Kb*BOR)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;
    
    DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*BOR)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
    A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
    B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
    C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
    PhiB=-1/(h*log(10) * ( B+A+C ) );
    BetaB=-h*log(10)*DIC/CO2*B*PhiB;
    
    DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*BOR)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
    A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
    B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
    C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
    PhiC=-1/(h*log(10) * ( B+A+C ) );
    BetaC=-h*log(10)*DIC/CO2*B*PhiC;
    
    D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
    D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
    D=D1+2*D2;
    PhiH=1/ (h*log(10)* (D +(-Kb*BOR/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ; 
    
    Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));

    # In situ values - use K0 computed with in situ T & in situ P (atm + hydrostatic)
    # PiH=((-h/K0)*log(10)*Pi)*PhiH;
    # PiB=CO2/(K0*DIC)*BetaB;
    # PiD=CO2/(K0*DIC)*BetaD;
    # PiC=CO2/(K0*DIC)*BetaC;

    # Potential values - use K0 computed with potential T & surface atm P only
    PiH=((-h/K0pot)*log(10)*Pi)*PhiH;
    PiB=CO2/(K0pot*DIC)*BetaB;
    PiD=CO2/(K0pot*DIC)*BetaD;
    PiC=CO2/(K0pot*DIC)*BetaC;

    col <- c("PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
    res <- data.frame(PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH)
    names(res) <- col 
    
    return(res)
  }
