# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
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
"Kspa" <-
function(S=35,T=25,P=0){

nK <- max(length(S), length(T), length(P))

    ##-------- Creation de vecteur pour toutes les entrees (si vectorielles)
    
    if(length(S)!=nK){S <- rep(S[1], nK)}
    if(length(T)!=nK){T <- rep(T[1], nK)}
    if(length(P)!=nK){P <- rep(P[1], nK)}
    
    	
    #---- issues de equic----
    tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
    TC = T + tk;           # TC [C]; T[K]
    
    # --------------------- Kspa (aragonite) ----------------------------
    #
    # apparent solubility product of aragonite
    #
    #  Kspa = [Ca2+]T [CO32-]T
    #
    #  where $[]_T$ refers to the equilibrium total 
    # (free + complexed) ion concentration.
    #
    #  Mucci 1983 mol/kg-soln
    
    tmp1 = -171.945-0.077993*TC+2903.293/TC+71.595*log10(TC);
    tmp2 = +(-0.068393+0.0017276*TC+88.135/TC)*sqrt(S);
    tmp3 = -0.10018*S+0.0059415*S^1.5;
    log10Kspa = tmp1 + tmp2 + tmp3;
    
    Kspa = 10^(log10Kspa);

    # ----------------- Pressure Correction ------------------
    	
    Kspa <- Pcorrect(Kvalue=Kspa, Ktype="Kspa", T=T, S=S, P=P)
    
		
    ##------------Warnings

    
    if (any (T>40 | S>44 | T<5 | S<5))
        warning("S and/or T is outside the range of validity of the formulation available for Kspa in seacarb.")
		
    attr(Kspa,"unit") = "mol/kg"
    return(Kspa)
}
