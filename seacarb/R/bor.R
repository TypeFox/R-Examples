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
"bor" <-
function(S=35, b="u74"){  #two formulations possible Lee et al., (2010) : "l10" or Uppstrom (1974) : "u74"

n <- max(c(length(S), length(b)))
if(length(S)!=n){ S <- rep(S[1],n)}
if(length(b)!=n){ b <- rep(b[1],n)}
bor <- rep(NA, n)
method <- rep(NA, n)

for(i in 1:n){
if(b[i]=="l10"){
bor[i] <- 0.1336*S[i] #boron in mg/kg
method[i] <- "Lee et al. (2010)"
}
if(b[i]=="u74"){
bor[i] <- 0.1284*S[i] #boron in mg/kg
method[i] <- "Uppstrom (1974)"
}
}
# conversion from mg/kg to mol/kg
bor <- bor*10^(-3)/10.811
attr(bor,"unit") <- "mol/kg"
attr(bor, "method") <- method

return(bor)
}
