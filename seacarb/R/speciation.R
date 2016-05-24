# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


"speciation" <- function(K1=K1(),K2=NULL,K3=NULL, # the dissociation constants
                         pH=8,               # the pH at which speciation is wanted 
                         conc=1)             # the concentration, mol/kg

{
univalent <- function(K1,H)
{ # concentration for univalent species (e.g. NH3,NH4+)
 C1<- H/(K1+H) # modificated by He\'loi\¨se Lavigne (C1=K1/(K1+H)-> C1=  H/(K1+H)) 2008
 C2<-1-C1
 return(list(C1=C1*conc,C2=C2*conc))
}

bivalent <- function(K1,K2,H)
{ # concentration for bivalent species(e.g. CO2,HCO3-,CO3--)
 den <-H*H + H*K1 + K1*K2
 C1  <-H*H/den
 C2  <-H*K1/den
 C3  <- 1-C1-C2
 return(list(C1=C1*conc,C2=C2*conc,C3=C3*conc))
}

trivalent <- function(K1,K2,K3,H)
{ # concentration for trivalent species(e.g. H3PO4,H2PO4-,HPO4--,PO4---)
 den <-H*H*H + H*H*K1 + H*K1*K2 + K1*K2*K3
 C1  <-H*H*H/den
 C2  <-H*H*K1/den
 C3  <-H*K1*K2/den
 C4  <- 1-C1-C2-C3
 return(list(C1=C1*conc,C2=C2*conc,C3=C3*conc,C4=C4*conc))
}

# estimates the speciation of the various ionic forms of a molecule as a function of pH 
 if (is.null(K1)) return()
 H   <- 10^(-pH)                  # proton concentration
 if (is.null(K2))      res<-univalent(as.double(K1),H)        
 else if (is.null(K3)) res<-bivalent (as.double(K1),as.double(K2),H)     
 else                  res<-trivalent(as.double(K1),as.double(K2),as.double(K3),H)  
 return(res)
}
