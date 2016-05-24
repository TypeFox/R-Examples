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


"bjerrum" <- function(K1=K1(),K2=NULL,K3=NULL,      # the dissociation constants
                    phmin=2, phmax=12, by=0.1, # the pH range and increment
                    conc=1,                    # the concentration, mol / kg
                    type="l",col="black",    # overruled default plotting options
                    ylab="Concentration (mol/kg)",
                    add=FALSE,               # false:start new, true: add to current
                    ...)                     # plotting options passed to matplot

{

# creates a bjerrum plot
 if (is.null(K1)) return()

 pH  <- seq (phmin,phmax,by=by)   # gradient of pH
 res <- speciation(K1, K2, K3, pH, conc)

 if (!add) matplot(pH, as.data.frame(res), ylab=ylab, type=type, col=col, ...)
 else matlines(pH, as.data.frame(res), type=type, col=col, ...)
}
