
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               PARAMETER ESTIMATION: 
#  'fGARCH'                S4: fGARCH Class representation   
#  garchFit                Fits GARCH and APARCH processes
################################################################################


test.garchFit.snorm <- 
    function()
{  
    # Conditional Densities:
    #   "norm", "snorm", "ged", "sged", "std", "sstd"
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # skewed normal GARCH(1, 1)
    model = list(omega = 1e-04, alpha = 0.1, beta = 0.8, skew = 0.9)
    spec = garchSpec(model, cond.dist = "snorm")
    x = garchSim(spec = spec, n = 250)
    
    # Fit:
    fit = garchFit( ~ garch(1,1), data = x, include.skew = TRUE, 
        cond.dist = "snorm", trace = FALSE)
    print(coef(fit))
    
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.garchFit.snorm.fixed <- 
    function()
{  
    # Conditional Densities:
    #   "norm", "snorm", "ged", "sged", "std", "sstd"
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # skewed normal GARCH(1, 1)
    model = list(omega = 1e-04, alpha = 0.1, beta = 0.8, skew = 0.9)
    spec = garchSpec(model, cond.dist = "snorm")
    x = garchSim(spec = spec, n = 250)
    
    # Fit: Skewed Normal GARCH(1, 1) with fixed skew ...
    fit = garchFit(~garch(1,1), data = x, skew = 0.9, 
        include.skew = FALSE, cond.dist = "snorm", trace = FALSE)
    print(coef(fit))
    
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.garchFit.ged <- 
    function()
{  
    # Conditional Densities:
    #   "norm", "snorm", "ged", "sged", "std", "sstd"
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")    
    
    # GED-GARCH(1, 1)
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8, shape = 2)
    spec = garchSpec(model, cond.dist = "ged")
    x = garchSim(spec = spec, n = 250)
    
    # Fit:
    fit = garchFit(~garch(1,1), data = x, 
        include.shape = TRUE, cond.dist = "ged", trace = FALSE)
    print(coef(fit))
    
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.garchFit.sged <- 
    function()
{  
    # Conditional Densities:
    #   "norm", "snorm", "ged", "sged", "std", "sstd"
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")    
    
    # Skewed GED-GARCH(1, 1)
    model = list(omega = 1e-05, alpha = 0.1, beta = 0.8, shape = 4, skew = 0.9)
    spec = garchSpec(model, cond.dist = "sged")
    x = garchSim(spec = spec, n = 250)
    
    # Fit
    fit = garchFit( ~ garch(1,1), data = x, 
        include.shape = TRUE, include.skew = TRUE, cond.dist = "sged", 
        trace = FALSE)
    print(coef(fit))
    
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.garchFit.std <- 
    function()
{  
    # Conditional Densities:
    #   "norm", "snorm", "ged", "sged", "std", "sstd"
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")    
    
    # Student t-GARCH(1, 1)
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8, shape = 5)
    spec = garchSpec(model, cond.dist = "std")
    x = garchSim(spec = spec, n = 250)
    
    # Fi
    fit = garchFit( ~ garch(1,1), data = x, 
        include.shape = TRUE, cond.dist = "std", trace = FALSE)
    print(coef(fit))
    
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.garchFit.sstd <- 
    function()
{  
    # Conditional Densities:
    #   "norm", "snorm", "ged", "sged", "std", "sstd"
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")    
    
    # Skewed Student t-GARCH(1, 1)
    model = list(omega = 1e-06, alpha = 0.1, beta = 0.8, 
        shape = 5, skew = 0.9)
    spec = garchSpec(model, cond.dist = "std")
    x = garchSim(spec = spec, n = 250)
    
    # Fit:
    fit = garchFit( ~ garch(1,1), data = x, 
        include.shape = TRUE, include.skew = TRUE, 
        cond.dist = "sstd", trace = FALSE)
    print(coef(fit))
    
    # Return Value:
    return()    
} 


################################################################################
    
