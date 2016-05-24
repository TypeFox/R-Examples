
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
# FUNCTION:               SIMULATION:
#  garchSim                Simulates a GARCH/APARCH process
################################################################################


test.garchSim.arch = 
function()
{   
    # Simulation of ARCH Models:
    
    # ARCH(1) - default omega and alpha
    spec = garchSpec(model = list())
    garchSim(n = 10, spec = spec)
    
    # ARCH(1) - default omega
    spec = garchSpec(model = list(alpha = 0.1))
    garchSim(n = 10, spec = spec)
    
    # ARCH(1)
    spec = garchSpec(model = list(omega = 1e-6, alpha = 0.1))
    garchSim(n = 10, spec = spec)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.garchSim.arma.arch = 
function()
{  
    # Simulation of ARMA-ARCH Models:
    
    # AR(1)-ARCH(1)
    spec = garchSpec(model = list(ar = 0.5)) 
    garchSim(n = 10, spec = spec)
    
    # AR([1,5])-ARCH(1)
    spec = garchSpec(model = list(ar = c(0.5, 0, 0, 0 ,0.1)))
    garchSim(n = 10, spec = spec)
    
    # ARMA(1,2)-ARCH(1)
    spec = garchSpec(model = list(ar = 0.5, ma = c(0.3,-0.3)))
    garchSim(n = 10, spec = spec)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.garchSim.dist.arch = 
function()
{     
    # Simulation of non-normal ARCH Models:
    
    # rsnorn-ARCH(2)
    spec = garchSpec(model = list(alpha = c(0.12, 0.04), dist = 2/3), 
        cond.dist = "snorm")
    garchSim(n = 10, spec = spec)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.garchSim.garch = 
function()
{    
    # Simulation of GARCH Models:
    
    # GARCH(1,1)
    spec = garchSpec()
    garchSim(n = 10, spec = spec)
    
    # GARCH(1,1)
    spec = garchSpec(model = list(alpha = 0.1, beta = 0.8))
    garchSim(n = 10, spec = spec)
    
    # GARCH(1,1)
    spec = garchSpec(model = list(omega = 1e-6, alpha = 0.1, beta = 0.8))
    garchSim(n = 10, spec = spec)
    
    # GARCH(1,2)
    spec = garchSpec(model = list(alpha = 0.1, beta = c(0.4, 0.4)))
    garchSim(n = 10, spec = spec)
    
    # GARCH(2,1)
    spec = garchSpec(model = list(alpha = c(0.12, 0.04), beta = 0.08))
    garchSim(n = 10, spec = spec)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.garchSim.dist.garch = 
function()
{     
    # Simulation of non-normal GARCH Models:    
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Normal - GARCH(1,1)   
    spec = garchSpec(model = list(), cond.dist = "norm")
    garchSim(n = 10, spec = spec)
    
    # Skewed Normal - GARCH(1,1)
    spec = garchSpec(model = list(parm = 2), cond.dist = "snorm")
    garchSim(n = 10, spec = spec)
    
    # GED - GARCH(1,1)
    spec = garchSpec(model = list(parm = 4), cond.dist = "ged")
    garchSim(n = 10, spec = spec)
    
    # Skewed GED - GARCH(1,1)
    spec = garchSpec(model = list(parm = c(4, 2)), cond.dist = "sged")
    garchSim(n = 10, spec = spec)
    
    # Normalized Student t - GARCH(1,1)
    spec = garchSpec(model = list(parm = 4), cond.dist = "std")
    garchSim(n = 10, spec = spec)
    
    # Skewed Normalized Student t - GARCH(1,1)
    spec = garchSpec(model = list(parm = c(4, 2)), cond.dist = "sstd")
    garchSim(n = 10, spec = spec)
 
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.garchSim.aparch = 
function()
{         
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # Taylor Schwert Normal GARCH(1,1)
    spec = garchSpec(list(alpha = 0.1, delta = 1))
    garchSim(n = 10, spec = spec)
    
    # Return Value:
    return()    
}


################################################################################

