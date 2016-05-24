
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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


##############################################################################
# FUNCTION:               DESCRIPTION:    
#  garchFitControl         Sets default values for Garch Optimizer
##############################################################################
  

garchFitControl <- 
function(
    llh = c("filter", "internal", "testing"),
    nlminb.eval.max = 2000, 
    nlminb.iter.max = 1500,
    nlminb.abs.tol = 1.0e-20, 
    nlminb.rel.tol = 1.0e-14, 
    nlminb.x.tol = 1.0e-14, 
    nlminb.step.min = 2.2e-14,
    nlminb.scale = 1,  
    nlminb.fscale = FALSE,
    nlminb.xscale = FALSE,     
    sqp.mit = 200,       
    sqp.mfv = 500,       
    sqp.met = 2,                                 
    sqp.mec = 2,                                
    sqp.mer = 1,                                  
    sqp.mes = 4,                                 
    sqp.xmax = 1.0e3,    
    sqp.tolx = 1.0e-16,    
    sqp.tolc = 1.0e-6,  
    sqp.tolg = 1.0e-6,
    sqp.told = 1.0e-6,   
    sqp.tols = 1.0e-4,  
    sqp.rpf = 1.0e-4,
    lbfgsb.REPORT = 10,
    lbfgsb.lmm = 20, 
    lbfgsb.pgtol = 1e-14, 
    lbfgsb.factr = 1, 
    lbfgsb.fnscale = FALSE,
    lbfgsb.parscale = FALSE,   
    nm.ndeps = 1e-14,
    nm.maxit = 10000, 
    nm.abstol = 1e-14,
    nm.reltol = 1e-14, 
    nm.alpha = 1.0, 
    nm.beta = 0.5, 
    nm.gamma = 2.0,
    nm.fnscale = FALSE, 
    nm.parscale = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets default values for Garch Optimizer
    
 
    # FUNCTION:
    
    # Generate Control List with Default Settings:
    control <- list(
    
       llh = llh,

       nlminb.eval.max = nlminb.eval.max, 
       nlminb.iter.max = nlminb.iter.max,
       nlminb.abs.tol = nlminb.abs.tol, 
       nlminb.rel.tol = nlminb.rel.tol, 
       nlminb.x.tol = nlminb.x.tol, 
       nlminb.step.min = nlminb.step.min,
       nlminb.scale = nlminb.scale,
       nlminb.fscale = nlminb.fscale,
       nlminb.xscale = nlminb.xscale,
         
       sqp.mit = sqp.mit,       
       sqp.mfv = sqp.mfv,       
       sqp.met = sqp.met,                                 
       sqp.mec = sqp.mec,                                
       sqp.mer = sqp.mer,                                  
       sqp.mes = sqp.mes,                                 
       sqp.xmax = sqp.xmax,    
       sqp.tolx = sqp.tolx, 
       sqp.tolc = sqp.tolc,    
       sqp.tolg = sqp.tolg,   
       sqp.told = sqp.told,  
       sqp.tols = sqp.tols,  
       sqp.rpf = sqp.rpf, 

       lbfgsb.REPORT = lbfgsb.REPORT,
       lbfgsb.lmm = lbfgsb.lmm, 
       lbfgsb.pgtol = lbfgsb.pgtol, 
       lbfgsb.factr = lbfgsb.factr,  
       lbfgsb.fnscale = lbfgsb.fnscale,
       lbfgsb.parscale = lbfgsb.parscale,  
   
       nm.ndeps = nm.ndeps,
       nm.maxit = nm.maxit, 
       nm.abstol = nm.abstol,
       nm.reltol = nm.reltol, 
       nm.alpha = nm.alpha, 
       nm.beta = nm.beta, 
       nm.gamma = nm.gamma,
       nm.fnscale = nm.fnscale, 
       nm.parscale = nm.parscale
       
       ) 
        
    # Return Value:
    control
}


################################################################################

