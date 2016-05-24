
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


.distStandardErrors <- 
function(fit, obj, x)
{
    # Add Standard Errors and t-Values:
    hessian = tsHessian(x = fit$par, fun = obj, y = x, trace = FALSE)
    colnames(hessian) = rownames(hessian) = names(fit$par)
    fit$cvar = solve(hessian)
    fit$se.coef = sqrt(diag(fit$cvar))
    if (fit$scale) 
        fit$se.coef = fit$se.coef / fit$scaleParams
    fit$tval = fit$par/fit$se.coef
    fit$matcoef = cbind(fit$par, fit$se.coef,
        fit$tval, 2*(1-pnorm(abs(fit$tval))))
    dimnames(fit$matcoef) = list(names(fit$tval), 
        c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
        
    # Return Value:
    fit
}
        
      
################################################################################
  
        