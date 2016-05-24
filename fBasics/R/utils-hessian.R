
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


################################################################################
# FUNCTION:               DESCRIPTION:
#  tsHessian               Computes Two Sided Hessian matrix
################################################################################


tsHessian <-
function(x, fun, ...)
{
    # Description:
    #   Computes two sided (TS) approximated Hessian
    
    # Source:
    #   A function borrowed from Kevin Sheppard's Matlab garch toolbox
    #   as implemented by Alexios Ghalanos in his rgarch package
    
    # Notes:
    #   optionally requires package Matrix (added as suggestion)
    
    # FUNCTION:
    
    # Settings:
    n = length(x)
    fx <- fun(x, ...)
    eps = .Machine$double.eps
    
    # Compute the stepsize:
    h = eps^(1/3) * 
        apply( as.data.frame(x), 1, FUN = function(z) max(abs(z), 1.0e-2))
    xh = x + h
    h = xh - x
    ee = diag(h) # Matrix(diag(h), sparse = TRUE)
    
    # Compute forward and backward steps:
    gp = vector(mode = "numeric", length = n)
    for(i in 1:n) gp[i] <- fun(x + ee[, i], ...)
    gm = vector(mode = "numeric", length = n)
    for(i in 1:n) gm[i] <- fun(x - ee[, i], ...)
    H = h %*% t(h)
    Hm = H
    Hp = H

    # Compute double forward and backward steps:
    for(i in 1:n){
        for(j in  i:n){
            Hp[i, j] <- fun(x + ee[, i] + ee[, j], ...)
            Hp[j, i] <- Hp[i, j]                
            Hm[i, j] <- fun(x - ee[, i] - ee[, j], ...)
            Hm[j, i] <- Hm[i, j]
        }
    }

    # Compute the Hessian:
    for(i in 1:n){
        for(j in  i:n){
            H[i, j] = ((Hp[i, j] - gp[i] - gp[j] + fx + fx - gm[i] - gm[j] +
                Hm[i, j]) / H[i, j]) / 2
            H[j, i] = H[i, j]
        }
    }
    
    # Return Value:
    return(H)
}


################################################################################

