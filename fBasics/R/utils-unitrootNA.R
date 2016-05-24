
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
# FUNCTION:                 DESCRIPTION:
#  .unirootNA                Computes zero of a function without error exit
################################################################################


.unirootNA <- 
function(f, interval, lower = min(interval), upper = max(interval),
    tol = .Machine$double.eps^0.25, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Searches the interval from lower to upper for a
    #   root (i.e., zero) of the function f with respect
    #   to its first argument.

    # Arguments:
    #   see 'uniroot'

    # Value:
    #   Returns the x value of f where the root is located. If
    #   now root exists NA will be returned. In the last case
    #   the function doesn't terminate with an error like in the
    #   case of the standard function uniroot.

    # Details:
    #   R:
    #   uniroot(f, interval, lower = min(interval), upper = max(interval),
    #       tol = .Machine$double.eps^0.25,
    #       maxiter = 1000, ...)
    #   uniroot(f, interval, lower = min(interval), upper = max(interval),
    #       tol = .Machine$double.eps^.25,
    #       keep.xy = F, f.lower = NA,  f.upper = NA, ...)

    # Example:
    #   .unirootNA(sin, c(1, 2)); .unirootNA(sin, c(-1, 1))

    # FUNCTION:

    # There is no Root:
    if (is.null(args(f))) {
        if (f(lower) * f(upper) >=0) return(NA)
    } else {
        if (f(lower, ...) * f(upper, ...) >= 0) return(NA)
    }

    # There is a Root:
    ans = uniroot(f = f, interval = interval, lower = lower,
        upper = upper, tol = tol, ...)

    # Return Value:
    ans$root
}


################################################################################

