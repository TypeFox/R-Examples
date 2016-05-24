
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
# FUNCTION:               PARAMETER ESTIMATION:
#  .garchRoptimhess        Uses R internal optimhess function
#  .garchRCDAHessian       Computes R coded CDA Hessian matrix
#  .garchRTSHessian        Computes R coded Two Sided Hessian matrix
#  .hessian2sided          Function called from .garchRTSHessian
################################################################################


.garchRoptimhess <-
    function(par, .params, .series, eps = 1.0e-4)
{
    # A function implemeted by Diethelm Wuertz

    # Description:
    #   Compute Hessian via R's function optimHess()

    # Arguments:
    #   par -
    #   .params -
    #   .series -
    #   eps -

    # FUNCTION:

    # Take Time:
    .StartHessian <- Sys.time()

    # Compute Hessian:
    H <- optimHess(par, .garchLLH)
    H <- 0.5 * (H + t(H))
    nm <- names(par)
    dimnames(H) <- list(nm, nm)

    # Elapsed Time:
    time = Sys.time() - .StartHessian
    attr(H, "time") = time

    # Return Value:
    H
}


# ------------------------------------------------------------------------------


.garchRCDAHessian <-
    function(par, .params, .series, eps = 1.0e-4)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute CDA (central difference approximated) Hessian

    # Arguments:
    #   par -
    #   .params -
    #   .series -
    #   eps -

    # Reference:
    #   http://status.sph.umich.edu/computing/manuals/sas8/stat/chap46/sect26.htm

    # FUNCTION:

    # Take start time
    .StartHessian <- Sys.time()

    # Algorithm:
    algorithm = .params$control$algorithm[1]
    .trace = FALSE

    # Compute Hessian:
    eps = eps * par
    n = length(par)
    H = matrix(0, ncol = n, nrow = n)
    for (i in 1:n) {
        for (j in 1:n) {
            x1 = x2 = x3 = x4 = par
            x1[i] = x1[i] + eps[i]
            x1[j] = x1[j] + eps[j]
            x2[i] = x2[i] + eps[i]
            x2[j] = x2[j] - eps[j]
            x3[i] = x3[i] - eps[i]
            x3[j] = x3[j] + eps[j]
            x4[i] = x4[i] - eps[i]
            x4[j] = x4[j] - eps[j]
            H[i, j] = (
             .garchLLH(x1, .trace) -
             .garchLLH(x2, .trace) -
             .garchLLH(x3, .trace) +
             .garchLLH(x4, .trace) ) / (4*eps[i]*eps[j])
        }
    }
    colnames(H) = rownames(H) = names(par)

    # Attribute execution time:
    time = Sys.time() - .StartHessian
    attr(H, "time") = time

    # Return Value:
    H
}


# ------------------------------------------------------------------------------


.garchTSHessian <-
    function(par, .params, .series, eps = NA)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute two sided (TS) approximated Hessian

    # Arguments:
    #   par -
    #   .params -
    #   .series -
    #   eps - not used

    # FUNCTION:

    # Take start time
    .StartHessian <- Sys.time()

    # Algorithm:
    algorithm = .params$control$algorithm[1]

    # Compute Hessian:
    H <- .hessian2sided(f = .garchLLH, x = par, trace = FALSE, fGarchEnv = FALSE)
    colnames(H) = rownames(H) = names(par)

    # Attribute execution time:
    time = Sys.time() - .StartHessian
    attr(H, "time") = time

    # Return Value:
    H
}


# ------------------------------------------------------------------------------


.hessian2sided <-
    function(f, x, ...)
{
    # A function adapted from Kevin Sheppard's Matlab garch toolbox
    # ... implemented by Alexios Ghalanos in his rgarch package
    # ... R port for Rmetrics' fGarch by Diethelm Wuertz

    # Description:
    #   Computes two sided (TS) approximated Hessian

    # Arguments:
    #   f -
    #   x -

    # Notes:
    #   requires package Matrix (added as suggestion)

    # FUNCTION:

    # Settings:
    n = length(x)
    fx <- f(x, ...)
    eps = .Machine$double.eps

    # Compute the stepsize (h)
    h = eps^(1/3) *
        apply( as.data.frame(x), 1, FUN = function(z) max(abs(z), 1.0e-2))
    xh = x + h
    h = xh - x
    ee = Matrix(diag(h), sparse = TRUE)

    # Compute forward and backward steps:
    gp = vector(mode = "numeric", length = n)
    for(i in 1:n) gp[i] <- f(x + ee[, i], ...)
    gm = vector(mode = "numeric", length = n)
    for(i in 1:n) gm[i] <- f(x - ee[, i], ...)
    H = h %*% t(h)
    Hm = H
    Hp = H

    # Compute "double" forward and backward steps:
    for(i in 1:n){
        for(j in  i:n){
            Hp[i, j] <- f(x + ee[, i] + ee[, j], ...)
            Hp[j, i] <- Hp[i, j]
            Hm[i, j] <- f(x - ee[, i] - ee[, j], ...)
            Hm[j, i] <- Hm[i, j]
        }
    }

    # Compute the hessian:
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

