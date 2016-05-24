
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


################################################################################
# FUNCTION:               PARAMETER ESTIMATION:
#  .garchRCDAGradient      Computes R coded CDA matrix of contributions
#				                   to the Gradient
################################################################################


.garchRCDAGradient <-
    function(par, .params, .series, eps = 1.0e-4)
{
    # A function implemented by Michal Miklovic & Yohan Chalabi

    # Description:
    #   Compute R coded CDA (central difference approximated) Gradient

    # Reference:
    #   http://status.sph.umich.edu/computing/manuals/sas8/stat/chap46/sect26.htm

    # FUNCTION:

    # Starttime
    .StartGradient <- Sys.time()

    # Algorithm
    algorithm = .params$control$algorithm[1]
    .trace = FALSE


    # LLH for the computation of matrix of contributions to the Gradient
    skew <- .params$skew
    shape <- .params$shape
    delta <- .params$delta
    deltainv = 1/delta
    llh.start = .series$llh.start
    N <- length(.series$x)
    .garchDist <- .getfGarchEnv(".garchDist")

    # Compute matrix of contributions to the Gradient:
    eps = eps * par
    n = N - llh.start + 1
    K = length(par)
    G = matrix(0, nrow = n, ncol = K)

    for (i in 1:K) {
        x1 = x2 = par
        x1[i] = x1[i] + eps[i]
        x2[i] = x2[i] - eps[i]
        #
        .garchLLH(x1, .trace, TRUE)
        h1 <- .getfGarchEnv(".series")$h
        z1 <- .getfGarchEnv(".series")$z
        hh1 = (abs(h1[(llh.start):N]))^deltainv
        zz1 = z1[(llh.start):N]
        llh.grad1 <-
            log(.garchDist(z = zz1, hh = hh1, skew = skew, shape = shape))
        #
        .garchLLH(x2, .trace, TRUE)
        h2 <- .getfGarchEnv(".series")$h
        z2 <- .getfGarchEnv(".series")$z
        hh2 = (abs(h2[(llh.start):N]))^deltainv
        zz2 = z2[(llh.start):N]
        llh.grad2 <-
            log(.garchDist(z = zz2, hh = hh2, skew = skew, shape = shape))
        #
        G[,i] = (llh.grad1 - llh.grad2) / (2*eps[i])
    }

    rownames(G) = c(1:n)
    colnames(G) = names(par)

    # make sure that h and z are ok
    .setfGarchEnv(.series = .series)

    time = Sys.time() - .StartGradient

    # Attribute Exdecution time
    attr(G, "time") = time

    # Return Value:
    G
}

################################################################################
