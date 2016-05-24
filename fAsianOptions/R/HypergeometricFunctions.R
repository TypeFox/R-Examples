
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:       KUMMER DESCRIPTION:
#  kummerM         Computes Confluent Hypergeometric Function of the 1st Kind
#  kummerU         Computes Confluent Hypergeometric Function of the 2nd Kind
# FUNCTION:       WHITTAKER DESCRIPTION:
#  whittakerM      Computes Whittaker's M Function
#  whittakerW      Computes Whittaker's M Function
# FUNCTION:       HERMITE POLYNOMIAL:
#  hermiteH        Computes the Hermite Polynomial
################################################################################


################################################################################
# KUMMER:


kummerM =
function(x, a, b, lnchf = 0, ip = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the Confluent Hypergeometric Function of the First
    #   Kind for complex argument "x" and complex indexes "a" and "b"

    # Arguments:
    #   x - complex function argument
    #   a, b - complex indexes
    #   lnchf -
    #   ip -

    # FUNCTION:

    # You can also input real arguments:
    if (!is.complex(x)) x = complex(real = x, imaginary = 0*x)
    if (!is.complex(a)) a = complex(real = a, imaginary = 0)
    if (!is.complex(b)) b = complex(real = b, imaginary = 0)

    # Calculate KummerM:
    chm = rep(complex(real = 0, imaginary = 0), length = length(x))
    value = .Fortran("chfm",
        as.double(Re(x)),
        as.double(Im(x)),
        as.double(Re(a)),
        as.double(Im(a)),
        as.double(Re(b)),
        as.double(Im(b)),
        as.double(Re(chm)),
        as.double(Im(chm)),
        as.integer(length(x)),
        as.integer(lnchf),
        as.integer(ip),
        PACKAGE = "fAsianOptions")
    result = complex(real = value[[7]], imaginary = value[[8]])

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


kummerU =
function(x, a, b, ip = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the Confluent Hypergeometric Function of the Second
    #   Kind for complex argument "x" and complex indexes "a" and "b"

    # Arguments:

    # FUNCTION:

    # Todo ...
    lnchf = 0

    # Test for complex arguments:
    if (!is.complex(x)) x = complex(real = x, imaginary = 0*x)
    if (!is.complex(a)) a = complex(real = a, imaginary = 0)
    if (!is.complex(b)) b = complex(real = b, imaginary = 0)

    # Calculate KummerU:
    # From KummerM:
    # Uses the formula ...
    #   pi/sin(pi*b) [ M(a,b,z) / (Gamma(1+a-b)*Gamma(b)) -
    #        x^(1-b) * M(1+a-b,2-b,z) / (Gamma(a)*Gamma(2-b)) ]
    ans = ( pi/sin(pi*b) ) * (
        kummerM(x, a = a, b = b, lnchf = lnchf, ip=ip) /
            ( cgamma(1+a-b)*cgamma(b) ) - (x^(1-b)) *
        kummerM(x, a = (1+a-b), b=2-b, lnchf = lnchf, ip = ip) /
            ( cgamma(a)*cgamma(2-b) ) )

    # Return Value:
    ans
}


################################################################################
# WHITTAKER:


whittakerM =
function(x, kappa, mu, ip = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Whittaker's M Function

    # Arguments:

    # FUNCTION:

    # Test for complex arguments:
    if (!is.complex(x)) x = complex(real = x, imaginary = 0*x)
    if (!is.complex(kappa)) kappa = complex(real = kappa, imaginary = 0)
    if (!is.complex(mu)) mu = complex(real = mu, imaginary = 0)

    # Calculate:
    ans = exp(-x/2) * x^(1/2+mu) * kummerM(x, 1/2+mu-kappa, 1+2*mu, ip = ip)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


whittakerW =
function(x, kappa, mu, ip = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Whittaker's M Function

    # Arguments:

    # FUNCTION:

    # Test for complex arguments:
    if (!is.complex(x)) x = complex(real = x, imaginary = 0*x)
    if (!is.complex(kappa)) kappa = complex(real = kappa, imaginary = 0)
    if (!is.complex(mu)) mu = complex(real = mu, imaginary = 0)

    # Calculate:
    ans = exp(-x/2) * x^(1/2+mu) * kummerU(x, 1/2+mu-kappa, 1+2*mu, ip = ip)

    # Return Value:
    ans
}


################################################################################
# HERMITE POLYNOMIAL:


hermiteH =
function(x, n, ip = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Hermite Polynomial

    # Arguments:
    #   n - the index of the Hermite polynomial.

    # FUNCTION:

    # Check
    stopifnot(n - round(n, 0) == 0)

    # Result:
    S = sign(x) + (1-sign(abs(x)))
    ans = (S*2)^n * Re ( kummerU(x^2, -n/2, 1/2, ip = ip) )

    # Return Value:
    ans
}


################################################################################

