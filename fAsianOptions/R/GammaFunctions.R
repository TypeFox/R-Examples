
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
# this R-port:
#   by Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
# for the code accessed (or partly included) from other R-ports:
#   R: see R's copyright and license file
# for Haug's Option Pricing Formulas:
#   Formulas are implemented along the book and the Excel spreadsheets of
#     E.G. Haug, "The Complete Guide to Option Pricing"; documentation
#     is partly taken from www.derivicom.com which implements
#     a C Library based on Haug. For non-academic and commercial use
#     we recommend the professional software from "www.derivicom.com".


################################################################################
# FUNCTION:       DESCRIPTION:
#  erf             Error function
#  [gamma]         Gamma function
#  [lgamma]        LogGamma function, returns log(gamma)
#  [digamma]       First derivative of of LogGamma, dlog(gamma(x))/dx
#  [trigamma]      Second derivative of of LogGamma, dlog(gamma(x))/dx
#  {tetragamma}    Third derivative of of LogGamma, dlog(gamma(x))/dx
#  {pentagamma}    Fourth derivative of LogGamma, dlog(gamma(x))/dx
#  [beta]*         Beta function
#  [lbeta]*        LogBeta function, returns log(Beta)
#  Psi             Psi(x) (Digamma) function
#  igamma          P(a,x) Incomplete Gamma Function
#  cgamma          Gamma function for complex arguments
#  Pochhammer      Pochhammer symbol
# NOTES:
#  Functions in [] paranthesis are part of the R's and SPlus' base distribution
#  Functions in {} paranthesis are only availalble in R
#  Function marked by []* are compute through the gamma function in SPlus
################################################################################


erf =
function(x)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Error function for real argument "x"

    # Arguments:
    #   x - a real numeric value or vector.

    # FUNCTION:

    # Result
    # DW 2005-05-04
    ans = 2 * pnorm(sqrt(2) * x) - 1

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


cgamma =
function(x, log = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Gamma Function for complex argument "x"

    # Arguments:
    #   z - a complex or real vector
    #   log - if TRUE the logarithm of the gamma is calculated
    #     otherwise if FALSE, the gamma function itself
    #     will be calculated.

    # Source:
    #   For the Fortran Routine:
    #   http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html

    # FUNCTION:

    # Test for complex arguments:
    if (!is.complex(x)) x = complex(real = x, imaginary = 0*x)

    # Calculate Gamma:
    KF = 1
    if (log) {
        KF = KF - 1
    }
    result = rep(NA, times = length(x))
    for ( i in 1:length(x) ) {
        value = .Fortran("cgama",
            as.double(Re(x[i])),
            as.double(Im(x[i])),
            as.integer(KF),
            as.double(0),
            as.double(0),
            PACKAGE = "fAsianOptions")
        result[i] = complex(real = value[[4]], imaginary = value[[5]])
    }

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


Psi =
function(x)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Psi or Digamma function for complex or real argument

    # Arguments:
    #   z - a complex numeric value or vector.

    # Details:
    #   [AS} formula 6.3.1
    #   $ \Psi(x) = d ln \Gamma(z) / dz = \Gamma prime (z) / \Gamma(z) $

    # Arguments:
    #   x - complex or real vector

    # Source:
    #   For the Fortran Routine:
    #   http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html

    # FUNCTION:

    # Psi:
    result = rep(NA, times = length(x))
    if (!is.complex(x)) {
        # Use R's digamma() function:
        result = digamma(x)
    } else {
        for ( i in 1:length(Re(x)) ) {
            value = .Fortran("cpsi",
                as.double(Re(x[i])),
                as.double(Im(x[i])),
                as.double(0),
                as.double(0),
                PACKAGE = "fAsianOptions")
            result[i] = complex(real = value[[3]], imaginary = value[[4]])
        }
    }

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


igamma =
function(x, a)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the Incomplete Gamma Function P(a, x) with
    #   Re(a) > 0 for complex or real argument "x" and for
    #   complex or real index "z"

    # Arguments:
    #   z - a complex or real vector
    #   a - a complex or real numeric value

    # Details:
    #   [AS] formula 6.5.1
    #   $ frac{1}{Gamma(a)}  * \int_0^x e^{-t} t^{a-1} dt $

    # FUNCTION:

    # igamma:
    if (!is.complex(x) && !is.complex(a)) {
        # Use R's pgamma() function:
        # if (a < 0) Not suppported ...
        result = pgamma(x, a)
    } else {
        # Why not derive the result from KummersM ?
        log = FALSE
        if (log) {
            # Not yet supported:
            result = kummerM(a, a + 1, -x, lnchf = 1) + a*log(x) - log(a)
        } else {
            result = kummerM(a, a + 1, -x, lnchf = 0) * x^a / a
        }
    }

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


Pochhammer =
function(x, n)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Pochhammer's Symbol

    # Arguments:
    #   x - a complex numeric value or vector.
    #   n - an integer n >=0. An notation used in the theory of special
    #       functions for the rising factorial, also known as the rising
    #       factorial power (Graham et al. 1994).

    # Details:
    #   as defined in [AS] by formula 6.1.22

    # FUNCTION:

    # Note:
    # $ (z)_0 = 1 $
    # $ (z)_n = z(z+1)(z+2) \dots (z+n-1) = frac{\Gamma(z+n)}{Gamma(z)} $

    # In case of wrong argument Type:
    Pochhammer = NA

    # For Complex Arguments:
    if (is.complex(x)) {
        Pochhammer = cgamma(x + n)/cgamma(x)
    }

    # For Real Arguments:
    # DW: 2006-05-10 is.real(z) replaced by is.real(x)
    # YC: is.real is deprecated -> replaced by is.double
    if (is.double(x)) {
        Pochhammer = gamma(x + n)/gamma(x)
    }

    # Return Value:
    Pochhammer
}


################################################################################
# SPlus Addon for beta() and lbeta()
# quick and dirty implementation ...


.S = FALSE


# ------------------------------------------------------------------------------


if (.S) {
    beta =
    function(a, b)
    {   # A function implemented by Diethelm Wuertz

        # Description:
        #   Computes the beta function

        # Result:
        ans = gamma(a) * gamma(b) / gamma(a+b)

        # Return Value:
        ans
    }


    lbeta =
    function(a, b)
    {   # A function implemented by Diethelm Wuertz

        # Description:
        #   Computes the beta function

        # Result:
        ans = lgamma(a) + lgamma(b) - lgamma(a+b)

        # Return Value:
        ans
    }
}


################################################################################
