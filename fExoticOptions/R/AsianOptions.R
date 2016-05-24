
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
# FUNCTION:                         DESCRIPTION:
# Asian Options:
#  GeometricAverageRateOption        Geometric Average Rate Option
#  TurnbullWakemanAsianApproxOption  Turnbull-Wakeman Approximated Asian Option
#  LevyAsianApproxOption             Levy Approximated Asian Option
################################################################################


GeometricAverageRateOption =
    function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Valuates geometric average rate options

    # References:
    #   Kemma and Vorst (1990)
    #   Haug, Chapter 2.12.1

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    b.A = 0.5 * (b - sigma^2 / 6)
    sigma.A = sigma / sqrt (3)
    GeometricAverageRate =
        GBSOption (TypeFlag = TypeFlag, S = S, X = X, Time = Time,
                   r = r, b = b.A, sigma = sigma.A)@price

    # Parameters:
    # TypeFlag = c("c", "p"), S, X, Time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Geometric Average Rate Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = GeometricAverageRate,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


TurnbullWakemanAsianApproxOption =
    function(TypeFlag = c("c", "p"), S, SA, X, Time, time, tau, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Valuates arithmetic average rate options by the
    #   Turnbull-Wakeman's Approximation

    # References:
    #   Haug, Chapter 2.12.2

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    m1 = (exp(b * Time) - exp(b * tau)) / (b * (Time - tau))
    m2 = (2 * exp((2 * b + sigma^2) * Time) / ((b + sigma^2) * (2*b +
          sigma^2) * (Time - tau)^2) + 2 * exp((2 * b + sigma^2) *
          tau) / (b * (Time - tau)^2) * (1/(2 * b + sigma^2) - exp(b *
          (Time - tau)) / (b + sigma^2)))
    b.A = log(m1) / Time
    sigma.A = sqrt(log(m2) / Time - 2*b.A)
    t1 = Time - time
    if (t1 > 0) {
        X = Time/time * X - t1/time * SA
        TurnbullWakemanAsianApprox =
            (GBSOption(TypeFlag, S, X, time, r, b.A, sigma.A)@price *
             time/Time) }
    else {
        TurnbullWakemanAsianApprox =
            GBSOption(TypeFlag, S, X, time, r, b.A, sigma.A)@price }

    # Parameters:
    # TypeFlag = c("c", "p"), S, SA, X, Time, time, tau, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$SA = SA
    param$X = X
    param$Time = Time
    param$time = time
    param$tau = tau
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Turnbull Wakeman Asian Approximated Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = TurnbullWakemanAsianApprox,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


LevyAsianApproxOption =
    function(TypeFlag = c("c", "p"), S, SA, X, Time, time, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Valuates arithmetic average rate options by the
    #   Levy Approximation

    # References:
    #   Haug, Chapter 2.12.2

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    SE = S / (Time*b) * (exp((b-r)*time) - exp(-r*time))
    m = 2 * S ^ 2 / (b + sigma ^ 2) * ((exp((2 *
    b + sigma^2) * time) - 1) / (2 * b + sigma^2) -
    (exp(b * time) - 1) / b)
    d = m / (Time^2)
    Sv = log (d) - 2 * (r * time + log(SE))
    XStar = X - (Time - time) / Time * SA
    d1 = 1 / sqrt (Sv) * (log(d) / 2 - log(XStar))
    d2 = d1 - sqrt (Sv)
    if (TypeFlag == "c") {
        LevyAsianApprox = (SE * CND (d1) - XStar * exp(-r*time) *
                           CND(d2))}
    if (TypeFlag == "p") {
        LevyAsianApprox = ((SE * CND(d1) - XStar * exp(-r*time) *
                            CND(d2)) - SE + XStar * exp (-r*time)) }

    # Parameters:
    # TypeFlag = c("c", "p"), S, SA, X, Time, time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$SA = SA
    param$X = X
    param$Time = Time
    param$time = time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Levy Asian Approximated Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = LevyAsianApprox,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


#CurranAsianApproxOption =
#function()
#{  # A function implemented by Diethelm Wuertz

# Description:
#   Arithmetic average rate option
#   Curran's Approximation

# References:
#   Haug, Chapter 2.12.2

# FUNCTION:

# Compute Price:
#   CurranAsianApprox = NA

# Return Value:
#   CurranAsianApprox
#}


################################################################################

