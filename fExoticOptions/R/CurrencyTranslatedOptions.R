
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
# FUNCTION:                     DESCRIPTION:
#  Currency Translated Options:
#   FEInDomesticFXOption          FX In Domestic Currency
#   QuantoOption                  Quanto Option
#   EquityLinkedFXOption          EquityLinked FX Option
#   TakeoverFXOption              Takeover FX Option
################################################################################


FEInDomesticFXOption =
    function(TypeFlag = c("c", "p"), S, E, X, Time, r, q, sigmaS, sigmaE, rho,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Foreign equity option struck in domestic currency

    # References:
    #   Haug, Chapter 2.13.1

    # FUNCTION:

    # Compute Settings:
    TypeFlag = TypeFlag[1]
    sigma = sqrt(sigmaE^2 + sigmaS^2 + 2*rho*sigmaE*sigmaS)
    d1 = (log(E*S/X) + (r-q+sigma^2/2) * Time) / (sigma*sqrt(Time))
    d2 = d1 - sigma * sqrt(Time)

    # Calculate Call and Put:
    if (TypeFlag == "c") {
        ForeignEquityInDomesticFX = (E * S * exp(-q*Time)*CND(d1) -
                                     X * exp(-r*Time)*CND(d2)) }
    if (TypeFlag == "p") {
        ForeignEquityInDomesticFX = (X * exp(-r*Time)*CND(-d2) -
                                     E * S * exp(-q*Time)*CND(-d1)) }

    # Parameters:
    # TypeFlag = c("c", "p"), S, E, X, Time, r, q, sigmaS, sigmaE, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$E = E
    param$X = X
    param$Time = Time
    param$r = r
    param$q = q
    param$sigmaS = sigmaS
    param$sigmaE = sigmaE
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "FE In Domestic FX Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = ForeignEquityInDomesticFX,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


QuantoOption =
    function(TypeFlag = c("c", "p"), S, Ep, X, Time, r, rf, q, sigmaS,
             sigmaE, rho, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fixed exchange rate foreign equity options

    # References:
    #   Haug, Chapter 2.13.2

    # FUNCTION:

    # Compute Settings:
    TypeFlag = TypeFlag[1]
    d1 = ((log(S/X) + (rf-q-rho*sigmaS*sigmaE + sigmaS^2/2) * Time) /
          (sigmaS*sqrt(Time)))
    d2 = d1 - sigmaS*sqrt (Time)

    # Calculate Call and Put:
    if (TypeFlag == "c") {
        Quanto = (Ep*(S*exp((rf-r-q-rho*sigmaS*sigmaE)*Time) *
                      CND(d1) - X*exp(-r*Time)*CND(d2))) }
    if (TypeFlag == "p") {
        Quanto = (Ep*(X*exp(-r*Time)*CND(-d2) -
                      S*exp((rf-r-q-rho*sigmaS*sigmaE)* Time)*CND(-d1))) }

    # Parameters:
    # TypeFlag = c("c", "p"), S, Ep, X, Time, r, rf, q, sigmaS, sigmaE, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$Ep = Ep
    param$X = X
    param$Time = Time
    param$r = r
    param$rf = rf
    param$q = q
    param$sigmaS = sigmaS
    param$sigmaE = sigmaE
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Quanto Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = Quanto,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


EquityLinkedFXOption =
    function(TypeFlag = c("c", "p"), E, S, X, Time, r, rf, q, sigmaS,
             sigmaE, rho, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Equity Linked FX Option -

    # References:
    #   Haug, Chapter 2.13.3

    # FUNCTION:

    # Compute Settings:
    TypeFlag = TypeFlag[1]
    vS = sigmaS
    vE = sigmaE
    d1 = ((log(E / X) + (r - rf + rho * vS * vE + vE ^ 2 / 2) * Time) /
          (vE * sqrt(Time)))
    d2 = d1 - vE * sqrt(Time)

    # Calculate Call and Put:
    if (TypeFlag == "c") {
        EquityLinkedFXO = (E * S * exp(-q * Time) * CND(d1) -
                           X * S * exp((rf - r - q - rho * vS * vE) * Time) *
                           CND(d2)) }
    if (TypeFlag == "p") {
        EquityLinkedFXO = (X * S * exp((rf - r - q - rho * vS * vE) * Time) *
                           CND(-d2) - E * S * exp(-q * Time) * CND(-d1)) }

    # Parameters:
    # TypeFlag = c("c", "p"), E, S, X, Time, r, rf, q, sigmaS, sigmaE, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$E = E
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$rf = rf
    param$q = q
    param$sigmaS = sigmaS
    param$sigmaE = sigmaE
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Equity Linked FX Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = EquityLinkedFXO,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


TakeoverFXOption =
    function(V, B, E, X, Time, r, rf, sigmaV, sigmaE, rho,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Takeover FX  Option -

    # References:
    #   Haug, Chapter 2.13.4

    # FUNCTION:

    # Compute Settings:
    v = V
    b = B
    vV = sigmaV
    vE = sigmaE
    a1 = ((log(v / b) + (rf - rho * vE * vV - vV ^ 2 / 2) * Time) /
          (vV * sqrt(Time)))
    a2 = ((log(E / X) + (r - rf - vE ^ 2 / 2) * Time) /
          (vE * sqrt(Time)))

    # Calculate:
    TakeoverFX = (b * (E * exp(-rf * Time) * CBND(a2 + vE *
                  sqrt(Time), -a1 - rho * vE * sqrt(Time), -rho) - X *
                  exp(-r * Time) * CBND(-a1, a2, -rho)))

    # Parameters:
    # V, B, E, X, Time, r, rf, sigmaV, sigmaE, rho
    param = list()
    param$V = V
    param$B = B
    param$E = E
    param$X = X
    param$Time = Time
    param$r = r
    param$rf = rf
    param$q = q
    param$sigmaV = sigmaV
    param$sigmaE = sigmaE
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Takeover FX Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = TakeoverFX,
        title = title,
        description = description
        )
}


################################################################################

