
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
# FUNCTION:                       DESCRIPTION:
# Multiple Asset Options:
#   TwoAssetCorrelationOption       Two Asset Correlation Option
#    [ExchangeOneForAnotherOption]    [Exchange One For Another Option]
#   EuropeanExchangeOption          European Exchange Optionn
#   AmericanExchangeOption          American Exchange Option
#   ExchangeOnExchangeOption        Exchange Exchange Option
#   TwoRiskyAssetsOption            Option On The MinMax
#   SpreadApproxOption              Spread Approximated Option
################################################################################


TwoAssetCorrelationOption =
    function(TypeFlag = c("c", "p"), S1, S2, X1, X2, Time, r, b1, b2,
             sigma1, sigma2, rho, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Two asset correlation options

    # References:
    #   Haug, Chapter 2.8.1

    # FUNCTION:

    # Compute Settings:
    TypeFlag = TypeFlag[1]
    y1 = (log(S1/X1) + (b1 - sigma1^2 / 2) * Time) / (sigma1*sqrt(Time))
    y2 = (log(S2/X2) + (b2 - sigma2^2 / 2) * Time) / (sigma2*sqrt(Time))

    # Calculate Call and Put:
    if (TypeFlag == "c")
        TwoAssetCorrelation = (S2 * exp ((b2 - r) * Time) * CBND(y2 +
                               sigma2 * sqrt(Time), y1 + rho * sigma2
                               * sqrt(Time), rho) - X2 * exp (-r *
                               Time) * CBND(y2, y1, rho))
    if (TypeFlag == "p")
        TwoAssetCorrelation = (X2 * exp (-r * Time) * CBND(-y2, -y1,
                               rho) - S2 * exp ((b2 - r) * Time) *
                               CBND(-y2 - sigma2 * sqrt(Time), -y1 -
                               rho * sigma2 * sqrt(Time), rho))

    # Parameters:
    # TypeFlag = c("c", "p"), S1, S2, X1, X2, Time, r, b1, b2, sigma1,
    #   sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$X1 = X1
    param$X2 = X2
    param$Time = Time
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Two Asset Correlation Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = TwoAssetCorrelation,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


EuropeanExchangeOption =
    function(S1, S2, Q1, Q2, Time, r, b1, b2, sigma1, sigma2, rho,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Exchange-One-Asset-for-Another-Asset options -
    #   European option to exchange one asset for another

    # References:
    #   Haug, Chapter 2.8.2 (European)

    # FUNCTION:

    # Compute Settings:
    sigma = sqrt (sigma1 ^ 2 + sigma2 ^ 2 - 2 * rho * sigma1 * sigma2)
    d1 = ((log(Q1*S1/(Q2 * S2)) + (b1-b2+sigma^2/2)*Time)/(sigma*sqrt(Time)))
    d2 = d1 - sigma * sqrt (Time)

    # calculate Price:
    EuropeanExchange = (Q1 * S1 * exp ((b1 - r) * Time) * CND(d1) -
                        Q2 * S2 * exp((b2 - r) * Time) * CND(d2))

    # Parameters:
    # S1, S2, Q1, Q2, Time, r, b1, b2, sigma1, sigma2, rho
    param = list()
    param$S1 = S1
    param$S2 = S2
    param$Q1 = Q1
    param$Q2 = Q2
    param$Time = Time
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "European Exchange Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = EuropeanExchange,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


AmericanExchangeOption =
    function(S1, S2, Q1, Q2, Time, r, b1, b2, sigma1, sigma2, rho,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Exchange-One-Asset-for-Another-Asset options -
    #   American option to exchange one asset for another

    # References:
    #   Haug, Chapter 2.8.2 (American)

    # FUNCTION:

    # Compute Settings:
    sigma = sqrt(sigma1^2 + sigma2^2 - 2 * rho * sigma1 * sigma2)

    # Calculate Price:
    AmericanExchange = BSAmericanApproxOption("c", Q1*S1, Q2*S2,
    Time, r-b2, b1-b2, sigma)

    # Parameters:
    # S1, S2, Q1, Q2, Time, r, b1, b2, sigma1, sigma2, rho
    param = list()
    param$S1 = S1
    param$s2 = S2
    param$Q1 = Q2
    param$Time = Time
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho
    param$TriggerPrice = AmericanExchange@parameters$TriggerPrice

    # Add title and description:
    if (is.null(title)) title = "American Exchange Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = AmericanExchange@price,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


ExchangeOnExchangeOption =
    function(TypeFlag = c("1", "2", "3", "4"), S1, S2, Q, time1, Time2, r,
             b1, b2, sigma1, sigma2, rho, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Exchange-One-Asset-for-Another-Asset options -

    # References:
    #   Haug, Chapter 2.8.3

    # FUNCTION:

    # Define Functions:
    TypeFlag = TypeFlag[1]
    q = Q

    # Compute:
    v = sqrt(sigma1 ^ 2 + sigma2 ^ 2 - 2 * rho * sigma1 * sigma2)
    I1 = S1 * exp((b1 - r) * (Time2 - time1)) /
        (S2 * exp((b2 - r) * (Time2 - time1)))
    if (TypeFlag == "1" || TypeFlag == "2") {
        id = 1 }
    else {
        id = 2 }
    I = .EOnEOption.CriticalPrice(id, I1, time1, Time2, v, q)

    d1 = (log(S1 / (I * S2)) + (b1 - b2 + v ^ 2 / 2) * time1) / (v * sqrt(time1))
    d2 = d1 - v * sqrt(time1)
    d3 = (log((I * S2) / S1) + (b2 - b1 + v ^ 2 / 2) * time1) / (v * sqrt(time1))
    d4 = d3 - v * sqrt(time1)
    y1 = (log(S1 / S2) + (b1 - b2 + v ^ 2 / 2) * Time2) / (v * sqrt(Time2))
    y2 = y1 - v * sqrt(Time2)
    y3 = (log(S2 / S1) + (b2 - b1 + v ^ 2 / 2) * Time2) / (v * sqrt(Time2))
    y4 = y3 - v * sqrt(Time2)

    # Calculate Price:
    if (TypeFlag == "1")
        ExchangeOnExchange = (-S2 * exp((b2 - r) * Time2) * CBND(d2,
                              y2, sqrt(time1/Time2)) + S1 * exp((b1-r)
                              * Time2) * CBND(d1, y1,
                              sqrt(time1/Time2)) - q * S2 * exp((b2-r)
                              * time1) * CND(d2))
    if (TypeFlag == "2")
        ExchangeOnExchange = (S2 * exp((b2 - r) * Time2) * CBND(d3,
                              y2, -sqrt(time1/Time2)) - S1 *
                              exp((b1-r) * Time2) * CBND(d4, y1,
                              -sqrt(time1/Time2)) + q * S2 * exp((b2 -
                              r) * time1) * CND(d3))
    if (TypeFlag == "3")
        ExchangeOnExchange = (S2 * exp((b2 - r) * Time2) * CBND(d3,
                              y3, sqrt(time1/Time2)) - S1 * exp((b1-r)
                              * Time2) * CBND(d4, y4,
                              sqrt(time1/Time2)) - q * S2 * exp((b2-r)
                              * time1) * CND(d3))
    if (TypeFlag == "4")
        ExchangeOnExchange = (-S2 * exp((b2 - r) * Time2) * CBND(d2,
                              y3, -sqrt(time1/Time2)) + S1 *
                              exp((b1-r) * Time2) * CBND(d1, y4,
                              -sqrt(time1/Time2)) + q * S2 *
                              exp((b2-r) * time1) * CND(d2))

    # Parameters:
    # TypeFlag = c("1", "2", "3", "4"), S1, S2, Q, time1, Time2, r,
    #   b1, b2, sigma1, sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$Q = Q
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Exchange On Exchange Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = ExchangeOnExchange,
        title = title,
        description = description
        )
}


.EOnEOption.CriticalPart3 <-
    function(id, I, time1, Time2, v)
{
    if (id == 1) {
        z1 = (log(I)+v^2/2*(Time2 - time1))/(v*sqrt(Time2-time1))
        z2 = (log(I)-v^2/2*(Time2 - time1))/(v*sqrt(Time2-time1))
        .EOnEOption.CriticalPart3 = I * CND(z1) - CND(z2) }
    if (id == 2) {
        z1 = (-log(I)+v^2/2*(Time2-time1))/(v*sqrt(Time2-time1))
        z2 = (-log(I)-v^2/2*(Time2-time1))/(v*sqrt(Time2-time1))
        .EOnEOption.CriticalPart3 = CND(z1) - I * CND(z2) }
    .EOnEOption.CriticalPart3
}


.EOnEOption.CriticalPart2 <-
    function(id, I, time1, Time2, v)
{
    if (id == 1) {
        z1 = (log(I)+v^2/2*(Time2-time1))/(v*sqrt(Time2-time1))
        .EOnEOption.CriticalPart2 = CND(z1) }
    if (id == 2) {
        z2 = (-log(I)-v^2/2*(Time2-time1))/(v*sqrt(Time2-time1))
        .EOnEOption.CriticalPart2 = -CND(z2) }
    .EOnEOption.CriticalPart2
}


.EOnEOption.CriticalPrice =
    function(id, I1, time1, Time2, v, q)
{
    # Numerical search algorithm to find critical price I
    Ii = I1
    yi = .EOnEOption.CriticalPart3(id, Ii, time1, Time2, v)
    # cat("\n.EOnEOption.CriticalPart3: ", yi)
    di = .EOnEOption.CriticalPart2(id, Ii, time1, Time2, v)
    # cat("\n.EOnEOption.CriticalPart2: ", di)
    epsilon = 0.00001
    while (abs(yi - q) > epsilon) {
        Ii = Ii - (yi - q) / di
        yi = .EOnEOption.CriticalPart3(id, Ii, time1, Time2, v)
        # cat("\n.EOnEOption.CriticalPart3: ", yi)
        di = .EOnEOption.CriticalPart2(id, Ii, time1, Time2, v)
        # cat("\n.EOnEOption.CriticalPart2: ", di)
    }
    .EOnEOption.CriticalPrice = Ii
    .EOnEOption.CriticalPrice
}


# ------------------------------------------------------------------------------


TwoRiskyAssetsOption =
    function(TypeFlag = c("cmin", "cmax", "pmin", "pmax"), S1, S2, X, Time,
             r, b1, b2, sigma1, sigma2, rho, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Option on two risky assets

    # References:
    #   Haug, Chapter 2.8.4

    # FUNCTION:

    # Compute Settings:
    TypeFlag = TypeFlag[1]
    v = sqrt(sigma1 ^ 2 + sigma2 ^ 2 - 2 * rho * sigma1 * sigma2)
    rho1 = (sigma1 - rho * sigma2) / v
    rho2 = (sigma2 - rho * sigma1) / v
    d = (log(S1 / S2) + (b1 - b2 + v ^ 2 / 2) * Time) / (v * sqrt(Time))
    y1 = (log(S1 / X) + (b1 + sigma1 ^ 2 / 2) * Time) / (sigma1 * sqrt(Time))
    y2 = (log(S2 / X) + (b2 + sigma2 ^ 2 / 2) * Time) / (sigma2 * sqrt(Time))

    # Calculate Price:
    OnTheMaxMin = NA
    if (TypeFlag == "cmin")
        OnTheMaxMin = (S1 * exp((b1 - r) * Time) * CBND(y1, -d, -rho1)
                      + S2 * exp((b2 - r) * Time) * CBND(y2, d - v *
                      sqrt(Time), -rho2) - X * exp(-r * Time) *
                      CBND(y1 - sigma1 * sqrt(Time), y2 - sigma2 *
                      sqrt(Time), rho))
    if (TypeFlag == "cmax")
        OnTheMaxMin = (S1 * exp((b1 - r) * Time) * CBND(y1, d, rho1) +
                      S2 * exp((b2 - r) * Time) * CBND(y2, -d + v *
                      sqrt(Time), rho2) - X * exp(-r * Time) * (1 -
                      CBND(-y1 + sigma1*sqrt(Time), -y2 + sigma2 *
                      sqrt(Time), rho)))
    if (TypeFlag == "pmin")
        OnTheMaxMin = (X * exp(-r * Time) - S1 * exp((b1 - r) * Time)
                       + EuropeanExchangeOption(S1, S2, 1, 1, Time, r,
                       b1, b2, sigma1, sigma2, rho)@price +
                       TwoRiskyAssetsOption("cmin", S1, S2, X, Time,
                       r, b1, b2, sigma1, sigma2, rho)@price)
    if (TypeFlag == "pmax")
        OnTheMaxMin = (X * exp(-r * Time) - S2 * exp((b2 - r) * Time)
                       - EuropeanExchangeOption(S1, S2, 1, 1, Time, r,
                       b1, b2, sigma1, sigma2, rho)@price +
                       TwoRiskyAssetsOption("cmax", S1, S2, X, Time,
                       r, b1, b2, sigma1, sigma2, rho)@price)

    # Parameters:
    # TypeFlag = c("cmin", "cmax", "pmin", "pmax"), S1, S2, X, Time, r,
    #   b1, b2, sigma1, sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$X = X
    param$Time = Time
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Two Risky Assets Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = OnTheMaxMin,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


SpreadApproxOption =
    function(TypeFlag = c("c", "p"), S1, S2, X, Time, r, sigma1, sigma2, rho,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Spread Option Approximation

    # References:
    #   Haug, Chapter 2.8.5

    # FUNCTION:

    # Compute Settings:
    TypeFlag = TypeFlag[1]
    F1 = S1
    F2 = S2
    sigma = sqrt(sigma1 ^ 2 + (sigma2 * F2 / (F2 + X)) ^ 2 - 2 * rho *
    sigma1 * sigma2 * F2 / (F2 + X))
    FF = F1 / (F2 + X)

    # Calculate Price
    SpreadApproximation <- (GBSOption(TypeFlag, FF, 1, Time, r, 0, sigma)@price *
                            (F2 + X)  * exp(-r * Time))

    # Parameters:
    # TypeFlag = c("c", "p"), S1, S2, X, Time, r, sigma1, sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$X = X
    param$Time = Time
    param$r = r
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Spread Approx Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = SpreadApproximation,
        title = title,
        description = description
        )
}


################################################################################

