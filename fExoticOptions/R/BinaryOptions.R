
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
# Binary Options:
#  GapOption                     Gap Option
#  CashOrNothingOption           Cash Or Nothing Option
#  TwoAssetCashOrNothingOption   Two Asset Cash-Or Nothing Option
#  AssetOrNothingOption          Asset Or Nothing Option
#  SuperShareOption              Super Share Option
#  BinaryBarrierOption           Binary Barrier Option
################################################################################


GapOption =
    function(TypeFlag = c("c", "p"), S, X1, X2, Time, r, b, sigma,
             title = NULL, description = NULL)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Gap Options

    # References:
    #   Haug, Haug Chapter 2.11.1

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    d1 = (log(S/X1) + (b + sigma^2 / 2) * Time) / (sigma * sqrt(Time))
    d2 = d1 - sigma*sqrt (Time)
    if (TypeFlag == "c")
        GapOption = S*exp((b-r)*Time)*CND(d1) - X2*exp(-r*Time)*CND(d2)
    if (TypeFlag == "p")
        GapOption = X2*exp(-r*Time)*CND(-d2) - S*exp((b-r)*Time)*CND(-d1)

    # Parameters:
    # TypeFlag = c("c", "p"), S, X1, X2, Time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X1 = X1
    param$X2 = X2
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Gap Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = GapOption,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


CashOrNothingOption =
    function(TypeFlag = c("c", "p"), S, X, K, Time, r, b, sigma,
             title = NULL, description = NULL)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Cash-Or-Nothing Options

    # References:
    #   Haug, Chapter 2.11.2

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    d = (log(S / X) + (r + b - sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time))
    if (TypeFlag == "c")
        CashOrNothing = K * exp (-r * Time) * CND(d)
    if (TypeFlag == "p")
        CashOrNothing = K * exp (-r * Time) * CND(-d)

    # Parameters:
    # TypeFlag = c("c", "p"), S, X, K, Time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$K = K
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Cash Or Nothing Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = CashOrNothing,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


TwoAssetCashOrNothingOption =
    function(TypeFlag = c("c", "p", "ud", "du"), S1, S2, X1, X2, K, Time, r,
             b1, b2, sigma1, sigma2, rho, title = NULL, description = NULL)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Two Asset Cash-Or-Nothing Options

    # References:
    #   Haug, Chapter 2.11.3

    # Arguments:
    #   1: Asset One, 2: Asset Two
    #   TypeFlag
    #       1   Call
    #       2   Put
    #       3   Up-Down
    #       4   Down-Up
    #   S=c(S1,S2)      Asset Prices
    #   K           Payout
    #   X=c(X1,X2)      Strikes
    #   b=c(b1,b2)      Cost-of-Carry
    #   sigma=c(sigma1,sigma2)  Volatilities
    #   rho         Correlation
    #

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    d11 = ((log(S1/X1) + (b1 - sigma1^2/2) * Time) /
           (sigma1*sqrt(Time)))
    d22 = ((log(S2/X2) + (b2 - sigma2^2/2) * Time) /
           (sigma2*sqrt(Time)))
    # Select:
    if (TypeFlag == "c")
        TwoAssetCashOrNothing = K * exp (-r * Time) * CBND( d11,  d22,  rho)
    if (TypeFlag == "p")
        TwoAssetCashOrNothing = K * exp (-r * Time) * CBND(-d11, -d22,  rho)
    if (TypeFlag == "ud")
        TwoAssetCashOrNothing = K * exp (-r * Time) * CBND( d11, -d22, -rho)
    if (TypeFlag == "du")
        TwoAssetCashOrNothing = K * exp (-r * Time) * CBND(-d11,  d22, -rho)

    # Parameters:
    # TypeFlag = c("c", "p", "ud", "du"), S1, S2, X1, X2, K, Time, r,
    #   b1, b2, sigma1, sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$X1 = X1
    param$X2 = X2
    param$K = K
    param$Time = Time
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Two Asset Cash Or Nothing Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = TwoAssetCashOrNothing,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


AssetOrNothingOption =
    function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma,
             title = NULL, description = NULL)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Asset-or-Nothing Options

    # Reference:
    #   Cox Rubinstein (1985)
    #   Haug, Chapter 2.11.4

    # FUNCTION:

    # Compute Price:
    TypeFlag = TypeFlag[1]
    d = (log(S/X) + (b + sigma^2 / 2) * Time) / (sigma * sqrt(Time))
    if (TypeFlag == "c")
        AssetOrNothing = S * exp ((b - r) * Time) * CND( d)
    if (TypeFlag == "p")
        AssetOrNothing = S * exp ((b - r) * Time) * CND(-d)

    # Parameters:
    # # TypeFlag = c("c", "p"), S, X, Time, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Asset Or Nothing Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = AssetOrNothing,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


SuperShareOption =
    function(S, XL, XH, Time, r, b, sigma, title = NULL, description = NULL)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Supershare Options

    # Reference:
    #   Hakansson (1976)
    #   Haug, Chapter 2.11.5

    # FUNCTION:

    # Compute Price:
    d1 = (log(S/XL) + (b + sigma^2 / 2) * Time) / (sigma * sqrt(Time))
    d2 = (log(S/XH) + (b + sigma^2 / 2) * Time) / (sigma * sqrt(Time))
    SuperShare = (S * exp((b-r)*Time) / XL) * (CND(d1) - CND(d2))

    # Parameters:
    # S, XL, XH, Time, r, b, sigma
    param = list()
    param$S = S
    param$XL = XL
    param$XH = XH
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Super Share Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = SuperShare,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


BinaryBarrierOption =
    function(TypeFlag = as.character(1:28), S, X, H, K, Time, r, b, sigma,
             eta, phi, title = NULL, description = NULL)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Binary Barrier Options

    # Reference:
    #   Reiner and Rubinstein (1991)
    #   Haug, Chapter 2.11.6

    # FUNCTION:

    # Compute Price:
    TypeFlag = as.integer(TypeFlag[1])
    eta = rep(c(+1,-1), 14)[TypeFlag]
    #         1  2  3  4  5  6  7  8  9 10 11 12 13 14
    #        15 16 17 18 19 20 21 22 23 24 25 26 27 28
    phi = c(+0,+0,+0,+0,-1,+1,-1,+1,+1,-1,+1,-1,+1,+1,
            +1,+1,-1,-1,-1,-1,+1,+1,+1,+1,-1,-1,-1,-1)[TypeFlag]
    v = sigma
    mu = (b - v ^ 2 / 2) / v ^ 2
    lambda = sqrt(mu ^ 2 + 2 * r / v ^ 2)
    X1 = log(S / X) / (v * sqrt(Time)) + (mu + 1) * v * sqrt(Time)
    X2 = log(S / H) / (v * sqrt(Time)) + (mu + 1) * v * sqrt(Time)
    y1 = log(H ^ 2 / (S * X)) / (v * sqrt(Time)) + (mu + 1) * v * sqrt(Time)
    y2 = log(H / S) / (v * sqrt(Time)) + (mu + 1) * v * sqrt(Time)
    Z = log(H / S) / (v * sqrt(Time)) + lambda * v * sqrt(Time)

    # Values:
    a1 = S * exp((b - r) * Time) * CND(phi * X1)
    b1 = K * exp(-r * Time) * CND(phi * X1 - phi * v * sqrt(Time))
    a2 = S * exp((b - r) * Time) * CND(phi * X2)
    b2 = K * exp(-r * Time) * CND(phi * X2 - phi * v * sqrt(Time))
    a3 = (S * exp((b - r) * Time) * (H / S) ^ (2 * (mu + 1)) *
          CND(eta * y1))
    b3 = (K * exp(-r * Time) * (H / S) ^ (2 * mu) *
          CND(eta * y1 - eta * v * sqrt(Time)))
    a4 = (S * exp((b - r) * Time) * (H / S) ^ (2 * (mu + 1)) *
          CND(eta * y2))
    b4 = (K * exp(-r * Time) * (H / S) ^ (2 * mu) *
          CND(eta * y2 - eta * v * sqrt(Time)))
    a5 = (K * ((H / S) ^ (mu + lambda) *
               CND(eta * Z) + (H / S) ^ (mu - lambda) *
               CND(eta * Z - 2 * eta * lambda * v * sqrt(Time))))
    # Select:
    BinaryBarrier = NA
    if (X > H) {
        if (TypeFlag ==  1) BinaryBarrier = a5
        if (TypeFlag ==  2) BinaryBarrier = a5
        if (TypeFlag ==  3) BinaryBarrier = a5
        if (TypeFlag ==  4) BinaryBarrier = a5
        if (TypeFlag ==  5) BinaryBarrier = b2 + b4
        if (TypeFlag ==  6) BinaryBarrier = b2 + b4
        if (TypeFlag ==  7) BinaryBarrier = a2 + a4
        if (TypeFlag ==  8) BinaryBarrier = a2 + a4
        if (TypeFlag ==  9) BinaryBarrier = b2 - b4
        if (TypeFlag == 10) BinaryBarrier = b2 - b4
        if (TypeFlag == 11) BinaryBarrier = a2 - a4
        if (TypeFlag == 12) BinaryBarrier = a2 - a4
        if (TypeFlag == 13) BinaryBarrier = b3
        if (TypeFlag == 14) BinaryBarrier = b3
        if (TypeFlag == 15) BinaryBarrier = a3
        if (TypeFlag == 16) BinaryBarrier = a1
        if (TypeFlag == 17) BinaryBarrier = b2 - b3 + b4
        if (TypeFlag == 18) BinaryBarrier = b1 - b2 + b4
        if (TypeFlag == 19) BinaryBarrier = a2 - a3 + a4
        if (TypeFlag == 20) BinaryBarrier = a1 - a2 + a3
        if (TypeFlag == 21) BinaryBarrier = b1 - b3
        if (TypeFlag == 22) BinaryBarrier = 0
        if (TypeFlag == 23) BinaryBarrier = a1 - a3
        if (TypeFlag == 24) BinaryBarrier = 0
        if (TypeFlag == 25) BinaryBarrier = b1 - b2 + b3 - b4
        if (TypeFlag == 26) BinaryBarrier = b2 - b4
        if (TypeFlag == 27) BinaryBarrier = a1 - a2 + a3 - a4
        if (TypeFlag == 28) BinaryBarrier = a2 - a4 }
    # Continue:
    if (X < H) {
        if (TypeFlag ==  1) BinaryBarrier = a5
        if (TypeFlag ==  2) BinaryBarrier = a5
        if (TypeFlag ==  3) BinaryBarrier = a5
        if (TypeFlag ==  4) BinaryBarrier = a5
        if (TypeFlag ==  5) BinaryBarrier = b2 + b4
        if (TypeFlag ==  6) BinaryBarrier = b2 + b4
        if (TypeFlag ==  7) BinaryBarrier = a2 + a4
        if (TypeFlag ==  8) BinaryBarrier = a2 + a4
        if (TypeFlag ==  9) BinaryBarrier = b2 - b4
        if (TypeFlag == 10) BinaryBarrier = b2 - b4
        if (TypeFlag == 11) BinaryBarrier = a2 - a4
        if (TypeFlag == 12) BinaryBarrier = a2 - a4
        if (TypeFlag == 13) BinaryBarrier = b1 - b2 + b4
        if (TypeFlag == 14) BinaryBarrier = b2 - b3 + b4
        if (TypeFlag == 15) BinaryBarrier = a1 - a2 + a4
        if (TypeFlag == 16) BinaryBarrier = a2 - a3 + a4
        if (TypeFlag == 17) BinaryBarrier = b1
        if (TypeFlag == 18) BinaryBarrier = b3
        if (TypeFlag == 19) BinaryBarrier = a1
        if (TypeFlag == 20) BinaryBarrier = a3
        if (TypeFlag == 21) BinaryBarrier = b2 - b4
        if (TypeFlag == 22) BinaryBarrier = b1 - b2 + b3 - b4
        if (TypeFlag == 23) BinaryBarrier = a2 - a4
        if (TypeFlag == 24) BinaryBarrier = a1 - a2 + a3 - a4
        if (TypeFlag == 25) BinaryBarrier = 0
        if (TypeFlag == 26) BinaryBarrier = b1 - b3
        if (TypeFlag == 27) BinaryBarrier = 0
        if (TypeFlag == 28) BinaryBarrier = a1 - a3 }

    # Parameters:
    # TypeFlag = as.character(1:28), S, X, H, K, Time, r, b, sigma, eta, phi
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$H = H
    param$K = K
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$eta = eta
    param$phi = phi

    # Add title and description:
    if (is.null(title)) title = "Binary Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = BinaryBarrier,
        title = title,
        description = description
        )
}


################################################################################

