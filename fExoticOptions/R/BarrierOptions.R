
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
# Barrier Options:
#   StandardBarrierOption         Standard Barrier Option
#   DoubleBarrierOption           Double Barrier Option
#   PTSingleAssetBarrierOption    Partial Time Barrier Option
#   TwoAssetBarrierOption         Two Asset Barrier
#   PTTwoAssetBarrierOption       Partial Time TwoAsset Barrier Option
#   LookBarrierOption             Look Barrier Option
#   DiscreteBarrierOption         Discrete Adjusted Barrier Option
#   SoftBarrierOption             Soft Barrier Option
################################################################################


StandardBarrierOption =
    function(TypeFlag = c("cdi", "cui", "pdi", "pui", "cdo", "cuo", "pdo", "puo"),
             S, X, H, K, Time, r, b, sigma, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Standard Barrier Options

    # References:
    #   Haug, Chapter 2.10.1

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    StandardBarrier = NA
    mu = (b - sigma ^ 2 / 2) / sigma ^ 2
    lambda = sqrt (mu ^ 2 + 2 * r / sigma ^ 2)
    X1 = log (S / X) / (sigma * sqrt(Time)) + (1 + mu) * sigma * sqrt(Time)
    X2 = log (S / H) / (sigma * sqrt(Time)) + (1 + mu) * sigma * sqrt(Time)
    y1 = (log (H ^ 2 / (S * X)) / (sigma * sqrt(Time)) + (1 + mu) * sigma *
          sqrt(Time))
    y2 = log (H / S) / (sigma * sqrt(Time)) + (1 + mu) * sigma * sqrt(Time)
    Z  = log (H / S) / (sigma * sqrt(Time)) + lambda * sigma * sqrt(Time)
    if (TypeFlag == "cdi" || TypeFlag == "cdo") { eta = +1; phi = +1 }
    if (TypeFlag == "cui" || TypeFlag == "cuo") { eta = -1; phi = +1 }
    if (TypeFlag == "pdi" || TypeFlag == "pdo") { eta = +1; phi = -1 }
    if (TypeFlag == "pui" || TypeFlag == "puo") { eta = -1; phi = -1 }
    f1 = (phi * S * exp ((b - r) * Time) * CND(phi * X1) -
          phi * X * exp(-r * Time) * CND(phi * X1 - phi * sigma * sqrt(Time)))
    f2 = (phi * S * exp ((b - r) * Time) * CND(phi * X2) -
          phi * X * exp(-r * Time) * CND(phi * X2 - phi * sigma * sqrt(Time)))
    f3 = (phi * S * exp ((b - r) * Time) * (H / S) ^ (2 * (mu + 1))  *
          CND(eta * y1) - phi * X * exp(-r * Time) * (H / S) ^ (2 * mu) *
          CND(eta * y1 - eta * sigma * sqrt(Time)))
    f4 = (phi * S * exp ((b - r) * Time) * (H / S) ^ (2 * (mu + 1)) *
          CND(eta * y2) - phi * X * exp(-r * Time) * (H / S) ^ (2 * mu) *
          CND(eta * y2 - eta * sigma * sqrt(Time)))
    f5 = (K * exp (-r * Time) * (CND(eta * X2 - eta * sigma *
          sqrt(Time)) - (H / S) ^ (2 * mu) * CND(eta * y2 - eta *
          sigma * sqrt(Time))))
    f6 = (K * ((H / S) ^ (mu + lambda) * CND(eta * Z) + (H / S)^(mu - lambda) *
          CND(eta * Z - 2 * eta * lambda * sigma * sqrt(Time))))
    if (X >= H) {
        if (TypeFlag == "cdi") StandardBarrier = f3 + f5
        if (TypeFlag == "cui") StandardBarrier = f1 + f5
        if (TypeFlag == "pdi") StandardBarrier = f2 - f3 + f4 + f5
        if (TypeFlag == "pui") StandardBarrier = f1 - f2 + f4 + f5
        if (TypeFlag == "cdo") StandardBarrier = f1 - f3 + f6
        if (TypeFlag == "cuo") StandardBarrier = f6
        if (TypeFlag == "pdo") StandardBarrier = f1 - f2 + f3 - f4 + f6
        if (TypeFlag == "puo") StandardBarrier = f2 - f4 + f6 }
    if (X < H) {
        if (TypeFlag == "cdi") StandardBarrier = f1 - f2 + f4 + f5
        if (TypeFlag == "cui") StandardBarrier = f2 - f3 + f4 + f5
        if (TypeFlag == "pdi") StandardBarrier = f1 + f5
        if (TypeFlag == "pui") StandardBarrier = f3 + f5
        if (TypeFlag == "cdo") StandardBarrier = f2 + f6 - f4
        if (TypeFlag == "cuo") StandardBarrier = f1 - f2 + f3 - f4 + f6
        if (TypeFlag == "pdo") StandardBarrier = f6
        if (TypeFlag == "puo") StandardBarrier = f1 - f3 + f6 }

    # Parameters:
    # TypeFlag = c("cdi", "cui", "pdi", "pui", "cdo", "cuo", "pdo", "puo"),
    #   S, X, H, K, Time, r, b, sigma
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

    # Add title and description:
    if (is.null(title)) title = "Standard Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = StandardBarrier,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


DoubleBarrierOption =
    function(TypeFlag = c("co", "ci", "po", "pi"), S, X, L, U, Time, r, b,
             sigma, delta1, delta2, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Double barrier options

    # References:
    #   Haug, Chapter 2.10.2

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    DoubleBarrier = NA
    FU = U * exp (delta1 * Time)
    E = L * exp (delta1 * Time)
    Sum1 = Sum2 = 0

    # Call:
    if (TypeFlag == "co" || TypeFlag == "ci") {
        for (n in -5:5) {
            d1 = ((log(S * U ^ (2 * n) / (X * L ^ (2 * n))) +
                   (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            d2 = ((log(S * U ^ (2 * n) / (FU * L ^ (2 * n))) +
                  (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            d3 = ((log(L ^ (2 * n + 2) / (X * S * U ^ (2 * n))) +
                  (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            d4 = ((log(L ^ (2 * n + 2) / (FU * S * U ^ (2 * n))) +
                  (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / sigma^2 + 1
            mu2 = 2 * n * (delta1 - delta2) / sigma^2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / sigma^2 + 1
            Sum1 = (Sum1 + (U^n / L ^ n) ^ mu1 * (L / S) ^ mu2 *
                    (CND(d1) - CND(d2)) - (L^(n + 1) / (U ^ n * S)) ^ mu3 *
                    (CND(d3) - CND(d4)))
            Sum2 = (Sum2 + (U^n / L ^ n) ^ (mu1 - 2) * (L/S)^mu2 *
                    (CND(d1 - sigma * sqrt(Time)) - CND(d2 - sigma * sqrt(Time))) -
                    (L^(n + 1) / (U ^ n * S))^(mu3 - 2) *
                    (CND(d3 - sigma * sqrt(Time)) - CND(d4 - sigma * sqrt(Time))))
        }
        OutValue = S * exp ((b-r)*Time) * Sum1 - X * exp(-r*Time) * Sum2 }

    # Put:
    if (TypeFlag == "po" || TypeFlag == "pi") {
        for ( n in (-5:5) ) {
            d1 = ((log(S * U ^ (2 * n) / (E * L ^ (2 * n))) +
                   (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            d2 = ((log(S * U ^ (2 * n) / (X * L ^ (2 * n))) +
                   (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            d3 = ((log(L ^ (2 * n + 2) / (E * S * U ^ (2 * n))) +
                   (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            d4 = ((log(L ^ (2 * n + 2) / (X * S * U ^ (2 * n))) +
                   (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time)))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / sigma ^ 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / sigma ^ 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / sigma ^ 2 + 1
            Sum1 = (Sum1 + (U^n / L^n)^mu1 * (L / S) ^ mu2 *
                    (CND(d1) - CND(d2)) -
                    (L ^ (n + 1) / (U ^ n * S)) ^ mu3 *
                    (CND(d3) - CND(d4)))
            Sum2 = (Sum2 + (U ^n / L^n)^(mu1 - 2) * (L/S)^mu2 *
                    (CND(d1 - sigma * sqrt(Time)) - CND(d2 - sigma * sqrt(Time))) -
                    (L^(n + 1) / (U ^ n * S))^(mu3 - 2) *
                    (CND(d3 - sigma * sqrt(Time)) - CND(d4 - sigma * sqrt(Time))))
        }
        OutValue = X * exp (-r*Time) * Sum2 - S * exp((b - r)*Time) * Sum1 }

    # Final Values:
    if (TypeFlag == "co" || TypeFlag == "po")
        DoubleBarrier = OutValue
    if (TypeFlag == "ci")
        DoubleBarrier = GBSOption("c", S, X, Time, r, b, sigma)@price - OutValue
    if (TypeFlag == "pi")
        DoubleBarrier = GBSOption("p", S, X, Time, r, b, sigma)@price - OutValue

    # Parameters:
    # TypeFlag = c("co", "ci", "po", "pi"), S, X, L, U, Time, r, b,
    #   sigma, delta1, delta2
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$L = L
    param$U = U
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$delta1 = delta1
    param$delta2 = delta2

    # Add title and description:
    if (is.null(title)) title = "Double Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = DoubleBarrier,
        title = title,
        description = description
        )
}



# ------------------------------------------------------------------------------


PTSingleAssetBarrierOption =
    function(TypeFlag = c("cdoA", "cuoA", "pdoA", "puoA", "coB1", "poB1",
             "cdoB2", "cuoB2"), S, X, H, time1, Time2, r, b, sigma, title = NULL,
             description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Partial-time single asset barrier options

    # References:
    #   Haug, Chapter 2.10.3

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    PartialTimeBarrier = NA
    t1 = time1
    T2 = Time2
    if (TypeFlag == "cdoA") eta = 1
    if (TypeFlag == "cuoA") eta = -1

    # Continue:
    d1 = (log(S/X) + (b + sigma^2/2) * T2) / (sigma * sqrt(T2))
    d2 = d1 - sigma * sqrt (T2)
    f1 = (log(S/X) + 2 * log(H/S) + (b + sigma^2/2) * T2) / (sigma * sqrt(T2))
    f2 = f1 - sigma * sqrt (T2)
    e1 = (log(S / H) + (b + sigma ^ 2 / 2) * t1) / (sigma * sqrt(t1))
    e2 = e1 - sigma * sqrt (t1)
    e3 = e1 + 2 * log (H / S) / (sigma * sqrt(t1))
    e4 = e3 - sigma * sqrt (t1)
    mu = (b - sigma ^ 2 / 2) / sigma ^ 2
    rho = sqrt (t1 / T2)
    g1 = (log(S / H) + (b + sigma ^ 2 / 2) * T2) / (sigma * sqrt(T2))
    g2 = g1 - sigma * sqrt (T2)
    g3 = g1 + 2 * log (H / S) / (sigma * sqrt(T2))
    g4 = g3 - sigma * sqrt (T2)
    z1 = CND (e2) - (H / S) ^ (2 * mu) * CND(e4)
    z2 = CND (-e2) - (H / S) ^ (2 * mu) * CND(-e4)
    z3 = CBND (g2, e2, rho) - (H / S) ^ (2 * mu) * CBND(g4, -e4, -rho)
    z4 = CBND (-g2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-g4, e4, -rho)
    z5 = CND (e1) - (H / S) ^ (2 * (mu + 1)) * CND(e3)
    z6 = CND (-e1) - (H / S) ^ (2 * (mu + 1)) * CND(-e3)
    z7 = CBND (g1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(g3, -e3, -rho)
    z8 = CBND (-g1, -e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(-g3, e3, -rho)

    if (TypeFlag == "cdoA" || TypeFlag == "cuoA") {
        # call down-and out and up-and-out type A
        PartialTimeBarrier =
            (S * exp ((b - r) * T2) * (CBND(d1, eta * e1, eta * rho) -
             (H / S) ^ (2 * (mu + 1)) * CBND(f1, eta * e3, eta * rho))
             - X * exp (-r * T2) * (CBND(d2, eta * e2, eta * rho) - (H
             / S) ^ (2 * mu) * CBND(f2, eta * e4, eta * rho))) }

    if (TypeFlag == "cdoB2" && X < H) {
        # call down-and-out type B2
        PartialTimeBarrier =
            (S * exp ((b - r) * T2) * (CBND(g1, e1, rho) - (H / S) ^
            (2 * (mu + 1)) * CBND(g3, -e3, -rho)) - X * exp (-r * T2)
            * (CBND(g2, e2, rho) - (H / S) ^ (2 * mu) * CBND(g4, -e4,
                                                             -rho))) }

    if (TypeFlag == "cdoB2" && X > H) {
        PartialTimeBarrier = (PTSingleAssetBarrierOption("coB1", S, X, H, t1,
                                                         T2, r, b, sigma)@price) }

    if (TypeFlag == "cuoB2" && X < H) {
        # call up-and-out type B2
        PartialTimeBarrier =
            (S * exp ((b - r) * T2) * (CBND(-g1, -e1, rho) - (H / S) ^
            (2 * (mu + 1)) * CBND(-g3, e3, -rho)) - X * exp (-r * T2)
            * (CBND(-g2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-g4,
            e4, -rho)) - S * exp ((b - r) * T2) * (CBND(-d1, -e1, rho)
            - (H / S) ^ (2 * (mu + 1)) * CBND(e3, -f1, -rho)) + X *
            exp (-r * T2) * (CBND(-d2, -e2, rho) - (H / S) ^ (2 * mu)
            * CBND(e4, -f2, -rho)))}

    if (TypeFlag == "coB1" && X > H) {
        # call out type B1
        PartialTimeBarrier =
            (S * exp ((b - r) * T2) * (CBND(d1, e1, rho) - (H / S) ^
            (2 * (mu + 1)) * CBND(f1, -e3, -rho)) - X * exp (-r * T2)
            * (CBND(d2, e2, rho) - (H / S) ^ (2 * mu) * CBND(f2, -e4,
            -rho))) }

    if (TypeFlag == "coB1" && X < H) {
        PartialTimeBarrier =
            (S * exp ((b - r) * T2) * (CBND(-g1, -e1, rho) - (H / S) ^
            (2 * (mu + 1)) * CBND(-g3, e3, -rho)) - X * exp (-r * T2)
            * (CBND(-g2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-g4,
            e4, -rho)) - S * exp ((b - r) * T2) * (CBND(-d1, -e1, rho)
            - (H / S) ^ (2 * (mu + 1)) * CBND(-f1, e3, -rho)) + X *
            exp (-r * T2) * (CBND(-d2, -e2, rho) - (H / S) ^ (2 * mu)
            * CBND(-f2, e4, -rho)) + S * exp ((b - r) * T2) *
            (CBND(g1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(g3,
            -e3, -rho)) - X * exp (-r * T2) * (CBND(g2, e2, rho) - (H
            / S) ^ (2 * mu) * CBND(g4, -e4, -rho))) }

    if (TypeFlag == "pdoA") {
        # put down-and out and up-and-out type A
        PartialTimeBarrier = (PTSingleAssetBarrierOption("cdoA",
                              S, X, H, t1, T2, r, b, sigma)@price -
                              S * exp((b - r) * T2) * z5 + X * exp(-r * T2) * z1)}

    if (TypeFlag == "puoA") {
        PartialTimeBarrier = (PTSingleAssetBarrierOption("cuoA",
                              S, X, H, t1, T2, r, b, sigma)@price -
                              S * exp((b - r) * T2) * z6 + X * exp(-r * T2) * z2) }

    if (TypeFlag == "poB1") {
        # put out type B1
        PartialTimeBarrier = (PTSingleAssetBarrierOption("coB1",
                              S, X, H, t1, T2, r, b, sigma)@price -
                              S * exp((b - r) * T2) * z8 + X * exp(-r * T2) * z4 -
                              S * exp((b - r) * T2) * z7 + X * exp(-r * T2) * z3) }

    if (TypeFlag == "pdoB2") {
        # put down-and-out type B2
        PartialTimeBarrier = (PTSingleAssetBarrierOption("cdoB2",
                              S, X, H, t1, T2, r, b, sigma)@price -
                              S * exp((b - r) * T2) * z7 + X * exp(-r * T2) * z3) }

    if (TypeFlag == "puoB2") {
        # put up-and-out type B2
        PartialTimeBarrier = (PTSingleAssetBarrierOption("cuoB2",
                              S, X, H, t1, T2, r, b, sigma)@price -
                              S * exp((b - r) * T2) * z8 + X * exp(-r * T2) * z4) }

    # Parameters:
    # TypeFlag = c("cdoA", "cuoA", "pdoA", "puoA", "coB1", "poB1",
    #   "cdoB2", "cuoB2"), S, X, H, time1, Time2, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$H = H
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Partial Time Single Asset Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = PartialTimeBarrier,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


TwoAssetBarrierOption =
    function(TypeFlag = c("cuo", "cui", "cdo", "cdi", "puo", "pui", "pdo", "pdi"),
             S1, S2, X, H, Time, r, b1, b2, sigma1, sigma2, rho, title = NULL,
             description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Partial-time singel asset barrier options

    # References:
    #   Haug, Chapter 2.10.4

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    v1 = sigma1
    v2 = sigma2
    mu1 = b1 - v1 ^ 2 / 2
    mu2 = b2 - v2 ^ 2 / 2
    d1 = (log(S1 / X) + (mu1 + v1 ^ 2 / 2) * Time) / (v1 * sqrt(Time))
    d2 = d1 - v1 * sqrt (Time)
    d3 = d1 + 2 * rho * log (H / S2) / (v2 * sqrt(Time))
    d4 = d2 + 2 * rho * log (H / S2) / (v2 * sqrt(Time))
    e1 = (log(H / S2) - (mu2 + rho * v1 * v2) * Time) / (v2 * sqrt(Time))
    e2 = e1 + rho * v1 * sqrt (Time)
    e3 = e1 - 2 * log (H / S2) / (v2 * sqrt(Time))
    e4 = e2 - 2 * log (H / S2) / (v2 * sqrt(Time))

    # Make Decisions:
    if (TypeFlag == "cuo" || TypeFlag == "cui") {
        eta = 1
        phi = 1 }
    if (TypeFlag == "cdo" || TypeFlag == "cdi") {
        eta = 1
        phi = -1 }
    if (TypeFlag == "puo" || TypeFlag == "pui") {
        eta = -1
        phi = 1 }
    if (TypeFlag == "pdo" || TypeFlag == "pdi") {
        eta = -1
        phi = -1 }

    # Calculate Knock Out Value:
    KnockOutValue =
        (eta * S1 * exp ((b1 - r) * Time) * (CBND ( eta * d1, phi *
         e1, -eta * phi * rho) - exp (2 * (mu2 + rho * v1 * v2) *
         log(H / S2) / v2 ^ 2) * CBND(eta * d3, phi * e3, -eta * phi *
         rho)) - eta * exp(-r * Time) * X * (CBND(eta * d2, phi * e2,
         -eta * phi * rho) - exp (2 * mu2 * log(H / S2) / v2 ^ 2) *
         CBND(eta * d4, phi * e4, -eta * phi * rho)))

    # Calculate Two Asset Barrier:
    if (TypeFlag == "cuo" || TypeFlag == "cdo" ||
        TypeFlag == "puo" || TypeFlag == "pdo")
        TwoAssetBarrier = KnockOutValue
    if (TypeFlag == "cui" || TypeFlag == "cdi")
        TwoAssetBarrier = (GBSOption("c", S1, X, Time, r, b1, v1)@price -
                           KnockOutValue)
    if (TypeFlag == "pui" || TypeFlag == "pdi")
        TwoAssetBarrier = (GBSOption("p", S1, X, Time, r, b1, v1)@price -
                           KnockOutValue)

    # Parameters:
    # TypeFlag = c("cuo", "cui", "cdo", "cdi", "puo", "pui", "pdo", "pdi"),
    #   S1, S2, X, H, Time, r, b1, b2, sigma1, sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$X = X
    param$H = H
    param$Time = Time
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho = rho

    # Add title and description:
    if (is.null(title)) title = "Two Asset Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = TwoAssetBarrier,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


PTTwoAssetBarrierOption =
    function(TypeFlag = c("cdo", "pdo", "cdi", "pdi", "cuo", "puo", "cui", "pui"),
             S1, S2, X, H, time1, Time2, r, b1, b2, sigma1, sigma2, rho, title = NULL,
             description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Partial-time two asset barrier options

    # References:
    #   Haug, Chapter 2.10.5

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    t1 = time1
    T2 = Time2
    v1 = sigma1
    v2 = sigma2
    if (TypeFlag == "cdo" || TypeFlag == "pdo" ||
        TypeFlag == "cdi" || TypeFlag == "pdi") {
        phi = -1 }
    else {
        phi = 1 }
    if (TypeFlag == "cdo" || TypeFlag == "cuo" ||
        TypeFlag == "cdi" || TypeFlag == "cui") {
        eta = 1 }
    else {
        eta = -1 }

    mu1 = b1 - v1 ^ 2 / 2
    mu2 = b2 - v2 ^ 2 / 2
    d1 = (log(S1 / X) + (mu1 + v1 ^ 2) * T2) / (v1 * sqrt(T2))
    d2 = d1 - v1 * sqrt (T2)
    d3 = d1 + 2 * rho * log (H / S2) / (v2 * sqrt(T2))
    d4 = d2 + 2 * rho * log (H / S2) / (v2 * sqrt(T2))
    e1 = (log(H / S2) - (mu2 + rho * v1 * v2) * t1) / (v2 * sqrt(t1))
    e2 = e1 + rho * v1 * sqrt (t1)
    e3 = e1 - 2 * log (H / S2) / (v2 * sqrt(t1))
    e4 = e2 - 2 * log (H / S2) / (v2 * sqrt(t1))

    OutBarrierValue =
        (eta * S1 * exp ((b1 - r) * T2) * (CBND(eta * d1, phi * e1,
         -eta * phi * rho * sqrt(t1 / T2)) - exp(2 * log(H / S2) *
         (mu2 + rho * v1 * v2) / (v2 ^ 2)) * CBND (eta * d3, phi * e3,
         -eta * phi * rho * sqrt(t1 / T2))) - eta * exp (-r * T2) * X
         * (CBND(eta * d2, phi * e2, -eta * phi * rho * sqrt(t1 / T2))
         - exp(2 * log(H / S2) * mu2 / (v2 ^ 2)) * CBND (eta * d4, phi
         * e4, -eta * phi * rho * sqrt(t1 / T2))))

    if (TypeFlag == "cdo" || TypeFlag == "cuo" ||
        TypeFlag == "pdo" || TypeFlag == "puo")
        PartialTimeTwoAssetBarrier = OutBarrierValue
    if (TypeFlag == "cui" || TypeFlag == "cdi")
        PartialTimeTwoAssetBarrier = (GBSOption("c", S1, X, T2, r, b1, v1)@price -
                                      OutBarrierValue)
    if (TypeFlag == "pui" || TypeFlag == "pdi")
        PartialTimeTwoAssetBarrier = (GBSOption("p", S1, X, T2, r, b1, v1)@price -
                                      OutBarrierValue)

    # Parameters:
    # TypeFlag = c("cdo", "pdo", "cdi", "pdi", "cuo", "puo", "cui", "pui"),
    #   S1, S2, X, H, time1, Time2, r, b1, b2, sigma1, sigma2, rho
    param = list()
    param$TypeFlag = TypeFlag
    param$S1 = S1
    param$S2 = S2
    param$X = X
    param$H = H
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b1 = b1
    param$b2 = b2
    param$sigma1 = sigma1
    param$sigma2 = sigma2
    param$rho

    # Add title and description:
    if (is.null(title)) title = "Partial Time Two Asset Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = PartialTimeTwoAssetBarrier,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


LookBarrierOption =
    function(TypeFlag = c("cuo", "cui", "pdo", "pdi"), S, X, H, time1, Time2,
             r, b, sigma, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Look-Barrier Options

    # References:
    #   Haug, Chapter 2.10.6

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    t1 = time1
    T2 = Time2
    # Take care of the limit t1 -> T2
    if (T2 == t1) t1 = t1*(1-1.0e-12)
    v = sigma
    hh = log(H/S)
    K = log(X/S)
    mu1 = b - sigma^2/2
    mu2 = b + sigma^2/2
    rho = sqrt(t1/T2)

    # Make Decisions - Settings:
    if (TypeFlag == "cuo" || TypeFlag == "cui") {
        eta = +1
        m = min (hh, K) }
    if (TypeFlag == "pdo" || TypeFlag == "pdi") {
        eta = -1
        m = max (hh, K) }

    # Compute the g Values:
    g1 = ((CND(eta * (hh - mu2 * t1) / (sigma * sqrt(t1))) - exp(2 *
          mu2 * hh / sigma ^ 2) * CND(eta * (-hh - mu2 * t1) / (sigma
          * sqrt(t1)))) - (CND(eta * (m - mu2 * t1) / (sigma *
          sqrt(t1))) - exp(2 * mu2 * hh / sigma ^ 2) * CND(eta * (m -
          2 * hh - mu2 * t1) / (sigma * sqrt(t1)))))
    g2 = ((CND(eta * (hh - mu1 * t1) / (sigma * sqrt(t1))) - exp(2 *
          mu1 * hh / sigma ^ 2) * CND(eta * (-hh - mu1 * t1) / (sigma
          * sqrt(t1)))) - (CND(eta * (m - mu1 * t1) / (sigma *
          sqrt(t1))) - exp(2 * mu1 * hh / sigma ^ 2) * CND(eta * (m -
          2 * hh - mu1 * t1) / (sigma * sqrt(t1)))))

    # Needed by Out Value:
    part1 = (S * exp((b-r)*T2) * (1+v^2/(2*b)) * (CBND(
            eta*(+m-mu2*t1)/(v*sqrt(t1)),
            eta*(-K+mu2*T2)/(v*sqrt(T2)), -rho) - exp(2*mu2*hh/v^2) *
            CBND( eta*(m-2*hh-mu2*t1)/(v*sqrt(t1)),
            eta*(2*hh-K+mu2*T2)/(v*sqrt(T2)), -rho) ))
    part2 = (- X * exp(-r*T2) * ( CBND( eta*(+m-mu1*t1)/(v*sqrt(t1)),
            eta*(-K+mu1*T2)/(v*sqrt(T2)), -rho) - exp(2*mu1*hh/v^2) *
            CBND( eta*(m-2*hh-mu1*t1)/(v*sqrt(t1)),
            eta*(2*hh-K+mu1*T2)/(v*sqrt(T2)), -rho) ))

    part3 = (-exp(-r*T2) * v^2/(2*b) * ( S*(S/X)^(-2*b/v^2) * CBND(
             eta * (m + mu1 * t1) / (v * sqrt(t1)), eta * (-K - mu1 *
             T2) / (v * sqrt(T2)), -rho) - H*(H/X)^(-2*b/v^2) * CBND(
             eta*(m - 2 * hh + mu1 * t1) / (v * sqrt(t1)), eta * (2 *
             hh - K - mu1 * T2) / (v * sqrt(T2)), -rho) ))
    part4 = (S * exp((b-r)*T2) * ((1+v^2/(2 * b)) *
            CND(eta*mu2*(T2-t1)/(v*sqrt(T2-t1))) +
            exp(-b*(T2-t1))*(1-v^2/(2*b)) *
            CND(eta*(-mu1*(T2-t1))/(v*sqrt(T2-t1))))*g1 -
            exp(-r*T2)*X*g2)

    # Calculate Out Value:
    OutValue = eta * (part1 + part2 + part3 + part4)

    # Option Price:
    if (TypeFlag == "cuo" || TypeFlag == "pdo")
        LookBarrier = OutValue
    if (TypeFlag == "cui")
        LookBarrier = (PTFixedStrikeLookbackOption("c", S, X, t1, T2, r, b,
                                                   sigma)@price - OutValue)
    if (TypeFlag == "pdi")
        LookBarrier = (PTFixedStrikeLookbackOption("p", S, X, t1, T2, r, b,
                                                   sigma)@price - OutValue)

    # Parameters:
    # TypeFlag = c("cuo", "cui", "pdo", "pdi"), S, X, H, time1, Time2,
    #   r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$H = H
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Look Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = LookBarrier,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


DiscreteBarrierOption =
    function(S, H, sigma, dt, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Discrete Barrier Options

    # References:
    #   Haug, Chapter 2.10.

    # FUNCTION:

    # Compute:
    DiscreteBarrier = NA
    if (H > S) {
        DiscreteBarrier = H * exp(0.5826 * sigma * sqrt(dt)) }
    if (H < S) {
        DiscreteBarrier = H * exp(-0.5826 * sigma * sqrt(dt)) }

    # Parameters:
    # S, H, sigma, dt
    param = list()
    param$S = S
    param$H = H
    param$sigma = sigma
    param$dt = dt

    # Add title and description:
    if (is.null(title)) title = "Discrete Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = DiscreteBarrier,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


SoftBarrierOption =
    function(TypeFlag = c("cdi", "cdo", "pdi", "pdo"), S, X, L, U, Time ,
             r, b, sigma, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Soft Barrier Option

    # References:
    #   Haug, Haug Chapter 2.10.8


    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    v = sigma

    # Make Decisions - Settings:
    if (TypeFlag == "cdi" || TypeFlag == "cdo") {
        eta = 1}
    else {
        eta = -1 }

    # Continue:
    mu = (b + v ^ 2 / 2) / v ^ 2
    lambda1 = exp(-1 / 2 * v ^ 2 * Time * (mu + 0.5) * (mu - 0.5))
    lambda2 = exp(-1 / 2 * v ^ 2 * Time * (mu - 0.5) * (mu - 1.5))
    d1 = log(U ^ 2 / (S * X)) / (v * sqrt(Time)) + mu * v * sqrt(Time)
    d2 = d1 - (mu + 0.5) * v * sqrt(Time)
    d3 = log(U ^ 2 / (S * X)) / (v * sqrt(Time)) + (mu - 1) * v * sqrt(Time)
    d4 = d3 - (mu - 0.5) * v * sqrt(Time)
    e1 = log(L ^ 2 / (S * X)) / (v * sqrt(Time)) + mu * v * sqrt(Time)
    e2 = e1 - (mu + 0.5) * v * sqrt(Time)
    e3 = log(L ^ 2 / (S * X)) / (v * sqrt(Time)) + (mu - 1) * v * sqrt(Time)
    e4 = e3 - (mu - 0.5) * v * sqrt(Time)

    # Compute Value:

    Value = (eta * 1 / (U - L) * (S * exp((b - r) * Time) * S ^ (-2 *
             mu) * (S * X) ^ (mu + 0.5) / (2 * (mu + 0.5)) * ((U ^ 2 /
             (S * X)) ^ (mu + 0.5) * CND(eta * d1) - lambda1 * CND(eta
             * d2) - (L ^ 2 / (S * X)) ^ (mu + 0.5) * CND(eta * e1) +
             lambda1 * CND(eta * e2)) - X * exp(-r * Time) * S ^ (-2 *
             (mu - 1)) * (S * X) ^ (mu - 0.5) / (2 * (mu - 0.5)) * ((U
             ^ 2 / (S * X)) ^ (mu - 0.5) * CND(eta * d3) - lambda2 *
             CND(eta * d4) - (L ^ 2 / (S * X)) ^ (mu - 0.5) * CND(eta
             * e3) + lambda2 * CND(eta * e4))))
    ### print(Value)

    # Continue:
    if (TypeFlag == "cdi" || TypeFlag == "pui") {
        SoftBarrier = Value }
    if (TypeFlag == "cdo") {
        SoftBarrier = GBSOption("c", S, X, Time, r, b, v)@price - Value }
    if (TypeFlag == "puo") {
        SoftBarrier = GBSOption("p", S, X, Time, r, b, v)@price - Value }

    # Parameters:
    # TypeFlag = c("cdi", "cdo", "pdi", "pdo"), S, X, L, U, Time ,
    #   r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$L = L
    param$U = U
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Soft Barrier Option"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = SoftBarrier,
        title = title,
        description = description
        )
}


################################################################################

