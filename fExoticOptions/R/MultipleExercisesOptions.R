
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
# FUNCTION:                       MULTIPLE EXERCISES OPTIONS:
#  ExecutiveStockOption            Executive Stock Option
#  ForwardStartOption              Forward Start Option
#  RatchetOption                   Ratchet [Compound] Option
#  TimeSwitchOption                Time Switch Option
#  SimpleChooserOption             Simple Chooser Option
#  ComplexChooserOption            Complex Chooser Option
#  OptionOnOption                  Options On Options
#  HolderExtendibleOption          Holder Extendible Option
#  WriterExtendibleOption          Writer Extendible Option
################################################################################


ExecutiveStockOption =
    function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, lambda,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Executive stock options

    # References:
    #   Jennergren and Naslund (1993)
    #   Haug, Chapter 2.1

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Calculate Price:
    result = (exp (-lambda * Time) *
              GBSOption(TypeFlag = TypeFlag, S = S, X = X,
                        Time = Time, r = r, b = b, sigma = sigma)@price)

    # Parameters:
    # TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, lambda
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$lambda = lambda

    # Add title and description:
    if (is.null(title)) title = "Executive Stock Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


ForwardStartOption =
    function(TypeFlag = c("c", "p"), S, alpha, time1, Time2, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Forward Start Options

    # References:
    #   Rubinstein (1990)
    #   Haug, Chapter 2.2

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Compute Settings:
    Time = time1
    time = Time2

    # Compute Price:
    result = (S * exp ((b - r) * time ) *
              GBSOption(TypeFlag, S = 1, X = alpha, Time = Time-time,
                        r = r, b = b, sigma = sigma)@price)

    # Parameters:
    # TypeFlag = c("c", "p"), S, alpha, time1, Time2, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$alpha = alpha
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Forward Start Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


RatchetOption =
    function(TypeFlag = c("c", "p"), S, alpha, time1, Time2, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Ratchet Option,
    #   other names are MovingStrikeOption or CliquetOption

    # References:
    #   Haug, Chapter 2.3

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Calculate Price
    result = 0
    for ( i in 1:length(Time2) ) {
        result = (result +
                  ForwardStartOption(TypeFlag = TypeFlag, S = S, alpha = alpha,
                                     time1 = time1[i], Time2 = Time2[i],
                                     r = r, b = b, sigma = sigma)@price) }

    # Parameters:
    # TypeFlag = c("c", "p"), S, alpha, time1, Time2, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$alpha = alpha
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Ratchet Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


TimeSwitchOption =
    function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, A, m, dt,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Discrete time switch options

    # References:
    #   Pechtl (1995)
    #   Haug, Chapter 2.4

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Compute Settings:
    n = Time / dt
    Sum = 0
    if (TypeFlag == "c") Z = +1
    if (TypeFlag == "p") Z = -1

    # Calculate Price:
    Sum = 0
    for (I in (1:n)) {
        d = (log(S/X) + (b - sigma^2/2) * I * dt) / (sigma * sqrt(I * dt))
        Sum = Sum + CND (Z * d) * dt }
    result = A * exp (-r * Time) * Sum + dt * A * exp(-r * Time) * m

    # Parameters:
    # TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, A, m, dt
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$A = A
    param$m = m
    param$d = dt

    # Add title and description:
    if (is.null(title)) title = "Time Switch Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


SimpleChooserOption =
    function(S, X, time1, Time2, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simple Chooser Options

    # References:
    #   Rubinstein (1991)
    #   Haug, Chapter 2.5.1

    # FUNCTION:

    # Compute Settings:
    d = (log(S/X) + (b + sigma ^ 2 / 2) * Time2) / (sigma * sqrt(Time2))
    y = ((log(S/X) + b * Time2 + sigma ^ 2 * time1 / 2) /
         (sigma * sqrt(time1)))

    # Calculate Price:
    result = (S * exp ((b - r) * Time2) * CND(d) - X * exp(-r * Time2)
              * CND(d - sigma * sqrt(Time2)) - S * exp ((b - r) *
              Time2) * CND(-y) + X * exp(-r * Time2) * CND(-y + sigma
              * sqrt(time1)))

    # Parameters:
    # S, X, time1, Time2, r, b, sigma
    param = list()
    param$S = S
    param$X = X
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Simple Chooser Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


ComplexChooserOption =
    function(S, Xc, Xp, Time, Timec, Timep, r, b, sigma, doprint = FALSE,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Complex Chooser Options

    # References:
    #   Haug, Chapter 2.5.2

    # FUNCTION:

    # Compute Settings:
    Tc = Timec
    Tp = Timep

    # Calculate Price:
    CriticalValueChooser =
        function(S, Xc, Xp, Time, Tc, Tp, r, b, sigma){
            Sv = S
            ci = GBSOption("c", Sv, Xc, Tc - Time, r, b, sigma)@price
            Pi = GBSOption("p", Sv, Xp, Tp - Time, r, b, sigma)@price
            dc = GBSGreeks("Delta", "c", Sv, Xc, Tc - Time, r, b, sigma)
            dp = GBSGreeks("Delta", "p", Sv, Xp, Tp - Time, r, b, sigma)
            yi = ci - Pi
            di = dc - dp
            epsilon = 0.001
            # Newton-Raphson:
            while (abs(yi) > epsilon) {
                Sv = Sv - (yi) / di
                ci = GBSOption("c", Sv, Xc, Tc - Time, r, b, sigma)@price
                Pi = GBSOption("p", Sv, Xp, Tp - Time, r, b, sigma)@price
                dc = GBSGreeks("Delta", "c", Sv, Xc, Tc - Time, r, b, sigma)
                dp = GBSGreeks("Delta", "p", Sv, Xp, Tp - Time, r, b, sigma)
                yi = ci - Pi
                di = dc - dp }
            result = Sv
            result}

    # Complex chooser options:
    I = CriticalValueChooser (S, Xc, Xp, Time, Tc, Tp, r, b, sigma)
    if (doprint) {
        cat("\nCritical Value:\n")
        print(I)
        cat ("\n")}
    d1 = (log(S / I) + (b + sigma ^ 2 / 2) * Time) / (sigma * sqrt(Time))
    d2 = d1 - sigma * sqrt (Time)
    y1 = (log(S / Xc) + (b + sigma ^ 2 / 2) * Tc) / (sigma * sqrt(Tc))
    y2 = (log(S / Xp) + (b + sigma ^ 2 / 2) * Tp) / (sigma * sqrt(Tp))
    rho1 = sqrt (Time / Tc)
    rho2 = sqrt (Time / Tp)

    result = (S * exp ((b - r) * Tc) * CBND(d1, y1, rho1) -
              Xc * exp(-r * Tc) * CBND(d2, y1 - sigma * sqrt(Tc), rho1) -
              S * exp((b - r) * Tp) * CBND(-d1, -y2, rho2) +
              Xp * exp(-r * Tp) * CBND(-d2, -y2 + sigma * sqrt(Tp), rho2))

    # Parameters:
    # S, Xc, Xp, Time, Timec, Timep, r, b, sigma
    param = list()
    param$S = S
    param$Xc = Xc
    param$Xp = Xp
    param$time = time
    param$Timec = Timec
    param$Timep = Timep
    param$r = r
    param$b = b
    param$sigma = sigma
    param$criticalValue = I

    # Add title and description:
    if (is.null(title)) title = "Complex Chooser Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


OptionOnOption =
    function(TypeFlag = c("cc", "cp", "pc", "pp"), S, X1, X2, time1, Time2, r,
             b, sigma, doprint = FALSE, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Option on Option

    # References:
    #   Geske (1977), Geske (1979b), Hodges and Selby (1987),
    #   Rubinstein (1991a) et al.
    #   Haug, Chpater 2.6

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Compute Settings:
    Time = time1
    time = Time2

    # Internal Function:
    CriticalValueOptionOnOption =
        function(TypeFlag, X1, X2, Time, r, b, sigma) {
            # Calculation of critical price options on options
            Si = X1
            ci = GBSOption(TypeFlag, Si, X1, Time, r, b, sigma)@price
            di = GBSGreeks("Delta", TypeFlag, Si, X1, Time, r, b, sigma)
            epsilon = 0.000001
            # Newton-Raphson algorithm:
            while (abs(ci - X2) > epsilon) {
                Si = Si - (ci - X2) / di
                ci = GBSOption(TypeFlag, Si, X1, Time, r, b, sigma)@price
                di = GBSGreeks("Delta", TypeFlag, Si, X1, Time, r, b, sigma) }
            result = Si
            result }

    # Option On Option:
    T2 = Time
    t1 = time
    TypeFlag2 = "p"
    if (TypeFlag == "cc" || TypeFlag == "pc") TypeFlag2 = "c"
    I = CriticalValueOptionOnOption(TypeFlag2, X1, X2, T2-t1, r, b, sigma)
    if (doprint) { cat("\nCriticalValue: ", I, "\n") }
    rho = sqrt (t1 / T2)
    y1 = (log(S / I) + (b + sigma ^ 2 / 2) * t1) / (sigma * sqrt(t1))
    y2 = y1 - sigma * sqrt (t1)
    z1 = (log(S / X1) + (b + sigma ^ 2 / 2) * T2) / (sigma * sqrt(T2))
    z2 = z1 - sigma * sqrt (T2)
    if (TypeFlag == "cc")
        result = (S * exp ((b - r) * T2) * CBND(z1, y1, rho) -
                  X1 * exp(-r * T2) * CBND(z2, y2, rho) - X2 * exp(-r * t1) *
                  CND(y2))
    if (TypeFlag == "pc")
        result = (X1 * exp (-r * T2) * CBND(z2, -y2, -rho) -
                  S * exp((b - r) * T2) * CBND(z1, -y1, -rho) + X2 *
                  exp(-r * t1) * CND(-y2))
    if (TypeFlag == "cp")
        result = (X1 * exp (-r * T2) * CBND(-z2, -y2, rho) -
                  S * exp((b - r) * T2) * CBND(-z1, -y1, rho) - X2 *
                  exp(-r * t1) * CND(-y2))
    if (TypeFlag == "pp")
        result = (S * exp ((b - r) * T2) * CBND(-z1, y1, -rho) -
                  X1 * exp(-r * T2) * CBND(-z2, y2, -rho) + exp(-r * t1) *
                  X2 * CND(y2))

    # Parameters:
    # TypeFlag = c("cc", "cp", "pc", "pp"), S, X1, X2, time1, Time2, r, b, sigma
    param = list()
    param$S = S
    param$X1 = X1
    param$X2 = X2
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma
    param$criticalValue = I

    # Add title and description:
    if (is.null(title)) title = "Option On Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


HolderExtendibleOption =
    function(TypeFlag = c("c", "p"), S, X1, X2, time1, Time2, r, b, sigma, A,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Options that can be extended by the Holder

    # References:
    #   Haug, Chapter 2.7.1

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Calculate Price:
    HolderExtendible = NA
    if (TypeFlag == "c") {
        result = (max(c(S-X1, GBSOption(TypeFlag = "c", S = S, X = X2,
                                        Time = Time2-time1, r = r, b = b,
                                        sigma = sigma)@price - A, 0))) }
    if (TypeFlag == "p") {
        result = (max(c(X1-S, GBSOption(TypeFlag = "p", S = S, X = X2,
                                        Time = Time2-time1, r = r, b = b,
                                        sigma = sigma)@price - A, 0))) }

    # Parameters:
    # TypeFlag = c("c", "p"), S, X1, X2, time1, Time2, r, b, sigma, A
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X1 = X1
    param$X2 = X2
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma
    param$A = A

    # Add title and description:
    if (is.null(title)) title = "Holder Extendible Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


# ------------------------------------------------------------------------------


WriterExtendibleOption =
    function(TypeFlag = c("c", "p"), S, X1, X2, time1, Time2, r, b, sigma,
             title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Writer Extendible Options

    # References:
    #   Haug, Chapter 2.7.2

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]
    rho = sqrt (time1 / Time2)
    z1 = (log(S/X2) + (b + sigma^2 / 2) * Time2) / (sigma * sqrt(Time2))
    z2 = (log(S/X1) + (b + sigma^2 / 2) * time1) / (sigma * sqrt(time1))

    # Calculate Price:
    if (TypeFlag == "c")
        result = (GBSOption(TypeFlag, S, X1, time1, r, b, sigma)@price
                 + S * exp((b - r) * Time2) * CBND(z1, -z2, -rho) - X2
                 * exp(-r * Time2) * CBND(z1 - sqrt(sigma^2 * Time2),
                 -z2 + sqrt(sigma^2 * time1), -rho))
    if (TypeFlag == "p")
        result = (GBSOption(TypeFlag, S, X1, time1, r, b, sigma)@price
                  + X2 * exp(-r * Time2) * CBND(-z1 + sqrt(sigma^2 *
                  Time2), z2 - sqrt(sigma^2 * time1), -rho) - S *
                  exp((b - r) * Time2) * CBND(-z1, z2, -rho))

    # Parameters:
    # TypeFlag = c("c", "p"), S, X1, X2, time1, Time2, r, b, sigma
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X1 = X1
    param$X2 = X2
    param$time1 = time1
    param$Time2 = Time2
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Writer Extendible Option Valuation"
    if (is.null(description)) description = as.character(date())

    # Return Value:
    new("fOPTION",
        call = match.call(),
        parameters = param,
        price = result,
        title = title,
        description = description
        )
}


################################################################################

