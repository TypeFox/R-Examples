
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
# MOMENT MATCHING:                   DESCRIPTION:
#  MomentMatchedAsianOption           Valuate moment matched option prices
#  .LevyTurnbullWakemanAsianOption     Log-Normal Approximation
#  .MilevskyPosnerAsianOption          Reciprocal-Gamma Approximation
#  .PosnerMilevskyAsianOption          Johnson Type I Approximation
#  MomentMatchedAsianDensity          Valuate moment matched option densities
#  .LevyTurnbullWakemanAsianDensity    Log-Normal Approximation
#  .MilevskyPosnerAsianDensity         Reciprocal-Gamma Approximation
#  .PosnerMilevskyAsianDensity         Johnson Type I Approximation
# GRAM CHARLIER SERIES EXPANSION:    DESCRIPTION:
#  GramCharlierAsianOption            Calculate Gram-Charlier option prices
#  .GramCharlierAsianDensity          NA
# STATE SPACE MOMENTS:               DESCRIPTION:
#  AsianOptionMoments                 Methods to calculate Asian Moments
#  .DufresneAsianOptionMoments         Moments from Dufresne's Formula
#  .AbrahamsonAsianOptionMoments       Moments from Abrahamson's Formula
#  .TurnbullWakemanAsianOptionMoments  First 2 Moments from Turnbull-Wakeman
#  .TolmatzAsianOptionMoments          Asymptotic Behavior after Tolmatz
# STATE SPACE DENSITIES:              DESCRIPTION:
#  StateSpaceAsianDensity              NA
#  .Schroeder1AsianDensity             NA
#  .Schroeder2AsianDensity             NA
#  .Yor1AsianDensity                   NA
#  .Yor2AsianDensity                   NA
#  .TolmatzAsianDensity                NA
#  .TolmatzAsianProbability            NA
# PARTIAL DIFFERENTIAL EQUATIONS:     DESCRIPTION:
#  PDEAsianOption                      PDE Asian Option Pricing
#   .ZhangAsianOption                   Asian option price by Zhang's 1D PDE
#    ZhangApproximateAsianOption
#   .VecerAsianOption                   Asian option price by Vecer's 1D PDE
# LAPLACE INVERSION:                  DESCRIPTION:
#   GemanYorAsianOption                Asian option price by Laplace Inversion
#   gGemanYor                          Function to be Laplace inverted
# SPECTRAL EXPANSION:                 DESCRIPTION:
#   LinetzkyAsianOption                Asian option price by Spectral Expansion
#   gLinetzky                          Function to be integrated
# BOUNDS ON OPTION PRICES:            DESCRIPTION:
#   BoundsOnAsianOption                 Lower and upper bonds on Asian calls
#    CurranThompsonAsianOption          From Thompson's continuous limit
#    RogerShiThompsonAsianOption        From Thompson's single integral formula
#    ThompsonAsianOption                Thompson's upper bound
# SYMMETRY RELATIONS:                 DESCRIPTION:
#   CallPutParityAsianOption           Call-Put parity Relation
#   WithDividendsAsianOption           Adds dividends to Asian Option Formula
# TABULATED RESULTS:                  DESCRIPTION:
#   FuMadanWangTable                   Table from Fu, Madan and Wang's paper
#   FusaiTaglianiTable                 Table from Fusai und tagliani's paper
#   GemanTable                         Table from Geman's paper
#   LinetzkyTable                      Table from Linetzky's paper
#   ZhangTable                         Table from Zhang's paper
#   ZhangLongTable                     Long Table from Zhang's paper
#   ZhangShortTable                    Short Table from Zhang's paper
################################################################################


################################################################################
# MOMENT MATCHING:


MomentMatchedAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, table = NA, method = c("LN", "RG", "JI"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates price for Asian options based on moment matching
    #   LN:  Levy-Turnbull-Wakeman Log-Normal Approximation
    #   RG:  Milevsky-Posner Reciprocal-Gamma Approximation
    #   JI:  Posner-Milevski Johnson Type I Approximation

    # FUNCTION:

    # Set Default TypeFlag and Method, if no other is selected:
    TypeFlag = TypeFlag[1]
    method = method[1]

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[,1]
        X = table[,2]
        Time = table[,3]
        r = table[,4]
        sigma = table[,5]
    }
    call = rep(0, length=length(S))

    # Log-Normal Approximation:
    if (method == "LN") {
        for ( i in 1:length(S) ) {
            moments = masian(Time = Time[i], r = r[i],
                sigma = sigma[i])$rawMoments
            moments = moments / Time[i]^(1:4)
            meanlog = ( 2*log(moments[1]) - log(moments[2])/2 )
            sdlog = ( sqrt ( log(moments[2]) - 2*log(moments[1]) ) )
            d2 = ( -log(X[i]/S[i]) + meanlog*Time[i] ) / ( sdlog*sqrt(Time[i]) )
            d1 = d2 + sdlog*sqrt(Time[i])
            call[i] = moments[1]*pnorm(d1)-(X[i]/S[i])*pnorm(d2)
        }
     }

     # Reciprocal-Gamma Approximation:
     if (method == "RG") {
        for ( i in 1:length(S) ) {
            moments = masian(Time = Time[i], r = r[i],
                sigma = sigma[i])$rawMoments
            moments = moments / Time[i]^(1:4)
            alpha = (2*moments[2] - moments[1]^2) / (moments[2] - moments[1]^2)
            beta = (moments[2] - moments[1]^2) / (moments[1]*moments[2])
            call[i] = moments[1]*pgamma(S[i]/X[i], shape=alpha-1, scale=beta) -
                (X[i]/S[i])*pgamma(S[i]/X[i], shape=alpha, scale=beta)
        }
     }

     # Johnson-Type-I Approximation:
     if (method == "JI") {
        for ( i in 1:length(S) ) {
            cmoments = masian(Time = Time[i], r = r[i],
                sigma = sigma[i])$centralMoments
            cmoments = cmoments / Time[i]^(1:4)
            mu = cmoments[1]
            varsigma = sqrt(cmoments[2])
            eta = cmoments[3] / varsigma^3
            omega = 0.5 * ( 8 + 4*eta^2 + 4*sqrt(4*eta^2+eta^4) )^(1/3)
            omega = omega + 1/omega - 1
            b = 1 / sqrt(log(omega))
            a = 0.5 * b * log(omega*(omega-1)/varsigma^2)
            d = sign(eta)
            c = d*mu - exp(  (0.5/b-a)/b  )
            Q = a + b*log((X[i]/S[i]-c)/d)
            call[i] = mu - X[i]/S[i] + (X[i]/S[i] - c) * pnorm( Q ) -
                d * exp ( (1-2*a*b)/(2*b^2) ) * pnorm ( Q - 1/b )
        }
    }

    # Call Price:
    Price = S* exp(-r*Time) * call

    # Put Price:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        Price = Price - Parity
    }

    # Return Value:
    option = list(
        price = Price,
        call = match.call() )
    class(option) = "option"
    option
}


# ------------------------------------------------------------------------------


.LevyTurnbullWakemanAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    MomentMatchedAsianOption(TypeFlag[1], S, X, Time, r, sigma, method = "LN")
}


# ------------------------------------------------------------------------------


.MilevskyPosnerAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    MomentMatchedAsianOption(TypeFlag[1], S, X, Time, r, sigma, method = "RG")
}


# ------------------------------------------------------------------------------


.PosnerMilevskyAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    MomentMatchedAsianOption(TypeFlag[1], S, X, Time, r, sigma, method = "JI")
}


# ------------------------------------------------------------------------------


MomentMatchedAsianDensity =
function(x, Time = 1, r = 0.09, sigma = 0.30, method = c("LN", "RG", "JI"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates price for Asian options based on moment matching
    #   LN:  Levy-Turnbull-Wakeman Log-Normal Approximation
    #   RG:  Milevsky-Posner Reciprocal-Gamma Approximation
    #   JI:  Posner-Milevski Johnson Type I Approximation

    # FUNCTION:

    # Set Default Method, if no other is selected:
    method = method[1]

    # Log-Normal Approximation:
    if (method == "LN") {
        moments = masian(Time = Time, r = r, sigma = sigma)$rawMoments
        moments = moments / Time^(1:4)
        meanlog = 2*log(moments[1]) - log(moments[2])/2
        sdlog = sqrt ( log(moments[2]) - 2*log(moments[1]) )
        density = dlnorm(x = x, meanlog = meanlog, sdlog = sdlog)
     }

     # Reciprocal-Gamma Approximation:
     if (method == "RG") {
        moments = masian(Time = Time, r = r, sigma = sigma)$rawMoments
        moments = moments / Time^(1:4)
        alpha = (2*moments[2] - moments[1]^2) / (moments[2] - moments[1]^2)
        beta = (moments[2] - moments[1]^2) / (moments[1]*moments[2])
        density = drgam(x = x, alpha = alpha, beta = beta)
     }

     # Johnson Type I Approximation:
     if (method == "JI") {
        cmoments = masian(Time = Time, r = r, sigma = sigma)$centralMoments
        cmoments = cmoments / Time^(1:4)
        mu = cmoments[1]
        varsigma = sqrt(cmoments[2])
        eta = cmoments[3] / varsigma^3
        omega = 0.5 * ( 8 + 4*eta^2 + 4*sqrt(4*eta^2+eta^4) )^(1/3)
        omega = omega + 1/omega - 1
        b = 1 / sqrt(log(omega))
        a = 0.5 * b * log(omega*(omega-1)/varsigma^2)
        d = sign(eta)
        c = d*mu - exp(  (0.5/b-a)/b  )
        density = djohnson(x = x, a = a, b = b, c = c, d = d)
    }

    # Return Value:
    density
}


# ------------------------------------------------------------------------------


.LevyTurnbullWakemanAsianDensity =
function(x, Time = 1, r = 0.09, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    MomentMatchedAsianDensity(x, Time, r, sigma, method = "LN")
}


# ------------------------------------------------------------------------------


.MilevskyPosnerAsianDensity =
function(x, Time = 1, r = 0.09, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    MomentMatchedAsianDensity(x, Time, r, sigma, method = "RG")
}


# ------------------------------------------------------------------------------


.PosnerMilevskyAsianDensity =
function(x, Time = 1, r = 0.09, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Return Value:
    MomentMatchedAsianDensity(x, Time, r, sigma, method = "JI")
}


################################################################################
# GRAM CHARLIER:


GramCharlierAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, table = NA, method = c("LN", "RG", "JI"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate arithmetic Asian options price using
    #   Gram Charlier Statistical Series Expansion around
    #   LN:  Log-Normal Distribution
    #   RG:  Reciprocal-Gamma Distribution
    #   JI:  Johnson-Type-I Distribution

    # FUNCTION:

    # Select Method:
    TypeFlag = TypeFlag[1]
    method = method[1]

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[,1]
        X = table[,2]
        Time = table[,3]
        r = table[,4]
        sigma = table[,5]
    }

    # Calculate Price:
    Price = MomentMatchedAsianOption("c", S = S, X = X, Time = Time, r = r,
        sigma = sigma, method=method)$price
    gc3 = gc4 = rep(0, length(Price))

    # Log-Normal Approximation:
    if (method == "LN") {
        for ( i in 1:length(S) ) {
            moments = masian(Time[i], r[i], sigma[i])$rawMoments/Time[i]^(1:4)
            meanlog = 2*log(moments[1]) - log(moments[2])/2
            sdlog = sqrt ( log(moments[2]) - 2*log(moments[1]) )
            asian.cm = masian(Time[i], r[i],
                sigma[i])$centralMoments/Time[i]^(1:4)
            lnorm.cm = mlognorm(meanlog, sdlog)$centralMoments/Time[i]^(1:4)
            kappa = (asian.cm-lnorm.cm)
            gc3[i] = kappa[3]*dlognorm(X[i]/S[i], meanlog, sdlog, deriv=1)/6
            gc4[i] = kappa[4]*dlognorm(X[i]/S[i], meanlog, sdlog, deriv=2)/24
        }
    }

    # Reciprocal-Gamma Approximation:
    if (method == "RG" ) {
        for ( i in 1:length(S) ) {
            moments = masian(Time[i], r[i], sigma[i])$rawMoments/Time[i]^(1:4)
            alpha = (2*moments[2] - moments[1]^2) / (moments[2] - moments[1]^2)
            beta = (moments[2] - moments[1]^2) / (moments[1]*moments[2])
            asian.cm = masian(Time[i], r[i],
                sigma[i])$centralMoments/Time[i]^(1:4)
            rgam.cm = mrgam(alpha, beta)$centralMoments/Time[i]^(1:4)
            kappa = (asian.cm-rgam.cm)
            gc3[i] = kappa[3]*drgam(X[i]/S[i], alpha, beta, deriv = 1)/6
            gc4[i] = kappa[4]*drgam(X[i]/S[i], alpha, beta, deriv = 2)/24
        }
    }

    # Johnson-Type-I Approximation:
    if (method == "JI" ) {
        for ( i in 1:length(S) ) {
            cmoments = masian(Time = Time[i], r = r[i],
                sigma = sigma[i])$centralMoments
            cmoments = cmoments / Time[i]^(1:4)
            mu = cmoments[1]
            varsigma = sqrt(cmoments[2])
            eta = cmoments[3] / varsigma^3
            omega = 0.5 * ( 8 + 4*eta^2 + 4*sqrt(4*eta^2+eta^4) )^(1/3)
            omega = omega + 1/omega - 1
            b = 1 / sqrt(log(omega))
            a = 0.5 * b * log(omega*(omega-1)/varsigma^2)
            d = sign(eta)
            c = d*mu - exp(  (0.5/b-a)/b  )
            asian.cm = cmoments
            johnson.cm = cmoments
            skewness = sqrt((omega-1)*(omega+2)^2)
            kurtosis = omega^4 + 2*omega^3 + 3* omega^2 - 3
            johnson.cm[3] = skewness * varsigma^3
            johnson.cm[4] = kurtosis * varsigma^4
            kappa = (asian.cm-johnson.cm)
            gc3[i] = kappa[3]*djohnson(X[i]/S[i], a, b, c, d, deriv = 1)/6
            gc4[i] = kappa[4]*djohnson(X[i]/S[i], a, b, c, d, deriv = 2)/24
        }
    }

    # Gram-Charlier Approximated Call Price:
    Price = Price - S * exp(-r*Time) * (gc3-gc4)

    # Put Price:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        Price = Price - Parity
    }

    # Return Value:
    option = list(
        price = Price,
        call = match.call() )
    class(option) = "option"
    option
}


.GramCharlierAsianDensity =
function(Time = 1, r = 0.09, sigma = 0.30, method = c("LN", "RG", "JI"))
{   # A function ported by Diethelm Wuertz

    # Return Value:
    NA
}


################################################################################
# STATE SPACE MOMENTS:


AsianOptionMoments =
function(M = 4, Time = 1, r = 0.045, sigma = 0.30, log = FALSE,
method = c("A", "D", "TW", "T"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Asian Moments using several approaches:
    #    A -  Moments from Abrahamson's Formula
    #    D -  Moments from Dufresne's Formula
    #    TW - First 2 Moments from Turnbull-Wakeman
    #    T - Asymptotic Behavior after Tolmatz

    # FUNCTION:

    # Settings:
    method = method[1]
    result = NA

    # Abrahamson Formula:
    if (method == "A") result =
        .AbrahamsonAsianOptionMoments(M = M, Time = Time, r = r,
            sigma = sigma)

    # Dufresne Formula:
    if (method == "D") result =
        .DufresneAsianOptionMoments(M = M, Time = Time, r = r,
            sigma = sigma)

    # Tolmatz Formula - Asymptotic Behavior:
    if (method == "T") result =
        .TolmatzAsianOptionMoments(M = M, Time = Time, r = r,
            sigma = sigma, log = log)

    # Turnbull Wakeman - Explicit 1st and Second Moment:
    if (method == "TW") result =
        .TurnbullWakemanAsianOptionMoments(M = M, Time = Time, r = r,
            sigma = sigma)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


.DufresneAsianOptionMoments =
function(M = 4, Time = 1, r = 0.045, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Moments of Asian Options Density
    #   according to the formula of Dufresne-GemanYor

    # FUNCTION:

    # Calculates: E[(A_\tau^{(\nu)})^n]
    moments = function (M, tau, nu) {
        d = function(j, n, beta) {
            d = 2^n
            for (i in 0:n) if (i != j) d = d  / ( (beta+j)^2 - (beta+i)^2 )
            d
        }
        moments = rep(0, length=M)
        for (n in 1:M) {
            moments[n] = 0
            for (j in 0:n) moments[n] = moments[n] +
                d(j, n, nu/2)*exp(2*(j^2+j*nu)*tau)
            moments[n] = prod(1:n) * moments[n] / (2^(2*n)) }
        moments
    }

    # Calculate for:
    tau = sigma^2*Time/4
    nu = 2*r/sigma^2-1

    # Return Value:
    (4/sigma^2)^(1:M) * moments(M, tau, nu)
}


# ------------------------------------------------------------------------------


.AbrahamsonAsianOptionMoments =
function (M = 4, Time = 1, r = 0.045, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Moments of Asian Options Density
    #   according to the formula of Abrahamson

    # FUNCTION:

    # Calculates: E[(A_\tau^{(\nu)})^n]
    moments = function (M, tau, nu) {
        moments = rep(0, times = M)
        for (N in 1:M) {
            d = c =  2 * ( (1:N)^2 + (1:N)*nu )
            for (m in 1:N) {
                for (j in 1:N) {
                    if (j!= m) d[m] = d[m]*(c[m]-c[j])
                }
                d[m] = exp(c[m]*tau) / d[m]
            }
            moments[N] = prod(1:N) * ( sum(d) + (-1)^N/prod(c) )
        }
        moments
    }

    # Calculate for:
    tau = sigma^2*Time/4
    nu = 2*r/sigma^2-1

    # Return Value:
    (4/sigma^2)^(1:M) * moments(M, tau, nu)
}


# ------------------------------------------------------------------------------


.TurnbullWakemanAsianOptionMoments =
function (M = 2, Time = 1, r = 0.045, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the first two moments as derived explicitly
    #   by Turnbull and Wakeman. It can serve as a test for
    #   other implementations.

    # Note:
    #   Maximum M is 2!

    # FUNCTION:

    # Moments:
    moments = rep(NA, times = M)
    if (M == 1 || M == 2)
        moments[1] = (exp(r*Time)-1)/(r*Time)
    if (M == 2)
        moments[2] =
            2*exp((2*r+sigma^2)*Time)/ ((r+sigma^2)*(2*r+sigma^2)*Time^2) +
            (2/(r*Time^2)) * ( 1/(2*r+sigma^2) - exp(r*Time)/(r+sigma^2) )

    # Return Value:
    moments
}


# ------------------------------------------------------------------------------


.TolmatzAsianOptionMoments =
function (M = 100, Time = 1, r = 0.045, sigma = 0.30, log = FALSE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Asymptotic Moments of Asian Options Density
    #   according to the formula of Tolmatz for nu=0 and Wuertz
    #   for nu different from zero - Log returns can be selected

    # FUNCTION:

    # Calculates: log { E[(A_\tau^{(\nu)})^n] }
    moments = function (M, tau, nu=0) {
        moments = rep(0, times=M)
        M = 1:M
        log.moments = -M*log(2) + lgamma(nu+M) - lgamma(nu+2*M) +
            2*(M^2+M*nu)*tau
        log.moments
    }

    # Calculate for:
    tau = sigma^2*Time/4
    nu = 2*r/sigma^2-1

    # Return Value:
    moments = (1:M)*log(4/sigma^2) + moments(M, tau, nu)

    # Return value:
    if (!log) moments = exp(moments)
    moments
}


################################################################################
# ASIAN DENSITY:


# STATE SPACE DENSITIES:              DESCRIPTION:
#  StateSpaceAsianDensity
#  .Schroeder1AsianDensity             S1
#  .Schroeder2AsianDensity             S2
#  .Yor1AsianDensity                   Y1
#  .Yor2AsianDensity                   Y2
#  .TolmatzAsianDensity                T
#  .TolmatzAsianProbability




################################################################################
# PDE SOLVER:


ZhangAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, table = NA, correction = TRUE, nint = 800, eps = 1.0e-8,
dt = 1.0e-10)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Valuates Asian options by Solving Zhang's one
    #   dimensional Partial Differential Equations

    # Source:
    #   For the Fortran Routine:
    #   TOMS ...

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[, 1]
        X = table[, 2]
        Time = table[, 3]
        r = table[, 4]
        sigma = table[, 5]
    }
    Price = rep(0, times = length(S))

    # Set Model Identifier:
    modsel = 2

    # Option Parameters:
    if (TypeFlag == "c") z = +1
    if (TypeFlag == "p") z = -1

    # PDE Parameters - Do not change:
    T0 = 0; Tout = 1
    np = 0; Price.by.S = 0
    mf = 12
    npde = 1; kord = 4; ncc = 2; maxder = 5

    # Fill Working Arrays:
    xbkpt = rep(0, times = nint+1)
    length.work = kord+npde*(4+9*npde)+(kord+(nint-1)*(kord-ncc)) *
      (3*kord+2+npde*(3*(kord-1)*npde+maxder+4))
    work = rep(0, times = length.work)
    length.iwork = (npde+1)*(kord+(nint-1)*(kord-ncc))
    iwork = rep(0, times = length.iwork)

    # Compute Prices:
    for ( i in 1:length(S) ) {
        result = .Fortran("asianval",
            as.double(z),
            as.double(S[i]),
            as.double(X[i]),
            as.double(X[i]),
            as.double(X[i]),
            as.double(Time[i]),
            as.double(r[i]),
            as.double(sigma[i]),
            as.double(T0),
            as.double(Tout),
            as.double(eps),
            as.double(dt),
            as.double(Price.by.S),
            as.integer(np),
            as.integer(modsel),
            as.integer(mf),
            as.integer(npde),
            as.integer(kord),
            as.integer(nint),
            as.integer(ncc),
            as.integer(maxder),
            as.double(X[i]/S[i]),
            as.double(xbkpt),
            as.double(work),
            as.integer(iwork),
            PACKAGE = "fAsianOptions"
            )
        Price[i] = result[[13]]*S[i]
    }

    # ?
    Price = Price +
        ZhangApproximateAsianOption(TypeFlag, S, X, Time, r, sigma, table)

    # Return Value:
    Price
}


ZhangApproximateAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, table = NA)
{
    # Settings:
    TypeFlag = TypeFlag[1]

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[, 1]
        X = table[, 2]
        Time = table[, 3]
        r = table[, 4]
        sigma = table[, 5]
    }

    # Compute:
    I = 0
    xi = (Time*X-I)*exp(-r*Time)/S - (1-exp(-r*Time))/r
    eta = (sigma^2/(4*r^3)) * (-3 + 2*r*Time + 4*exp(-r*Time) - exp(-2*r*Time))

    # Call:
    price = (S/Time) *
        ( -xi * pnorm(-xi/sqrt(2*eta)) + sqrt(eta/pi)*exp(-xi^2/(4*eta)) )
    if (TypeFlag == "c") {
        ans = price
    } else {
        ans = CallPutParityAsianOption(TypeFlag = "c", price,
            S, X, Time, sigma, r, table = table)
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


VecerAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, table = NA, nint = 800, eps = 1.0e-8, dt = 1.0e-10)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Valuates Asian options by Solving Vecer's one
    #   dimensional Partial Differential Equations

    # FUNCTION:

    # Vecer's PDE modeled by: modsel == 1
    # SUBROUTINE ASIANVAL(
    #   ZZ, SS, XS, XSMIN, XSMAX, TIME, RR, SIGMA,
    #   T0, TOUT, EPS, DT, PRICEBYS, NP, MODSEL,
    #   MF1, NPDE1, KORD1, MX1, NCC1, MAXDER1,
    #   XBYS, XBKPT, WORK, IWORK)

    # Settings:
    TypeFlag = TypeFlag[1]

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[, 1]
        X = table[, 2]
        Time = table[, 3]
        r = table[, 4]
        sigma = table[, 5]
    }
    Price = rep(0, times = length(S))

    # Set Model Identifier:
    modsel = 1

    # Option Parameters:
    if (TypeFlag == "c") z = +1
    if (TypeFlag == "p") z = -1

    # PDE Parameters - Do not change:
    T0 = 0
    Tout = 1
    np = 0
    Price.by.S = 0
    mf = 12
    npde = 1
    kord = 4
    ncc = 2
    maxder = 5

    # Fill Working Arrays:
    xbkpt = rep(0, times = nint+1)
    length.work = kord+npde*(4+9*npde)+(kord+(nint-1)*(kord-ncc)) *
      (3*kord+2+npde*(3*(kord-1)*npde+maxder+4))
    work = rep(0, times = length.work)
    length.iwork = (npde+1)*(kord+(nint-1)*(kord-ncc))
    iwork = rep(0, times = length.iwork)

    # Compute Prices:
    for ( i in 1:length(S) ) {
        result = .Fortran("asianval",
            as.double(z),
            as.double(S[i]),
            as.double(X[i]),
            as.double(X[i]),
            as.double(X[i]),
            as.double(Time[i]),
            as.double(r[i]),
            as.double(sigma[i]),
            as.double(T0),
            as.double(Tout),
            as.double(eps),
            as.double(dt),
            as.double(Price.by.S),
            as.integer(np),
            as.integer(modsel),
            as.integer(mf),
            as.integer(npde),
            as.integer(kord),
            as.integer(nint),
            as.integer(ncc),
            as.integer(maxder),
            as.double(X[i]/S[i]),
            as.double(xbkpt),
            as.double(work),
            as.integer(iwork),
            PACKAGE = "fAsianOptions"
            )
        Price[i] = result[[13]]*S[i]
    }

    # Return Value:
    Price
}


################################################################################
# LAPLACE INVERSION:


gGemanYor =
function(lambda, S = 100, X = 100, Time = 1, r = 0.05, sigma = 0.30,
log = FALSE, doplot = FALSE)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Calculates function to be Laplace inverted

    # Arguments:
    #   lambda - complex vector

    # Notes:
    #   Equation 4.9 with notation as in
    #       Sudler G.F. [1999], "Asian Options: Inverse Laplace
    #       Transform and Martingale Methods Revisited".

    # FUNCTION:

    # Settings:
    x = lambda
    g = rep(complex(real = 0, imaginary = 0), length = length(x))

    # Calculate for each lambda value from Kummer Function:
    # Note Kummer function is not vectorized in Indexes !
    for ( i in 1:length(x) ) {
        # Settings:
        nu = 2*r/(sigma^2) - 1
        mu = sqrt(2*lambda[i] + nu^2)
        q = (sigma^2)*X*Time/(4*S)
        gamma1 = (mu-nu)/2
        gamma2 = (mu+nu)/2
        # Convergence Parameters:
        a = gamma1-2 + 1
        b = gamma2+1 + a + 1
        z = -1/(2*q)
        # From Kummer Function [one of ...]:
        # Use logarithmic Kummer and Gamma Functions to prevent
        # from numerical overflow!
        g[i] = kummerM(-z, b-a, b, lnchf = 1) +
            a*log(-z) + z + cgamma(b-a, log = TRUE) - cgamma(b, log = TRUE) -
            log (lambda[i]*(lambda[i] - 2 - 2 * nu))
        if (!log) g[i] = exp(g[i])
    }

    # Plot function if desired:
    if (doplot) {
        if (!is.complex(lambda)) {
            lam = lambda
            xlab = "lambda"
            ylab = "g"
        } else {
            lam = Im(lambda)
            xlab = "Im(lambda)"
            ylab = "Re(g)"
        }
        lambda.min = 4*r/sigma^2
        cat("\nmin lambda:", lambda.min, "\n")
        # Function to be Laplace inverted:
        print(cbind(lambda, g))
        plot(lam, Re(g), type = "l", main = "Laplace Inverse",
            xlab = xlab, ylab = ylab)
        lines(lam, 0*Re(g))
        # Convergence Indexes:
        mu = sqrt(2*lambda + nu^2)
        gamma1 = (mu-nu)/2; a = Re(gamma1-2 + 1)
        gamma2 = (mu+nu)/2; b = Re(gamma2+1 + a + 1)
        plot(lam, b, ylim = c(min(c(a,b)), max(c(a,b))), type = "n",
            xlab = xlab, ylab = "a  b", main = "Convergence Indexes")
        lines(lam, a, col = "red")
        lines(lam, b, col = "blue")
        lines(lam, 0*b, type = "l", col = "black")
        lines(x = lambda.min*c(1,1), y = c( min(c(a,b)), max(c(a,b)) ) )
    }

    # Return Value:
    g
}


# ------------------------------------------------------------------------------


GemanYorAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, doprint = FALSE)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Valuate Asian options by Laplace inversion.

    # FUNCTION:

    # Parameters:
    TypeFlag = TypeFlag[1]
    nu = (2*r)/(sigma^2) - 1
    alpha = 1 / sigma^2 * (2*nu + 2)
    h = sigma^2*Time/4

    # Function to be inverted:
    fx = function(x, S, X, Time, roh, sigma) {
        nu = (2*roh)/(sigma^2) - 1
        alpha = 1 / sigma^2 * (2*nu + 2)
        h = sigma^2*Time/4
        zi = complex(real = 0, imaginary = x)
        zc = complex(real = alpha, imaginary = x)
        g2 = gGemanYor(lambda = zc, S = S, X = X,
            Time = Time, r = roh, sigma = sigma, log = TRUE)
        # Return Value:
        Re ( exp(zi*h + g2) / (2*pi) )
    }

    # Call:
    # Integrate stepwise until (hopefully) convergence is reached:
    q = (sigma^2)*X*Time/(4*S)
    delta = 10/(2*q)
    eps = 1.0e-20
    if (doprint) {
        cat("\nDelta:", delta)
        cat("\nS:", S, "X:", X)
        cat("\nTime:", Time, "r:", r, "sigma:", sigma, "\n")
    }
    i = 1
    I = integrate(fx, lower = 0, upper = delta,
        S = S, X = X, Time=Time, roh = r, sigma = sigma)
    Price = Increment = exp(-r*Time)*(S/h)*exp(alpha*h)*2*I$value
    while (abs(Increment)/abs(Price) > eps) {
        i = i+1
        I = integrate(fx, lower = (i-1)*delta, upper = i*delta,
            S = S, X = X, Time = Time, roh = r, sigma = sigma)
        Increment = exp(-r*Time)*(S/h)*exp(alpha*h)*2*I$value
        Price = Price + Increment
        if (doprint) print(c(i*delta, Price, Increment))
    }

    # Put:
    # Use Call-Put Parity:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        Price = Price - Parity
    }

    # Return Value:
    option = list(
        price = Price,
        call = match.call() )
    class(option) = "option"
    option
}


################################################################################
# SPECTRAL EXPANSION:


gLinetzky =
function(x, y, tau, nu, ip = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the function to be integrated for the Put Price
    #   in the expression for the spectral representation of
    #   $ P^{(\nu)} (k, \tau) $. Proposition 2 described bu eq. (16)

    # Note:
    #   Requires Confluent Hypergeometric Functions

    # Reference:
    #   [L] V. Linetzky, Spectral Expansions for Asian (Average Price)
    #   Options, Preprint, revised Version from October 2002

    # Function:
    result = V = rep(0, length = length(x))
    for (i in 1:length(x)) {
        p = x[i]
        if (p == 0) {
            result[i] = 0
        } else {
            z = 1/(2*y)
            kappa = -(nu+3)/2
            mu = complex(real = 0, imaginary = p/2)
            # V1 = exp( -(nu^2+p^2)*tau/2 ) * z^kappa * exp(-z/2)
              logV1 = -(nu^2+p^2)*tau/2 + kappa*log(z) - z/2
            # V2 = (abs(cgamma(nu/2+mu)))^2
            if (p < 100) {
                logV2 = log ( (abs(cgamma(nu/2+mu)))^2 )
            } else {
                # Shift by Pi - Take of the proper Phi
                g = cgamma(nu/2+mu, log = TRUE)
                r = abs(g)
                phi = atan(Im(g)/Re(g)) + pi
                logV2 = 2*r*cos(phi) }
            # V3 = sinh(pi*p) * p
            logV3 = pi*p + log(1/2 - exp(-2*pi*p)/2) + log(p)
            # Whittaker:
            if (p < 100) {
                V4 = Re ( whittakerW(z, kappa, mu, ip ) )
            } else {
                # Use: 2 * Re ( (cgamma(-2*mu)/cgamma(1/2-mu-kappa)) *
                # exp(-z/2) * z^(1/2+mu) * kummerM(z, 1/2+mu-kappa, 1+2*mu)
                g = log(2) + cgamma(-2*mu, log=TRUE) -
                    cgamma(1/2-mu-kappa, log=TRUE) - z/2 +(1/2+mu)*log(z) +
                    kummerM(z, 1/2+mu-kappa, 1+2*mu, lnchf=1, ip=ip)
                r = abs(g)
                # Shift by Pi - Take of the proper Phi
                phi = atan(Im(g)/Re(g)) + pi
                logV4 = r*cos(phi)
                argV4 = cos(r*sin(phi))
            }
            # Collect all terms:
            if ( p < 100) {
                result[i] = exp(logV1+logV2+logV3)*V4/(8*pi^2)
            } else {
                result[i] = exp(logV1+logV2+logV3+logV4)*argV4/(8*pi^2)
            }
            # print(c(p, logV5, logV5a, argV5, argV5a))
        }
    }

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


LinetzkyAsianOption =
function(TypeFlag = c("c", "p"), S = 2, X = 2, Time = 1, r = 0.02,
sigma = 0.1, table = NA, lower = 0, upper = 100, method = "adaptive",
subdivisions = 100, ip = 0, doprint = TRUE, doplot = TRUE,...)
{   # A function implemented by Diethelm Wuertz

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[,1]
        X = table[,2]
        Time = table[,3]
        r = table[,4]
        sigma = table[,5] }
    if (doprint) print(cbind(S, X, Time, r, sigma))
    Price = rep(0, length=length(S))

    # Settings:
    tau = sigma^2*Time/4
    k = tau*X/S
    nu = 2*r/sigma^2 - 1
    if (doprint) print(cbind(k, tau, nu))

    # Parameters:
    z = 1 / (2*k)
    kappa = -(nu+3)/2
    mu.max = upper/2
    if (doprint) print(cbind(z, kappa, mu.max))

    # Calculate Spectral Measure:
    P = P1 = rep(0, times=length(S))
    for (i in 1:length(S)) {
        if (method == "adaptive") {
            P[i] = integrate(gLinetzky, lower = lower, upper = upper,
                y = k[i], tau = tau[i], nu = nu[i], ip = ip,
                subdivisions = subdivisions,
                rel.tol = .Machine$double.eps^0.25,
                abs.tol = .Machine$double.eps^0.25,
                stop.on.error = FALSE)$value
        }
        if (method == "trapez") {
            x = seq(lower, upper, length = subdivisions+1)
            delta = (upper-lower)/subdivisions
            F = gLinetzky(x = x, y = k[i], tau[i], nu[i])
            # print(c(F[1], F[length(F)], min(F), max(F)))
            P[i] = ( sum(F)-(F[1]+F[length(F)])/2 ) * delta
        }
        if (method == "simpson") {
            x = seq(lower, upper, length = subdivisions+1)
            delta = (upper-lower)/subdivisions
            F = gLinetzky(x = x, y = k[i], tau = tau[i], nu = nu[i], ip = ip)
            # print(c(F[1], F[length(F)], min(F), max(F)))
            FF = matrix(F[2:length(F)], byrow = TRUE, ncol = 2)
            P[i] = (F[1]+4*sum(FF[,1])+2*sum(FF[,2])-F[length(F)]) *
                delta[i]/3 }
        # For nu < 0 add:
        P1[i] = 0
        if (nu[i] < 0) {
            z = 1/(2*k[i])
            P1[i] = (2*k[i]*pgamma(abs(nu[i]), z) - pgamma(abs(nu[i])-1, z)) /
                ( 2 * gamma(abs(nu[i])) )
            P[i] = P[i] + P1[i]
        }
        # Plot:
        if (doplot) {
            x = seq(lower, upper, length = subdivisions)
            F = gLinetzky(x = x, y = k[i], tau = tau[i], nu = nu[i])
            plot(x, F, type = "l")
            lines(x, 0*x, col = "red")
            lines(x, F)
        }
    }
    if (doprint) print(cbind(P, P1))

    # Derive Call/Put Price:

    # Put Price:
    Linetzky = exp(-r*Time) * (S/tau) * P

    # Call Price:
    if (TypeFlag == "c") {
        Linetzky = Linetzky + S*(1-exp(-r*Time))/(r*Time) - X*exp(-r*Time) }

    # Return Value:
    option = list(
        price = Linetzky,
        call = match.call() )
    class(option) = "option"
    option
}



################################################################################
# ASIAN BOUNDS:

# Note:
#   We have not implemented the formula for the upper bound derived
#   by Rogers and Shi. Thompson's upper bound formula is much more
#   precise and therefore we have concentrated ourself on their
#   approach.


BoundsOnAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30, table = NA, method = c("CT", "RST", "T"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Bounds on Asian Option Prices
    #   CT  - Curran-Thompson Lower Bound
    #   RST - Roger-Shi-Thompson Lower Bound
    #   T   - Thompson Upper Bound

    # FUNCTION:

    # Set Default Method, if no other is selected:
    TypeFlag = TypeFlag[1]
    if (length(method) == 3) method = "T"

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[,1]
        X = table[,2]
        Time = table[,3]
        r = table[,4]
        sigma = table[,5]
    }
    Price = rep(NA, length = length(S))

    # Curran-Thompson Lower Bound:
    if (method == "CT") {
        for ( i in 1:length(S) ) {
        Price[i] = CurranThompsonAsianOption(TypeFlag=TypeFlag,
            S=S[i], X=X[i], Time[i], r=r[i], sigma[i])$price
        }
    }

    # Roger-Shi-Thompson Lower Bound:
    if (method == "RST") {
        for ( i in 1:length(S) ) {
        Price[i] = RogerShiThompsonAsianOption(TypeFlag=TypeFlag,
            S = S[i], X = X[i], Time[i], r = r[i], sigma[i])$price
        }
    }

    # Thompson Upper Bound:
    if (method == "T") {
        for ( i in 1:length(S) ) {
        Price[i] = ThompsonAsianOption(TypeFlag = TypeFlag,
            S = S[i], X = X[i], Time[i], r = r[i], sigma[i])$price
        }
    }

    # Return Value:
    option = list(
        price = Price,
        call = match.call() )
    class(option) = "option"
    option
}


# ------------------------------------------------------------------------------


CurranThompsonAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates "lower bound" for Asian Call Option from
    #   Thompson's formula describing the continuous limit
    #   of Curran's approximation.

    # Note:
    #   Rescale sigma:
    #   Note the formula of Thompson work for Time=1 only!
    #   Thus the easiest way to cover times to maturity different
    #     from unity can be achieved by scale the volatility and
    #     interest rate! - Just do it

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]
    sigma = sigma*sqrt(Time)
    r = r*Time
    Time = 1

    # Settings:
    alpha = r - 0.5*sigma^2

    # Solve for gamma star:
    f1 =  function(x, S, X, alpha, sigma) {
        exp( 3 * ( log(X/S)-alpha/2 ) * x *(1-x/2) +
            alpha * x + 0.5 * sigma^2 * (x-3*x^2*(1-x/2)^2) )
    }

    # Integrate:
    gs = integrate(f1, lower = 0, upper = 1, S = S, X = X, alpha = alpha,
        sigma = sigma, subdivisions = 1000, rel.tol = .Machine$double.eps^0.5,
        abs.tol = .Machine$double.eps^0.5)$value

    # Final Calculate:
    g.star = ( log(2*X/S-gs) - alpha/2 ) / sigma

    # Solve for lower bound:
    f =  function(x, g.star, alpha, sigma) {
        time = x
        arg = (-g.star + sigma*time*(1-time/2))*sqrt(3)
        f = exp( (alpha+sigma^2/2)*time ) * pnorm(arg)
        f
    }

    # Integrate:
    value = integrate(f, lower=0, upper=1, g.star=g.star,
        alpha=alpha, sigma=sigma, subdivisions=1000,
        rel.tol=.Machine$double.eps^0.5,
        abs.tol=.Machine$double.eps^0.5)$value

    # Call Price:
    CurranThompson = exp(-r*Time) * ( S*value - X*pnorm(-g.star*sqrt(3)) )

    # Put Price:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        CurranThompson = CurranThompson - Parity
    }

    # Return Value:
    option = list(
        price = CurranThompson,
        call = match.call() )
    class(option) = "option"
    option
}


# ------------------------------------------------------------------------------


RogerShiThompsonAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates "lower bound" for Asian Call Option from
    #   Thompson's formula. Thompson's result is the same
    #   as can be obtained from Roger and Shi's formula.
    #   However, Thompson's formula is numerically more
    #   efficient since it requires a single integration
    #   only, whereas Roger and Shi's formula requires
    #   double integration.

    # Note:
    #   Rescale sigma:
    #   Note the formula of Thompson work for Time=1 only!
    #   Thus the easiest way to cover times to maturity different
    #     from unity can be achieved by scale the volatility and
    #     interest rate! - Just do it

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]
    sigma = sigma*sqrt(Time)
    r = r*Time
    Time = 1

    # Function from which to calculate gamma star:
    gamma.star =  function(S, X, r, sigma, lower = -99, upper = 99) {
        func =  function(gamma, S, X, r, sigma) {
            f =  function(x, gamma, roh, sigma) {
                time = x
                alpha = roh - sigma*sigma/2
                f = exp( 3 * gamma * sigma * time * (1-time/2) +
                    alpha * time +
                    sigma * sigma * (time-3*time^2*(1-time/2)^2)/2 )
                f
            }
            # Integrate:
            integrate(f, lower = 0, upper = 1, gamma = gamma, roh = r,
                sigma = sigma, subdivisions = 1000,
                rel.tol = .Machine$double.eps^0.5,
                abs.tol = .Machine$double.eps^0.5)$value - X/S
        }
        # Find Root Value:
        uniroot(func, lower = lower, upper = upper, S = S, X = X, r = r,
            sigma = sigma)$root
    }
    g.star = gamma.star(S, X, r, sigma)

    # Function to be integrated:
    f =  function(x, g.star, roh, sigma) {
        time = x
        alpha = roh - sigma*sigma/2
        arg = (-g.star + sigma*time*(1-time/2))*sqrt(3)
        f = exp( (alpha+sigma^2/2)*time ) * pnorm(arg)
        f
    }

    # Integrate:
    value = integrate(f, lower=0, upper=1, g.star=g.star, roh=r,
        sigma=sigma, subdivisions=1000, rel.tol=.Machine$double.eps^0.5,
        abs.tol=.Machine$double.eps^0.5)$value

    # Call Price:
    RogerShiThompson = exp(-r*Time) * ( S*value - X*pnorm(-g.star*sqrt(3)) )

    # Put Price:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        RogerShiThompson = RogerShiThompson - Parity }

    # Return Value:
    option = list(
        price = RogerShiThompson,
        call = match.call() )
    class(option) = "option"
    option
}


# ------------------------------------------------------------------------------


ThompsonAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates "upper bound" for Asian Call Option from
    #   Thompson's formula.

    # Note:
    #   Rescale sigma:
    #   Note the formula of Thompson work for Time=1 only!
    #   Thus the easiest way to cover times to maturity different
    #     from unity can be achieved by scale the volatility and
    #     interest rate! - Just do it

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]
    sigma = sigma*sqrt(Time)
    r = r*Time
    Time = 1

    # Internal Functions:
    sqrtvt =  function(x, S, X, alpha, sigma) {
        t = x
        ct = S*exp(alpha*t)*sigma - X*sigma
        vt = ct^2*t + 2*(X*sigma)*ct*t*(1-t/2) + (X*sigma)^2/3
        sqrt(vt) }

    fmu =  function(t, S, X, r, sigma) {
        alpha = r -sigma^2/2
        gint = integrate(sqrtvt, lower=0, upper=1, S=S, X=X,
            alpha=alpha, sigma=sigma, subdivisions=1000,
            rel.tol=.Machine$double.eps^0.5,
            abs.tol=.Machine$double.eps^0.5)$value
        gamma = (X - S * (exp(alpha)-1) / alpha ) / gint
        mu = (S*exp(alpha*t) +
            gamma*sqrtvt(x=t, S=S, X=X, alpha=alpha, sigma=sigma) ) / X
        mu }

    fatx =  function(t, x, S, X, r, sigma) {
        alpha = r - sigma^2/2
        mu = fmu(t=t, S=S, X=X, r=r, sigma=sigma)
        atx = S*exp(sigma*x+alpha*t) - X*(mu + sigma*x) + X*sigma*(1-t/2)*x
        atx }

    fbtx =  function(t, x, S, X, r, sigma) {
        alpha = r - sigma^2/2
        btx = X*sigma*sqrt(1/3 -t*(1-t/2)^2)
        btx }

    fw =  function(x, v, S2, X2, r2, sigma2) {
        w = x
        atx = fatx(t=v^2, x=w*v, S=S2, X=X2, r=r2, sigma=sigma2)
        btx = fbtx(t=v^2, x=w*v, S=S2, X=X2, r=r2, sigma=sigma2)
        2 * v * dnorm(w) * (atx*pnorm(atx/btx) + btx*dnorm(atx/btx)) }

    fv =  function(x, S1, X1, r1, sigma1) {
        fv = rep(0, length=length(x))
        for (i in 1:length(x))
            fv[i] = integrate(fw, lower=-20, upper=20,
                v=x[i], S2=S1, X2=X1, r2=r1, sigma2=sigma1,
                subdivisions=1000, rel.tol=.Machine$double.eps^0.5,
                abs.tol=.Machine$double.eps^0.5)$value
        exp(-r)*fv }

    # Integrate - Call Price:
    Thompson = integrate(fv, lower=0, upper=1,  S1=S, X1=X, r1=r,
        sigma1=sigma, subdivisions=1000,
        rel.tol=.Machine$double.eps^0.5,
        abs.tol=.Machine$double.eps^0.5)$value

    # Put Price:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        Thompson = Thompson - Parity }

    # Return Value:
    option = list(
        price = Thompson,
        call = match.call() )
    class(option) = "option"
    option
}


# ------------------------------------------------------------------------------


TolmatzAsianOption =
function(TypeFlag = c("c", "p"), S = 100, X = 100, Time = 1, r = 0.09,
sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates "lower bound" for Asian Call Option from
    #   the asymptotic behavior derived by Tolmatz

    # FUNCTION:
    TypeFlag = TypeFlag[1]
    Tolmatz = NA

    # Put Price:
    if (TypeFlag == "p") {
        Parity = (1/(r*Time))*(1-exp(-r*Time))*S - exp(-r*Time)*X
        Tolmatz = Tolmatz - Parity }

    # Return Value:
    option = list(
        price = Tolmatz,
        call = match.call() )
    class(option) = "option"
    option
}


################################################################################
# SYMMETRY AND EQUIVALENCE RELATIONS:


CallPutParityAsianOption =
function(TypeFlag = "p", Price = 8.828759, S = 100, X = 100, Time = 1,
r = 0.09, sigma = 0.3, table = NA)
{   # A function implemented by Diethelm Wuertz

    # Test for Table:
    if (is.data.frame(table)) {
        S = table[,1]
        X = table[,2]
        Time = table[,3]
        r = table[,4]
        sigma = table[,5] }

    # Call from Put:
    if (TypeFlag == "c") {
        # print(Price)
        Parity = S*(1-exp(-r*Time))/(r*Time) - X*exp(-r*Time)
        # print(Parity)
        result = Price + Parity }

    # Put from Call:
    if (TypeFlag == "p") {
        Parity = S*(1-exp(-r*Time))/(r*Time) - X*exp(-r*Time)
        result = Price - Parity }

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


WithDividendsAsianOption =
function(TypeFlag = "c", Dividends = 0.45, S = 100, X = 100, Time = 1,
r = 0.09, sigma = 0.3, calculator = MomentMatchedAsianOption, method = "LN")
{   # A function implemented by Diethelm Wuertz

    # Add Dividends:
    q = Dividends = 0.05
    r.q = r - q
    X.q = X * exp(-q*Time)
    S.q<- S * exp(-q*Time)

    # Call Price:
    if (TypeFlag == "c")
        Price = calculator(TypeFlag = TypeFlag, S = S.q, X = X.q,
            Time = Time, r = r.q, sigma = sigma, method = method)$price
    # Put Price:
    if (TypeFlag == "p" )
        Price = NA

    # Return Value:
    option = list(
        price = Price,
        call = match.call() )
    class(option) = "option"
    option
}


################################################################################
# TABULATED RESULTS:


FuMadanWangTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Fu-Madan-Wang's results from Table 1

    # Source:

    # Settings:
    X = rep(c(90,95,100,105,110), times = 6)
    S = 100*rep(1, times = length(X))
    sigma = rep(0.20, length = length(X))
    r = rep(0.09, length = length(X))

    # Call Prices:
    CallEulerFMW = c(
         11.5293,  7.2131,  3.8087,  1.6465,  0.5761,
         11.9247,  7.7249,  4.3696,  2.1175,  0.8734,
         13.8372,  9.9998,  6.7801,  4.2982,  2.5473,
         17.1212, 13.6763, 10.6319,  8.0436,  5.9267,
         19.8398, 16.6740, 13.7974, 11.2447,  9.0316,
         24.0861, 21.3774, 18.8399, 16.4917, 14.3442)
     CPUEulerFMW = c(
        53,44,42,45,38, 34,33,33,32,32, 26,26,26,25,25,
        24,23,23,23,22, 21,21,21,20,20, 19,22,19,21,21)
     CallPostWidderFMW = c(
         11.5176,  7.1981,  3.8196,  1.6623,  0.5728,
         11.9241,  7.7185,  4.3759,  2.1290,  0.8753,
         13.8439, 10.0029,  6.7823,  4.3010,  2.5450,
         17.1297, 13.6830, 10.6370,  8.0474,  5.9295,
         19.8495, 16.6822, 13.8042, 11.2502,  9.0360,
         24.0861, 21.3774, 18.8399, 16.4917, 14.3442)
     CPUPostWidderFMW = c(
        640,631,628,623,625, 613,591,601,585,595, 518,513,514,503,503,
        475,473,473,469,474, 468,467,467,462,465, 453,442,449,441,453)
     CallZhang = c(
        NA,         NA,         NA,        NA,        NA,
        NA,         NA,         NA,        NA,        NA,
        13.8314996, 9.99566567, 6.7773481, 4.2964626, 2.5462209,
        NA,         NA,         NA,        NA,        NA,
        NA,         NA,         NA,        NA,        NA,
        NA,         NA,         NA,        NA,        NA)

     # Return Value:
     data.frame(cbind(
        S, X, r, sigma,
        CallEulerFMW, CPUEulerFMW, CallPostWidderFMW, CPUPostWidderFMW,
        CallZhang))
}


# ------------------------------------------------------------------------------


FusaiTaglianiTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Fusai and Tagliani's results from Table 6a - 6c

    # Source:
    #   G. Fusai and A. Tagliani [2002]
    #   An Accurate Valuation of Asian Options Using Moments

    # Settings:
    S = 100*rep(1, times = 45)
    X = rep(rep(c(90, 95, 100, 105, 110), times = 3), times = 3)
    Time = rep(1, times = 45)
    r = rep(sort(rep(c(0.05, 0.09, 0.15), times = 5)), times = 3)
    sigma = rep(0.10, times = 15); sigma = c(sigma, 3*sigma, 5*sigma)

    # Moment Matched Call Prices:
    CallLNa = c(
        11.95333,  7.41517,  3.64748,  1.30684,  0.32399,
        13.38629,  8.91721,  4.92310,  2.07045,  0.62338,
        15.39906, 11.12196,  7.03485,  3.61869,  1.41080)
    CallLNb = c(
        14.03812, 10.74879,  7.99251,  5.77417,  4.05724,
        15.06704, 11.73287,  8.88576,  6.54628,  4.69511,
        16.59082, 13.22661, 10.27814,  7.78408,  5.74800)
    CallLNc = c(
        17.67918, 14.91574, 12.48954, 10.38535,  8.58075,
        18.43698, 15.66486, 13.21198, 11.06751,  9.21323,
        19.55391, 16.78229, 14.30234, 12.10904, 10.18997)
    CallFTLN = c(CallLNa, CallLNb, CallLNc)

    # Gram-Charlier Call Prices:
    CallFTa = c(
        11.95127,  7.40754,  3.64091, 1.31097,  0.33156,
        13.38535,  8.91185,  4.91459, 2.07002,  0.63072,
        15.39885, 11.11966,  7.02746, 3.61216,  1.41384)
    CallFTb = c(
        13.89562, 10.63393,  7.92948, 5.76852,  4.09909,
        14.92543, 11.60605,  8.80190, 6.51749,  4.71737,
        16.45856, 13.09056, 10.16897, 7.72184,  5.73774)
    CallFTc = c(
        17.00382, 14.46257, 12.25646, 10.34604,  8.69620,
        17.69986, 15.13557, 12.89937, 10.95396,  9.26553,
        18.73925, 16.14761, 13.87252, 11.88030, 10.13926)
    CallFTLNGC = c(CallFTa, CallFTb, CallFTc)

    # Return Value:
    data.frame(S, X, Time, r, sigma, CallFTLN, CallFTLNGC)
}


# ------------------------------------------------------------------------------


GemanTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Geman's Table

    # Source:
    #   H. Geman [],
    #   Functionals of Brownian Motion in Path Dependent
    #   Option Valuation

    # Settings:
    S = c(1.9, 2.0, 2.1, 2.0, 2.0, 2.0)
    X = rep(2, times = 6)
    Time = c(1, 1, 1, 1, 2, 2)
    r = c(5, 5, 5, 2, 1.25, 5)/100
    sigma = c(5, 5, 5, 1, 2.5, 5)/10

    # Call Prices:
    CallGY = c(
        0.195, 0.248, 0.308, 0.058, 0.1772, 0.352)
    CallMC = c(
        0.191, 0.248, 0.306, 0.056, 0.1771, 0.347)

    # Return Value:
    data.frame(cbind(S, X, Time, r, sigma, CallGY, CallMC))
}


# ------------------------------------------------------------------------------


LinetzkyTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Linetzky's Table 3

    # Source:
    #   V. Linetzky [2002]
    #   Spectral Expansions for Asian (Average Price) Options

    # Settings:
    S     = c(2.00, 2.00, 2.00,   1.90, 2.00, 2.10, 2.00)
    X     = c(2.00, 2.00, 2.00,   2.00, 2.00, 2.00, 2.00)
    Time  = c(1.00, 1.00, 2.00,   1.00, 1.00, 1.00, 2.00)
    r     = c(0.02, 0.18, 0.0125, 0.05, 0.05, 0.05, 0.05)
    sigma = c(10.0, 30.0, 25.0,   50.0, 50.0, 50.0, 50.0) / 100

    # Call Prices:
    CallEE = c(
        0.0559860415,  0.2183875466,  0.1722687410,  0.1931737903,
        0.2464156905,  0.3062203648,  0.3500952190)
    CallSLT = c(
        0.055986,      0.218388,      0.172269,      0.193174,
        0.246416,      0.306220,      0.350095)
    CallMC = c(
        0.05602,       0.2185,        0.1725,        0.1933,
        0.2465,        0.3064,        0.3503)
    CallTLB = c(
        0.055985,      0.218366,      0.172226,      0.193060,
        0.246298,      0.306094,      0.349779)
    CallTUB = c(
        0.055989,      0.218473,      0.172451,      0.193799,
        0.247054,      0.306904,      0.352556)

    # Return Value:
    data.frame(S, X, Time, r, sigma, CallEE, CallSLT, CallMC, CallTLB, CallTUB)
}


# ------------------------------------------------------------------------------


ZhangTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Zhangs's results from Table 1

    # Source:
    #   J.E. Zhang [2002],
    #   A semi-analytical method for pricing and hedging
    #   continuously sampled arithmetic average rate options

    # Settings:
    S = 100*rep(1, times = 36)
    X = c(rep(c(95, 100, 105), times = 3), rep(c(90, 100, 110), times = 9))
    Time = rep(1, times = 36)
    r = rep(c(5, 5, 5, 9, 9, 9, 15, 15, 15)/100, times = 4)
    sigma = sort(rep(c(0.05, 0.10, 0.20, 0.30), times = 9))

    # Call Prices:
    CallZ = c(
         7.1777275, 2.7161745, 0.3372614,  8.8088302,  4.3082350, 0.9583841,
        11.0940944, 6.7943550, 2.7444531, 11.9510927,  3.6413864, 0.3312030,
        13.3851974, 4.9151167, 0.6302713, 15.3987687,  7.0277081, 1.4136149,
        12.5959916, 5.7630881, 1.9898945, 13.8314996,  6.7773481, 2.5462209,
        15.6417575, 8.4088330, 3.5556100, 13.9538233,  7.9456288, 4.0717942,
        14.9839595, 8.8287588, 4.6967089, 16.5129113, 10.2098305, 5.7301225)

    # Return Value:
    data.frame(cbind(S, X, Time, r, sigma, CallZ))
}


# ------------------------------------------------------------------------------


ZhangShortTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Zhangs's short-tenor (ST) results from Table 7

    # Source:
    #   J.E. Zhang [2002],
    #   A semi-analytical method for pricing and hedging
    #   continuously sampled arithmetic average rate options

    # Settings:
    S = 100*rep(1, times = 18)
    X = c(rep(c(95,100,105), times = 6))
    Time = rep(0.08, times = 18)
    r = rep(0.09, times = 18)
    sigma = sort(rep(c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50), times = 3))

    # Call Prices:
    CallPM = c(
         5.3224059, 0.5343219, 0.0000000,  5.3225549, 0.8431357, 0.0014921,
         5.3816262, 1.4835707, 0.1292560,  5.6304554, 2.1291126, 0.4880727,
         6.0222502, 2.7758194, 0.9753549,  6.4947818, 3.4228569, 1.5274729)
    CallCT = c(
         5.3224059, 0.5343040, 0.0000000,  5.3225550, 0.8431302, 0.0014924,
         5.3816340, 1.4835272, 0.1292529,  5.6304480, 2.1289662, 0.4879999,
         6.0221516, 2.7754740, 0.9750983,  6.4944850, 3.4221863, 1.5268836)

    # Return Value:
    data.frame(cbind(S, X, Time, r, sigma, CallPM, CallCT))
}


# ------------------------------------------------------------------------------


ZhangLongTable =
function()
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Display Zhangs's long-tenor (LT) results from Table 7

    # Source:
    #   J.E. Zhang [2002],
    #   A semi-analytical method for pricing and hedging
    #   continuously sampled arithmetic average rate options

    # Settings:
    S = 100*rep(1, times = 18)
    X = c(rep(c(95,100,105), times = 6))
    Time = rep(3, times = 18)
    r = rep(0.09, times = 18)
    sigma = sort(rep(c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50), times = 3))

    # Call Prices:
    CallPM = c(
         15.11626, 11.30359,  7.55327,   15.21331, 11.63716,  8.39113,
         16.63363, 13.76592, 11.22170,   19.01304, 16.58331, 14.39723,
         21.70848, 19.57038, 17.62156,   24.47434, 22.55837, 20.79507)
    CallCT = c(
         15.11626, 11.30361,  7.55328,   15.21376, 11.63752,  8.39084,
         16.63611, 13.76547, 11.21783,   19.01856, 16.58101, 14.38708,
         21.72991, 19.57593, 17.61228,   24.54788, 22.60645, 20.81811)

    # Return Value:
    data.frame(cbind(S, X, Time, r, sigma, CallPM, CallCT))
}


################################################################################

