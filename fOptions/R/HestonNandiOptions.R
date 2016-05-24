
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


################################################################################
# FUNCTION:           DESCRIPTION:
#  HNGOption           Computes Option Price from the HN-GARCH Formula
#  HNGGreeks           Calculates one of the Greeks of the HN-GARCH Formula
#  HNGCharacteristics  Computes Option Price and all Greeks of HN-GARCH Model
################################################################################


HNGOption =
function(TypeFlag = c("c", "p"), model, S, X, Time.inDays, r.daily)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the price of a HN-GARCH option.

    # Details:
    #   The function calculates the price of a Heston-Nandi GARCH(1,1)
    #   call or put option.

    # FUNCTION:

    # Option Type:
    TypeFlag = TypeFlag[1]

    # Integrate:
    call1 = integrate(.fstarHN, 0, Inf, const = 1, model = model,
        S = S, X = X, Time.inDays = Time.inDays, r.daily = r.daily)
    # For SPlus Compatibility:
    if (is.null(call1$value)) call1$value = call1$integral
    call2 = integrate(.fstarHN, 0, Inf, const = 0, model = model,
        S = S, X = X, Time.inDays = Time.inDays, r.daily = r.daily)
    # For SPlus Compatibility:
    if (is.null(call2$value)) call2$value = call2$integral

    # Compute Call Price:
    call.price = S/2 + exp(-r.daily*Time.inDays) * call1$value -
        X * exp(-r.daily*Time.inDays) * ( 1/2 + call2$value )

    # Select Option Price:
    price = NA
    if (TypeFlag == "c" ) price = call.price
    if (TypeFlag == "p" ) price = call.price + X*exp(-r.daily*Time.inDays) - S

    # Return Value:
    option = list(
        price = price,
        call = match.call())
    class(option) = "option"
    option
}


.fstarHN <-
function(phi, const, model, S, X, Time.inDays, r.daily)
{
    # Internal Function:

    # Model Parameters:
    lambda = -1/2
    omega = model$omega
    alpha = model$alpha
    gamma = model$gamma + model$lambda + 1/2
    beta = model$beta
    sigma2 = (omega + alpha)/(1 - beta - alpha * gamma^2)
    # Function to be integrated:
    cphi0 = phi*complex(real = 0, imaginary = 1)
    cphi = cphi0 + const
    a = cphi * r.daily
    b = lambda*cphi + cphi*cphi/2
    for (i in 2:Time.inDays) {
        a = a + cphi*r.daily + b*omega - log(1-2*alpha*b)/2
        b = cphi*(lambda+gamma) - gamma^2/2 + beta*b +
            0.5*(cphi-gamma)^2/(1-2*alpha*b) }
    f = Re(exp(-cphi0*log(X)+cphi*log(S)+a+b*sigma2 )/cphi0)/pi

    # Return Value:
    f
}


# ------------------------------------------------------------------------------


HNGGreeks =
function(Selection = c("Delta", "Gamma"), TypeFlag = c("c", "p"), model,
S, X, Time.inDays, r.daily)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the Greeks of a HN-GARCH option.

    # Details:
    #   The function calculates the delta and gamma Greeks of
    #   a Heston Nandi GARCH(1,1) call or put option.

    # FUNCTION:

    # Type Flags:
    Selection = Selection[1]
    TypeFlag = TypeFlag[1]

    # Delta:
    if (Selection == "Delta") {
        # Integrate:
        delta1 = integrate(.fdeltaHN, 0, Inf, const = 1, model = model,
            S = S, X = X, Time.inDays = Time.inDays, r.daily = r.daily)
        # For SPlus Compatibility:
        if (is.null(delta1$value)) delta1$value = delta1$integral
        delta2 = integrate(.fdeltaHN, 0, Inf, const = 0, model = model,
            S = S, X = X, Time.inDays = Time.inDays, r.daily = r.daily)
        # For SPlus Compatibility:
        if (is.null(delta2$value)) delta2$value = delta2$integral
        # Compute Call and Put Delta :
        greek = 1/2 + exp(-r.daily*Time.inDays) * delta1$value -
            X * exp(-r.daily*Time.inDays) * delta2$value
        if (TypeFlag == "p") greek = greek - 1 }

    # Gamma:
    if (Selection == "Gamma") {
        # Integrate:
        gamma1 = integrate(.fgammaHN, 0, Inf, const = 1, model = model,
            S = S, X = X, Time.inDays = Time.inDays, r.daily = r.daily)
        # For SPlus Compatibility:
        if (is.null(gamma1$value)) gamma1$value = gamma1$integral
        gamma2 = integrate(.fgammaHN, 0, Inf, const = 0, model = model,
            S = S, X = X, Time.inDays = Time.inDays, r.daily = r.daily)
        # For SPlus Compatibility:
        if (is.null(gamma2$value)) gamma2$value = gamma2$integral
        # Compute Call and Put Gamma :
        greek = put.gamma = exp(-r.daily*Time.inDays) * gamma1$value -
            X * exp(-r.daily*Time.inDays) * gamma2$value }

    # Return Value:
    greek
}


.fdeltaHN <-
function(phi, const, model, S, X, Time.inDays, r.daily)
{
    # Function to be integrated:
    cphi0 = phi * complex(real = 0, imaginary = 1)
    cphi = cphi0 + const
    fdelta = cphi *
        .fHN(phi, const, model, S, X, Time.inDays, r.daily) / S

    # Return Value:
    Re(fdelta)
}


.fgammaHN <-
function(phi, const, model, S, X, Time.inDays, r.daily)
{
    # Function to be integrated:
    cphi0 = phi * complex(real = 0, imaginary = 1)
    cphi = cphi0 + const
    fgamma = cphi * ( cphi - 1 ) *
        .fHN(phi, const, model, S, X, Time.inDays, r.daily) / S^2

    # Return Value:
    Re(fgamma)
}



.fHN <-
function(phi, const, model, S, X, Time.inDays, r.daily)
{
    # Internal Function:

    # Model Parameters:
    lambda = -1/2
    omega = model$omega
    alpha = model$alpha
    gamma = model$gamma + model$lambda + 1/2
    beta = model$beta
    sigma2 = (omega + alpha)/(1 - beta - alpha * gamma^2)
    # Function to be integrated:
    cphi0 = phi*complex(real = 0, imaginary = 1)
    cphi = cphi0 + const
    a = cphi * r.daily
    b = lambda*cphi + cphi*cphi/2
    for (i in 2:Time.inDays) {
        a = a + cphi*r.daily + b*omega - log(1-2*alpha*b)/2
        b = cphi*(lambda+gamma) - gamma^2/2 + beta*b +
            0.5*(cphi-gamma)^2/(1-2*alpha*b) }
    fun = exp(-cphi0*log(X)+cphi*log(S)+a+b*sigma2)/cphi0/pi

    # Return Value:
    fun
}


# ------------------------------------------------------------------------------


HNGCharacteristics =
function(TypeFlag = c("c", "p"), model, S, X, Time.inDays, r.daily)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   The function calculates the option price for the Heston
    #   Nandi Garch(1,1) option model together with the delta
    #   and gamma option sensitivies.

    # FUNCTION:

    # Premium and Function Call to all Greeks
    TypeFlag = TypeFlag[1]
    premium = HNGOption(TypeFlag, model, S, X, Time.inDays, r.daily)
    delta = HNGGreeks("Delta", TypeFlag, model, S, X, Time.inDays, r.daily)
    gamma = HNGGreeks("Gamma", TypeFlag, model, S, X, Time.inDays, r.daily)

    # Return Value:
    list(premium = premium, delta = delta, gamma = gamma)
}


################################################################################

