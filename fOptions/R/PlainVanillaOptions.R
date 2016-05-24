
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
# FUNCTION:                  DESCRIPTION:
#  'fOPTION'                  S4 Class Representation
# FUNCTION:                  DESCRIPTION:
#  NDF                        Normal distribution function
#  CND                        Cumulative normal distribution function
#  CBND                       Cumulative bivariate normal distribution
# FUNCTION:                  DESCRIPTION:
#  GBSOption                  Computes Option Price from the GBS Formula
#  GBSCharacteristics         Computes Option Price and all Greeks of GBS Model
#   BlackScholesOption         Synonyme Function Call to GBSOption
#  GBSGreeks                  Computes one of the Greeks of the GBS formula
# FUNCTION:                  DESCRIPTION:
#  Black76Option              Computes Prices of Options on Futures
#  MiltersenSchwartzOption    Pricing a Miltersen Schwartz Option
# S3 METHODS:                DESCRIPTION:
#  print.option               Print Method
#  summary.otion              Summary Method
################################################################################


setClass("fOPTION",
    representation(
        call = "call",
        parameters = "list",
        price = "numeric",
        title = "character",
        description = "character"
    )
)


################################################################################


NDF =
function(x)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the normal distribution function.

    # FUNCTION:

    # Compute:
    result = exp(-x*x/2)/sqrt(8*atan(1))

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


CND =
function(x)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the cumulated normal distribution function.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Compute:
    k  = 1 / ( 1 + 0.2316419 * abs(x) )
    a1 =  0.319381530; a2 = -0.356563782; a3 =  1.781477937
    a4 = -1.821255978; a5 =  1.330274429
    result = NDF(x) * (a1*k + a2*k^2 + a3*k^3 + a4*k^4 + a5*k^5) - 0.5
    result = 0.5 - result*sign(x)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


CBND =
function(x1, x2, rho)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the cumulative bivariate normal distribution function.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Compute:
    # Take care for the limit rho = +/- 1
    a = x1
    b = x2
    if (abs(rho) == 1) rho = rho - (1e-12)*sign(rho)
    # cat("\n a - b - rho :"); print(c(a,b,rho))
    X = c(0.24840615, 0.39233107, 0.21141819, 0.03324666, 0.00082485334)
    y = c(0.10024215, 0.48281397, 1.0609498, 1.7797294, 2.6697604)
    a1 = a / sqrt(2 * (1 - rho^2))
    b1 = b / sqrt(2 * (1 - rho^2))
    if (a <= 0 && b <= 0 && rho <= 0) {
       Sum1 = 0
       for (I in 1:5) {
            for (j in 1:5) {
            Sum1 = Sum1 + X[I] * X[j] *
              exp(a1*(2*y[I]-a1) + b1*(2*y[j]-b1) +
              2*rho*(y[I]-a1)*(y[j]-b1)) } }
       result = sqrt(1 - rho^2) / pi * Sum1
       return(result) }
    if (a <= 0 && b >= 0 && rho >= 0) {
        result = CND(a) - CBND(a, -b, -rho)
        return(result) }
    if (a >= 0 && b <= 0 && rho >= 0) {
        result = CND(b) - CBND(-a, b, -rho)
        return(result) }
    if (a >= 0 && b >= 0 && rho <= 0) {
        result = CND(a) + CND(b) - 1 + CBND(-a, -b, rho)
        return(result) }
    if (a * b * rho >= 0 ) {
        rho1 = (rho*a - b) * sign(a) / sqrt(a^2 - 2*rho*a*b + b^2)
        rho2 = (rho*b - a) * sign(b) / sqrt(a^2 - 2*rho*a*b + b^2)
        delta = (1 - sign(a) * sign(b)) / 4
        result = CBND(a, 0, rho1) + CBND(b, 0, rho2) - delta
        return(result) }

    # Return Value:
    invisible()
}


# ******************************************************************************


GBSOption =
function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, title = NULL,
description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the Generalized Black-Scholes option
    #   price either for a call or a put option.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Compute:
    TypeFlag = TypeFlag[1]
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    d2 = d1 - sigma*sqrt(Time)
    if (TypeFlag == "c")
        result = S*exp((b-r)*Time)*CND(d1) - X*exp(-r*Time)*CND(d2)
    if (TypeFlag == "p")
        result = X*exp(-r*Time)*CND(-d2) - S*exp((b-r)*Time)*CND(-d1)

    # Parameters:
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Black Scholes Option Valuation"
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


GBSCharacteristics =
function(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the Options Characterisitics (Premium
    #   and Greeks for a Generalized Black-Scholes option
    #   either for a call or a put option.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Premium and Function Call to all Greeks
    TypeFlag = TypeFlag[1]
    premium = GBSOption(TypeFlag, S, X, Time, r, b, sigma)@price
    delta = GBSGreeks("Delta", TypeFlag, S, X, Time, r, b, sigma)
    theta = GBSGreeks("Theta", TypeFlag, S, X, Time, r, b, sigma)
    vega = GBSGreeks("Vega", TypeFlag, S, X, Time, r, b, sigma)
    rho = GBSGreeks("Rho", TypeFlag, S, X, Time, r, b, sigma)
    lambda = GBSGreeks("Lambda", TypeFlag, S, X, Time, r, b, sigma)
    gamma = GBSGreeks("Gamma", TypeFlag, S, X, Time, r, b, sigma)

    # Return Value:
    list(premium = premium, delta = delta, theta = theta,
        vega = vega, rho = rho, lambda = lambda, gamma = gamma)
}


# ------------------------------------------------------------------------------


BlackScholesOption =
function(...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   A synonyme for GBSOption

    # FUNCTION:

    # Return Value:
    GBSOption(...)
}



# ******************************************************************************


GBSGreeks =
function(Selection = c("Delta", "Theta", "Vega", "Rho", "Lambda", "Gamma",
"CofC"), TypeFlag = c("c", "p"), S, X, Time, r, b, sigma)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate the Options Greeks for a Generalized
    #   Black-Scholes option either for a call or a put option.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]
    Selection = Selection[1]

    # Function Call to all Greeks via selection parameter
    result = NA
    if (Selection == "Delta" || Selection == "delta")
            result = .GBSDelta (TypeFlag, S, X, Time, r, b, sigma)
    if (Selection == "Theta" || Selection == "theta")
            result = .GBSTheta (TypeFlag, S, X, Time, r, b, sigma)
    if (Selection == "Vega" || Selection == "vega")
            result = .GBSVega  (TypeFlag, S, X, Time, r, b, sigma)
    if (Selection == "Rho" || Selection == "rho")
            result = .GBSRho   (TypeFlag, S, X, Time, r, b, sigma)
    if (Selection == "Lambda" || Selection == "lambda")
            result = .GBSLambda(TypeFlag, S, X, Time, r, b, sigma)
    if (Selection == "Gamma" || Selection == "gamma")
            result = .GBSGamma (TypeFlag, S, X, Time, r, b, sigma)
    if (Selection == "CofC" || Selection == "cofc")
            result = .GBSCofC  (TypeFlag, S, X, Time, r, b, sigma)

    # Return Value:
    result
}


# Internal Functions:

.GBSDelta <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    if (TypeFlag == "c") result = exp((b-r)*Time)*CND(d1)
    if (TypeFlag == "p") result = exp((b-r)*Time)*(CND(d1)-1)
    result
}


.GBSTheta <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    d2 = d1 - sigma*sqrt(Time)
    Theta1 = -(S*exp((b-r)*Time)*NDF(d1)*sigma)/(2*sqrt(Time))
    if (TypeFlag == "c") result = Theta1 -
        (b-r)*S*exp((b-r)*Time)*CND(+d1) - r*X*exp(-r*Time)*CND(+d2)
    if (TypeFlag == "p") result = Theta1 +
        (b-r)*S*exp((b-r)*Time)*CND(-d1) + r*X*exp(-r*Time)*CND(-d2)
    result
}


.GBSVega <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    result = S*exp((b-r)*Time)*NDF(d1)*sqrt(Time) # Call,Put
    result
}


.GBSRho <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    d2 = d1 - sigma*sqrt(Time)
    CallPut = GBSOption(TypeFlag, S, X, Time, r, b , sigma)@price
    if (TypeFlag == "c") {
        if (b != 0) {result =  Time * X * exp(-r*Time)*CND( d2)}
        else {result = -Time * CallPut } }
    if (TypeFlag == "p") {
        if (b != 0) {result = -Time * X * exp(-r*Time)*CND(-d2)}
        else { result = -Time * CallPut } }
    result
}


.GBSLambda <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    CallPut = GBSOption(TypeFlag,S,X,Time,r,b,sigma)@price
    if (TypeFlag == "c") result = exp((b-r)*Time)* CND(d1)*S / CallPut
    if (TypeFlag == "p") result = exp((b-r)*Time)*(CND(d1)-1)*S / CallPut
    result
}


.GBSGamma <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    result = exp((b-r)*Time)*NDF(d1)/(S*sigma*sqrt(Time)) # Call,Put
    result
}


.GBSCofC <-
function(TypeFlag, S, X, Time, r, b, sigma)
{
    d1 = ( log(S/X) + (b+sigma*sigma/2)*Time ) / (sigma*sqrt(Time))
    if (TypeFlag == "c") result = Time*S*exp((b-r)*Time)*CND(d1)
    if (TypeFlag == "p") result = -Time*S*exp((b-r)*Time)*CND(-d1)
    result }


# ------------------------------------------------------------------------------


Black76Option =
function(TypeFlag = c("c", "p"), FT, X, Time, r, sigma, title = NULL,
description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculate Options Price for Black (1977) Options
    #   on futures/forwards

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Settings:
    TypeFlag = TypeFlag[1]

    # Result:
    result = GBSOption(TypeFlag = TypeFlag, S = FT, X = X, Time = Time,
        r = r, b = 0, sigma = sigma)@price

    # Parameters:
    param = list()
    param$TypeFlag = TypeFlag
    param$FT = FT
    param$X = X
    param$Time = Time
    param$r = r
    param$sigma = sigma

    # Add title and description:
    if (is.null(title)) title = "Black 76 Option Valuation"
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


# ******************************************************************************


MiltersenSchwartzOption =
function (TypeFlag = c("c", "p"), Pt, FT, X, time, Time, sigmaS, sigmaE,
sigmaF, rhoSE, rhoSF, rhoEF, KappaE, KappaF, title = NULL,
description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Miltersen Schwartz (1997) commodity option model.

    # References:
    #   Haug E.G., The Complete Guide to Option Pricing Formulas

    # FUNCTION:

    # Settings:
    TyoeFlag = TypeFlag[1]

    # Compute:
    vz = sigmaS^2*time+2*sigmaS*(sigmaF*rhoSF*1/KappaF*(time-1/KappaF*
      exp(-KappaF*Time)*(exp(KappaF*time)-1))-sigmaE*rhoSE*1/KappaE*
      (time-1/KappaE*exp(-KappaE*Time)*(exp(KappaE*time)-1)))+sigmaE^2*
      1/KappaE^2*(time+1/(2*KappaE)*exp(-2*KappaE*Time)*(exp(2*KappaE*time)-
      1)-2*1/KappaE*exp(-KappaE*Time)*(exp(KappaE*time)-1))+sigmaF^2*
      1/KappaF^2*(time+1/(2*KappaF)*exp(-2*KappaF*Time)*(exp(2*KappaF*time)-
      1)-2*1/KappaF*exp(-KappaF*Time)*(exp(KappaF*time)-1))-2*sigmaE*
      sigmaF*rhoEF*1/KappaE*1/KappaF*(time-1/KappaE*exp(-KappaE*Time)*
      (exp(KappaE*time)-1)-1/KappaF*exp(-KappaF*Time)*(exp(KappaF*time)-
      1)+1/(KappaE+KappaF)*exp(-(KappaE+KappaF)*Time)*(exp((KappaE+KappaF)*
      time)-1))
    vxz = sigmaF*1/KappaF*(sigmaS*rhoSF*(time-1/KappaF*(1-exp(-KappaF*
      time)))+sigmaF*1/KappaF*(time-1/KappaF*exp(-KappaF*Time)*(exp(KappaF*
      time)-1)-1/KappaF*(1-exp(-KappaF*time))+1/(2*KappaF)*exp(-KappaF*
      Time)*(exp(KappaF*time)-exp(-KappaF*time)))-sigmaE*rhoEF*1/KappaE*
      (time-1/KappaE*exp(-KappaE*Time)*(exp(KappaE*time)-1)-1/KappaF*(1-
      exp(-KappaF*time))+1/(KappaE+KappaF)*exp(-KappaE*Time)*
      (exp(KappaE*time)-exp(-KappaF*time))))
    vz = sqrt(vz)
    d1 = (log(FT/X)-vxz+vz^2/2)/vz
    d2 = (log(FT/X)-vxz-vz^2/2)/vz

    # Call/Put:
    if (TypeFlag == "c") {
        result = Pt*(FT*exp(-vxz)*CND(d1)-X*CND(d2)) }
    if (TypeFlag == "p") {
        result = Pt*(X*CND(-d2)-FT*exp(-vxz)*CND(-d1)) }

    # Parameters:
    param = list()
    param$TypeFlag = TypeFlag
    param$Pt = Pt
    param$FT = FT
    param$X = X
    param$time = time
    param$Time = Time
    param$sigmaS = sigmaS
    param$sigmaE = sigmaE
    param$sigmaF = sigmaF
    param$rhoSE = rhoSE
    param$rhoSF = rhoSF
    param$rhoEF = rhoEF
    param$KappaE = KappaE
    param$KappaF = KappaF

    # Add title and description:
    if (is.null(title)) title = "Miltersen Schwartz Option Valuation"
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


# ******************************************************************************


GBSVolatility = function(price, TypeFlag = c("c", "p"), S, X, Time, r, b,
tol = .Machine$double.eps, maxiter = 10000)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute implied volatility

    # Example:
    #   sigma = GBSVolatility(price=10.2, "c", S=100, X=90, Time=1/12, r=0, b=0)
    #   sigma
    #   GBSOption("c", S=100, X=90, Time=1/12, r=0, b=0, sigma=sigma)@price

    # FUNCTION:

    # Option Type:
    TypeFlag = TypeFlag[1]

    # Search for Root:
    volatility = uniroot(.fGBSVolatility, interval = c(-10,10), price = price,
        TypeFlag = TypeFlag, S = S, X = X, Time = Time, r = r, b = b,
        tol = tol, maxiter = maxiter)$root

    # Return Value:
    volatility
}



# Internal Function:
.fGBSVolatility <-
function(x, price, TypeFlag, S, X, Time, r, b, ...)
{
    GBS = GBSOption(TypeFlag = TypeFlag, S = S, X = X, Time = Time,
        r = r, b = b, sigma = x)@price
    price - GBS
}


# ------------------------------------------------------------------------------


print.option =
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Print method for objects of class "option".

    # FUNCTION:

    # Print Method:
    object = x
    cat("\nCall:", deparse(object$call), "", sep = "\n")
    cat("Option Price:\n")
    cat(object$price, "\n")
}


# ------------------------------------------------------------------------------


summary.option =
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Summary method for objects of class "option".

    # FUNCTION:

    # Summary Method:
    print(object, ...)
}


################################################################################



setMethod("show", "fOPTION",
    function(object)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Print method for objects of class "fOPTION".

    # FUNCTION:

    # Print Method:
    Parameter = unlist(object@parameters)
    Names = names(Parameter)
    Parameter = cbind(as.character(Parameter))
    rownames(Parameter) = paste("", Names)
    colnames(Parameter) = "Value:"

    # Title:
    cat("\nTitle:\n ")
    cat(object@title, "\n")

    # Call:
    cat("\nCall:", paste("", deparse(object@call)), "", sep = "\n")

    # Parameters:
    cat("Parameters:\n")
    print(Parameter, quote = FALSE)

    # Price:
    cat("\nOption Price:\n ")
    cat(object@price, "\n")

    # Description:
    cat("\nDescription:\n ")
    cat(object@description, "\n\n")

    # Return Value:
    invisible()
})


# ------------------------------------------------------------------------------


summary.fOPTION =
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Summary method for objects of class "option".

    # FUNCTION:

    # Summary Method:
    print(object, ...)
}


################################################################################

