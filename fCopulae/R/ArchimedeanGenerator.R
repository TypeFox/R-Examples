
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
# FUNCTION:                 ARCHIMEDEAN COPULAE PARAMETER:
#  archmList                 Returns list of implemented Archimedean copulae
#  archmParam                Sets Default parameters for an Archimedean copula
#  archmRange                Returns the range of valid alpha values
#  archmCheck                Checks if alpha is in the valid range
# FUNCTION:                 ARCHIMEDEAN COPULAE PHI GENERATOR:
#  Phi                       Computes Archimedean Phi, inverse and derivatives
#  PhiSlider                 Displays interactively generator function
#  .Phi                      Computes Archimedean generator Phi
#  .Phi0                     Utility Function
#  .PhiFirstDer              Computes first derivative of Phi
#  .PhiSecondDer             Computes second derivative of Phi
#  .invPhi                   Computes inverse of Archimedean generator
#  .invPhiFirstDer           Computes first derivative of inverse Phi
#  .invPhiSecondDer          Computes second derivative of inverse Phi
# FUNCTION:                 ARCHIMEDEAN DENSITY K GENERATOR:
#  Kfunc                     Computes Archimedean Density Kc and its Inverse
#  KfuncSlider               Displays interactively the density and concordance
#  .Kfunc                    Computes Density for Archimedean Copulae
#  .invK                     Computes Inverse of Density
#  .invK2                    Utility Function
#  .ALPHA                    Utility Function
#  .TAU                      Utility Function
#  .RHO                      Utility Function
################################################################################


################################################################################
# FUNCTION:                  ARCHIMEDEAN COPULAE PARAMETER:
#  archmList                  Returns list of implemented Archimedean copulae
#  archmParam                 Sets default parameters for an Archimedean copula
#  archmCheck                 Checks if alpha is in the valid range
#  archmRange                 Returns the range of valid alpha values


archmList <- 
    function()
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns list of implemented Archimedean copulae

    # Compose List:
    ans <- paste(1:22)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


archmParam <- 
    function(type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets default parameters for Archimedean copulae

    # Arguments:
    #   type - a character string or integer value naming the copula.
    #       By default the first copula will be chosen.

    # Value:
    #   returns a list with two elements, 'param' sets the parameters
    #       which may be a vector, 'range' the range with minimum and
    #       maximum values for each of the parameters.

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Parameter Values:
    B = Inf
    lower=c(-1, 1,-1, 1,-B, 1,  0, 1, 0, 0, 0, 1, 0, 1, 1, 0,-B, 2, 0, 0, 1, 0)
    upper=c( B, B, 1, B, B, B,  1, B, 1, 1,.5, B, B, B, B, B, B, B, B, B, B, 1)
    Alpha=c( 1, 2,.5, 2, 1, 2, .5, 2,.5,.5,.2, 2, 1, 2, 2, 1,.5, 3, 1, 1, 2,.5)

    # Parameter List:
    ans = list(copula = type)
    ans$param = c(alpha = Alpha[Type])
    ans$range = c(lower = lower[Type], upper = upper[Type])

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


archmRange <- 
    function(type = archmList(), B = Inf)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the range of valid alpha values

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Range:
    lower = c(-1, 1,-1, 1,-B, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0,-B, 2, 0, 0, 1, 0)
    upper = c( B, B, 1, B, B, B, 1, B, 1, 1,.5, B, B, B, B, B, B, B, B, B, B, 1)

    # Return Value:
    ans = cbind(lower[Type], upper[Type])
    rownames(ans) = type
    colnames(ans) = c("lower", "upper")
    ans
}


# ------------------------------------------------------------------------------


archmCheck <- 
    function(alpha, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks if alpha is in the valid range

    # FUNCTION:

    # Type:
    type = match.arg(type)

    # Check:
    ans = TRUE
    range = as.vector(archmRange(type))
    if (alpha < range[1] | alpha > range[2]) {
        print(c(alpha = alpha))
        print(c(range = range))
        stop("alpha is out of range")
    }

    # Return Value:
    invisible(TRUE)
}


################################################################################
# FUNCTION:                  ARCHIMEDEAN COPULAE PHI GENERATOR:
#  Phi                        Computes Archimedean Phi, inverse and derivatives
#  PhiSlider                  Displays interactively generator function
#  .Phi                       Computes Archimedean generator Phi
#  .Phi0                      Utility Function
#  .PhiFirstDer               Computes first derivative of Phi
#  .PhiSecondDer              Computes second derivative of Phi
#  .invPhi                    Computes inverse of Archimedean generator
#  .invPhiFirstDer            Computes first derivative of inverse Phi
#  .invPhiSecondDer           Computes second derivative of inverse Phi


Phi <- 
    function(x, alpha = NULL, type = archmList(), inv = FALSE, 
    deriv = paste(0:2))
{   
    # A function implemented by Diethelm Wuertz

    # Type:
    type = match.arg(type)
    Type = as.integer(type)
    deriv = match.arg(deriv)

    # Default alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # Phi Generator:
    if (inv) {
        if (deriv == "0") {
            ans = .invPhi(x, alpha, type)
            names(ans) = "invPhi"
        }
        if (deriv == "1") {
            ans = .invPhiFirstDer(x, alpha, type)
            names(ans) = "invPhiFirstDer"
        }
        if (deriv == "2") {
            ans = .invPhiSecondDer(x, alpha, type)
            names(ans) = "invPhiSecondDer"
        }
    } else {
        if (deriv == "0") {
            ans = .Phi(x, alpha, type)
            names(ans) = "Phi"
        }
        if (deriv == "1") {
            ans = .PhiFirstDer(x, alpha, type)
            names(ans) = "PhiFirstDer"
        }
        if (deriv == "2") {
            ans = .PhiSecondDer(x, alpha, type)
            names(ans) = "PhiSecondDer"
        }
    }

    # Add Control Attribute:
    attr(ans, "control")<-cbind.data.frame(alpha = alpha, type = type,
        inv = inv, deriv = deriv, row.names = "")

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


PhiSlider <- 
    function(B = 5)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively the dependence function

    # FUNCTION:

    # Graphic Frame:
    par(mfcol = c(2, 2), cex = 0.7)

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 10) return ()

        # Sliders:
        Copula = as.integer(.sliderMenu(no = 1))
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)

        # Skip:
        if (Copula == 13 & alpha == 0) return(invisible())

        # Do we have a strict Copula?
        strict = c(
            "Yes","No","Yes","Yes","Yes","Yes","No","No","Yes","Yes",
            "No","Yes","Yes","Yes","No","Yes","Yes","No","Yes","Yes",
            "No","Yes")[Copula]
        if (alpha < 0 & Copula == 1) strict[1] = "No"
        if (alpha == 0 & Copula == 16) strict[16] = "No"

        # What is the Range?
        RANGE = c(
            "-1|Inf", "1|Inf", "-1|1", "-Inf|inf", "0|1", "0|0.5",
            "0|Inf", "2|Inf")[No]

        # Which one is the Limit Copula?
        limitTitle = rep("NA", times = 22)
        if (alpha == -1)
            limitTitle = c(
                "W ", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "NA", "NA", "NA", "NA", "NA", "NA", "Pi", "NA", "NA", "NA",
                "NA", "NA")
        if (alpha == 0)
            limitTitle = c(
                "Pi", "NA", "Pi", "NA", "Pi", "NA", "W ", "NA", "Pi", "Pi",
                "Pi", "NA", "NA", "NA", "NA", "W ", "NA", "NA", "L ", "Pi",
                "NA", "Pi")
        if (alpha == 1)
            limitTitle = c(
                "L ", "W ", "L ", "Pi", "NA", "Pi", "Pi", "W ", "NA", "NA",
                "NA", "L ", "Pi", "L ", "W ", "NA", "NA", "NA", "NA", "NA",
                "W ", "NA")
        limitTitle = limitTitle[Copula]
        if (limitTitle == "NA") {
            limitTitle = " "
        } else {
            limitTitle = paste("  Copula = ", limitTitle[1])
        }

        # Plot phi:
        x = (0:N)/N
        Title = paste("Generator Phi - Copula No:", as.character(Copula),
            "\nalpha = ", as.character(alpha), "  Strict = ", strict,
            limitTitle)
        phi.0 = .Phi(x = 0, alpha = alpha, type = as.character(Copula))
        y = .Phi(x = x, alpha = alpha, type = as.character(Copula))
        x = x[y < 1e6]
        y = y[y < 1e6]
        if (is.finite(y[1])) ylim = c(0, y[1]) else ylim = c(0, y[2])
        plot(x = x, y = y, type = "l", ylim = ylim, main = Title[1],
            xlab = "t", ylab = paste("Phi |", RANGE))
        if (N < 100) points(x = x, y = y, pch = 19, cex = 0.5)
        y.inv = .invPhi(x = y, alpha = alpha, type = as.character(Copula))
        lines(x = y.inv, y = y, col = "red", lty = 3)
        abline(h = 0, lty = 3)
        points(0, phi.0, col = "red", pch = 19)

        # Plot phi first and second Derivative:
        y1 = .PhiFirstDer(x = x, alpha = alpha,
            type = as.character(Copula))
        y2 = .PhiSecondDer(x = x, alpha = alpha,
            type = as.character(Copula))
        r1 = max(abs(y1[is.finite(y1)]))
        r2 = max(abs(y2[is.finite(y2)]))
        if (r2 == 0) r2 = 1
        plot(x = x, y = y1/r1, ylim = c(-1, 1), type = "l", xlab = "t",
            ylab = "Derivatives", main = "Phi first and second Derivative",
            col = "blue")
        if (N < 100) points(x = x, y = y1/r1, pch = 19, cex = 0.5)
        lines(x = x, y = y2/r2, col = "red")
        if (N < 100) points(x = x, y = y2/r2, pch = 19, cex = 0.5)
        abline(h = 0, lty = 3)
        mtext("First                  ", 4, col = "blue", cex = 0.75)
        mtext("                 Second", 4, col = "red ", cex = 0.75)
        mtext(paste("x", as.character(round(r1, digits = 2))), 1,
            line = -2, col = "blue", cex = 0.75)
        mtext(paste("x", as.character(round(r2, digits = 2))), 3,
            line = -2, col = "red", cex = 0.75)

        # Plot invPhi:
        Title = paste( "Inverse Phi\n Phi(0) =",
            as.character(round(phi.0, digits = 3)))
        plot(x = y, y = y.inv, type = "l", main = Title,
            xlab = paste("Phi |", RANGE), ylab = "t")
        if (N < 100) points(x = y, y = y.inv, pch = 19, cex = 0.5)
        abline(h = 0, lty = 3)
        points(phi.0, 0, col = "red", pch = 19)

        # Plot invPhi first & second Derivative:
        y = y[y < .Phi0(alpha, Copula)]
        Title = "Inverse Phi 1st Derivative"
        y1.inv = .invPhiFirstDer(x = y, alpha = alpha,
            type = as.character(Copula))
        y2.inv = .invPhiSecondDer(x = y, alpha = alpha,
            type = as.character(Copula))
        r1 = max(abs(y1.inv[is.finite(y1.inv)]))
        r2 = max(abs(y2.inv[is.finite(y2.inv)]))
        if (r2 == 0) r2 = 1
        plot(x = y, y = y1.inv/r1, ylim = c(-1, 1),
            type = "l", xlim = range(y), xlab = paste("Phi |", RANGE),
            ylab = "dewrivatives",
            main = "Inv Phi first and second Derivative", col = "blue")
        if (N < 100) points(x = y, y = y1.inv/r1, pch = 19, cex = 0.5)
        lines(x = y, y = y2.inv/r2, col = "red")
        if (N < 100) points(x = y, y = y2.inv/r2, pch = 19, cex = 0.5)
        abline(h = 0, lty = 3)
        mtext("First                  ", 4, col = "blue", cex = 0.75)
        mtext("                 Second", 4, col = "red ", cex = 0.75)
        mtext(paste("x", as.character(round(r1, digits = 2))), 1,
            line = -2, col = "blue", cex = 0.75)
        mtext(paste("x", as.character(round(r2, digits = 2))), 3,
            line = -2, col = "red", cex = 0.75)

        # Reset Frame:
        par(mfcol = c(2, 2), cex = 0.7)
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 10)
    C1 = "1: [-1,Inf]"
    C2 = "2-4-6-8-12-14-15-21: [1,Inf)"
    C3 = "3: [-1,1)"
    C4 = "5-17: (-Inf,Inf)|{0}"
    C5 = "7-9-10-22: (0,1]"
    C6 = "11: (0, 1/2]"
    C7 = "13-16-19-20: (0,Inf)"
    C8 = "18: [2, Inf)"
    C = c(   C1, C2,   C3,   C4,   C5,   C6,   C7,  C8 )
    L = c(   -1,  1,   -1,   -B,    0,    0,    0,   2 )
    U = c(3*B/5,  B,    1,    B,    1,  0.5,  B/2, 2*B )
    A = c(  0.5,  2,  0.5,    1,  0.5,  0.2,    1,   3 )
    V = rep(0.01, 20)
    .sliderMenu(refresh.code,
        names       = c("Copula",  "N", C),
        minima      = c(       1,   10, L),
        maxima      = c(      22, 1000, U),
        resolutions = c(       1,   10, V),
        starts      = c(       1,  100, A))
}


# ------------------------------------------------------------------------------


.Phi <- 
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Archimedean generator "phi"

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # As listed in Nelsen:
    N = length(x)
    Type = "NA"
    if (type ==  1)
        if (alpha == -1) Type = "W"
        else if (alpha == 0) Type = "Pi"
        else if (alpha == 1) Type = "L"
        else f = 1/alpha*(x^(-alpha)-1)                          # Clayton
    if (type == 2)
        if (alpha == 1) Type = "W"
        else f = (1-x)^alpha
    if (type == 3)
        if (alpha == 0) Type = "Pi"
        else if (alpha == 1) Type = "L"
        else f = log((1-alpha*(1-x))/x)                  # Ali-Mikhail-Haq
    if (type == 4)
        if (alpha == 1) Type = "Pi"
        else f = (-log(x))^alpha                          # Gumbel-Hougard
    if (type == 5)
        if (alpha == 0) Type = "Pi"
        else f = -log((exp(-alpha*x)-1)/(exp(-alpha)-1))           # Frank
    if (type == 6)
        if (alpha == 1) Type = "Pi"
        else f = -log(1-(1-x)^alpha)                                 # Joe
    if (type == 7)
        if (alpha == 0) Type = "W"
        else if (alpha == 1) Type = "Pi"
        else f = -log(alpha*x+(1-alpha))
    if (type == 8)
        if (alpha == 0) Type = "Pi"
        else f = (1-x)/(1+x*(alpha-1))
    if (type == 9)
        if (alpha == 0) Type = "Pi"
        else f = log(1-alpha*log(x))                      # Gumbel-Barnett
    if (type == 10)
        if (alpha == 0) Type = "Pi"
        else f = log(2*x^(-alpha)-1)
    if (type == 11)
        if (alpha == 0) Type = "Pi"
        else f = log(2-x^alpha)
    if (type == 12)
        if (alpha == 1) Type = "L"
        else f = (1/x-1)^alpha
    if (type == 13)
        if (alpha == 1) Type = "Pi"
        else f = (1-log(x))^alpha-1
    if (type == 14)
        if (alpha == 1) Type = "L"
        else f = (x^(-1/alpha)-1)^alpha
    if (type == 15)
        if (alpha == 1) Type = "W"
        else f = (1-x^(1/alpha))^alpha
    if (type == 16)
        if (alpha == 0) Type = "W"
        else f = (alpha/x+1)*(1-x)
    if (type == 17)
        if (alpha == -1) Type = "Pi"
        else f = -log(((1+x)^(-alpha)-1)/(2^(-alpha)-1))
    if (type == 18)
        f = exp(alpha/(x-1))
    if (type == 19)
        if (alpha == 0) Type = "L"
        else f = exp(alpha/x)-exp(alpha)
    if (type == 20)
        if (alpha == 0) Type = "Pi"
        else f = exp(x^(-alpha))-exp(1)
    if (type == 21) if (alpha == 1) Type = "W"
        else f = (1-(1-(1-x)^alpha)^(1/alpha))
    if (type == 22)
        if (alpha == 0) Type = "Pi"
        else f = asin(1-x^alpha)

    if (Type == "Pi") f = -log(x)
    if (Type == "W") f = 1-x
    if (Type == "L") f = 1/x - 1

    f[x == 0] = .Phi0(alpha, type)

    # Return Value:
    f
}


# ------------------------------------------------------------------------------


.Phi0 <-
    function(alpha, type)
{   
    # A function implemented by Diethelm Wuertz

    # Phi(0):
    type <- as.integer(type)
    if (type == 1) phi0 = if (alpha < 0) -1/alpha else Inf
    else if (type ==  2) phi0 = 1
    else if (type ==  3) phi0 = Inf
    else if (type ==  4) phi0 = Inf
    else if (type ==  5) phi0 = Inf
    else if (type ==  6) phi0 = Inf
    else if (type ==  7) phi0 = if (alpha == 0) 1 else -log(1 - alpha)
    else if (type ==  8) phi0 = 1
    else if (type ==  9) phi0 = Inf
    else if (type == 10) phi0 = Inf
    else if (type == 11) phi0 = if (alpha == 0) Inf else log(2)
    else if (type == 12) phi0 = Inf
    else if (type == 13) phi0 = Inf
    else if (type == 14) phi0 = Inf
    else if (type == 15) phi0 = 1
    else if (type == 16) phi0 = if (alpha == 0) 1 else Inf
    else if (type == 17) phi0 = Inf
    else if (type == 18) phi0 = exp(-alpha)
    else if (type == 19) phi0 = Inf
    else if (type == 20) phi0 = Inf
    else if (type == 21) phi0 = 1
    else if (type == 22) phi0 = if (alpha == 0) Inf else pi/2

    # Return Value:
    phi0
}


# ------------------------------------------------------------------------------


.PhiFirstDer <-
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Derivative of Archimedean generator.

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # FUNCTION:

    # The functions were created by MAPLE:
    N = length(x)
    cType = "NA"
    if (Type == 1)
        if (alpha == -1) cType = "W"
        else if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else f1 = -x^(-alpha-1)
    if (Type == 2)
        if (alpha == 1) cType = "W"
        else f1 = -(1-x)^alpha*alpha/(1-x)
    if (Type == 3)
        if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else f1 = (alpha/x-(1-alpha*(1-x))/x^2)/(1-alpha*(1-x))*x
    if (Type == 4)
        if (alpha == 1) cType = "Pi"
        else f1 = (-log(x))^alpha*alpha/x/log(x)
    if (Type == 5)
        if (alpha == 0) cType = "Pi"
        else f1 = alpha*exp(-alpha*x)/(exp(-alpha*x)-1)
    if (Type == 6)
        if (alpha == 1) cType = "Pi"
        else f1 = -(1-x)^alpha*alpha/(1-x)/(1-(1-x)^alpha)
    if (Type == 7)
        if (alpha == 0) cType = "W"
        else if (alpha == 1) cType = "Pi"
        else f1 = -alpha/(alpha*x+1-alpha)
    if (Type == 8)
        if (alpha == 1) cType = "W"
        else f1 = -1/(1+x*(-1+alpha))-(1-x)/(1+x*(-1+alpha))^2*(-1+alpha)
    if (Type == 9)
        if (alpha == 0) cType = "Pi"
        else f1 = -alpha/x/(1-alpha*log(x))
    if (Type == 10)
        if (alpha == 0) cType = "Pi"
        else f1 = -2*x^(-alpha)*alpha/x/(2*x^(-alpha)-1)
    if (Type == 11)
        if (alpha == 0) cType = "Pi"
        else f1 = -x^alpha*alpha/x/(2-x^alpha)
    if (Type == 12)
        if (alpha == 1) cType = "L"
        else f1 = -(1/x-1)^alpha*alpha/x^2/(1/x-1)
    if (Type == 13)
        if (alpha == 1) cType = "Pi"
        else f1 = -(1-log(x))^alpha*alpha/x/(1-log(x))
    if (Type == 14)
        if (alpha == 1) cType = "L"
        else f1 = -(x^(-1/alpha)-1)^alpha*x^(-1/alpha)/x/(x^(-1/alpha)-1)
    if (Type == 15)
        if (alpha == 1) cType = "W"
        else f1 = -(1-x^(1/alpha))^alpha*x^(1/alpha)/x/(1-x^(1/alpha))
    if (Type == 16)
        if (alpha == 0) cType = "W"
        else f1 = -alpha/x^2*(1-x)-alpha/x-1
    if (Type == 17)
        if (alpha == -1) cType = "Pi"
        else f1 = (1+x)^(-alpha)*alpha/(1+x)/((1+x)^(-alpha)-1)
    if (Type == 18)
        f1 = -alpha/(-1+x)^2*exp(alpha/(-1+x))
    if (Type == 19)
        if (alpha == 0) cType = "L"
        else f1 = -alpha/x^2*exp(alpha/x)
    if (Type == 20)
        if (alpha == 0) cType = "Pi"
        else f1 = -x^(-alpha)*alpha/x*exp(x^(-alpha))
    if (Type == 21)
        if (alpha == 1) cType = "W"
        else f1 = -(1-(1-x)^alpha)^(-(-1+alpha)/alpha)*(1-x)^(-1+alpha)
    if (Type == 22)
        if (alpha == 0) cType = "Pi"
        else f1 = -x^(-1+alpha)*alpha/(2*x^alpha-x^(2*alpha))^(1/2)

    if (cType == "Pi") f1 = -1/x
    if (cType == "W") f1 = rep(-1, times = N)
    if (cType == "L") f1 = -1/x^2

    # Return Value:
    f1
}


# ------------------------------------------------------------------------------


.PhiSecondDer <- 
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Derivative of Archimedean generator.

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # FUNCTION:

    # The functions were created by MAPLE:
    a = alpha
    N = length(x)
    cType = "NA"
    if (Type == 1)
        if (alpha == -1) cType = "W"
        else if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else f2 = x^(-a-2)*a+x^(-a-2)
    if (Type == 2)
        if (alpha == 1) cType = "W"
        else f2 = (1-x)^(a-2)*a^2-(1-x)^(a-2)*a
    if (Type == 3)
        if (alpha == 0) cType = "Pi"
        else if (alpha == 1) Type = "L"
        else f2 = -1/x^2*(a-1)*(1-a+2*x)/(1-a+x)^2
    if (Type == 4)
        if (alpha == 1) cType = "Pi"
        else f2 = a*((-log(x))^(a-2)*a+(-log(x))^(a-1)-(-log(x))^(a-2))/x^2
    if (Type == 5)
        if (alpha == 0) cType = "Pi"
        else f2 = a^2*exp(-a*x)/(exp(-a*x)-1)^2
    if (Type == 6)
        if (alpha == 1) cType = "Pi"
        else f2 = a*((1-x)^(a-2)*a-(1-x)^(a-2)+(1-x)^(2*a-2))/(-1+(1-x)^a)^2
    if (Type == 7)
        if (alpha == 0) cType = "W"
        else if (alpha == 1) cType = "Pi"
        else f2 = alpha^2/(alpha*x+1-alpha)^2
    if (Type == 8)
        if (alpha == 1) cType = "W"
        else f2 = 2*(a-1)*a/(1+a*x-x)^3
    if (Type == 9)
        if (alpha == 0) cType = "Pi"
        else f2 = -a*(-1+a*log(x)+a)/x^2/(-1+a*log(x))^2
    if (Type == 10)
        if (alpha == 0) cType = "Pi"
        else f2 = -2*a*(x^a*a-2+x^a)/(-2+x^a)^2/x^2
    if (Type == 11)
        if (alpha == 0) cType = "Pi"
        else f2 = -a*(2*x^(a-2)*a-2*x^(a-2)+x^(2*a-2))/(-2+x^a)^2
    if (Type == 12)
        if (alpha == 1) cType = "L"
        else f2 = -(-(x-1)/x)^a*a*(-a+2*x-1)/x^2/(x-1)^2
    if (Type == 13)
        if (alpha == 1) cType = "Pi"
        else f2 = a*((1-log(x))^(a-2)*a+(1-log(x))^(a-1)-(1-log(x))^(a-2))/x^2
    if (Type == 14)
        if (alpha == 1) cType = "L"
        else f2 = ((x^(-1/a)-1)^(a-2)*x^(-2*(a+1)/a)*a+(x^(-1/a)-1)^(a-1) *
            x^(-(1+2*a)/a)+(x^(-1/a)-1)^(a-1)*x^(-(1+2*a)/a) *
            a-(x^(-1/a)-1)^(a-2)*x^(-2*(a+1)/a))/a
    if (Type == 15)
        if (alpha == 1) cType = "W"
        else f2 = ((1-x^(1/a))^(a-2)*x^(-2*(a-1)/a)*a-(1-x^(1/a))^(a-1) *
            x^(-(-1+2*a)/a)+(1-x^(1/a))^(a-1)*x^(-(-1+2*a)/a) *
            a-(1-x^(1/a))^(a-2)*x^(-2*(a-1)/a))/a
    if (Type == 16)
        if (alpha == 0) cType = "W"
        else f2 = 2*a/x^3
    if (Type == 17)
        if (alpha == -1) cType = "Pi"
        else f2 = a*((1+x)^(a-2)*a+2*(1+x)^(a-2)*a*x+(1+x)^(a-2)*a*x^2 -
            1+(1+x)^(a-2)+2*(1+x)^(a-2)*x+(1+x)^(a-2)*x^2) /
            (-1+(1+x)^a)^2/(1+x)^2
    if (Type == 18)
        f2 = a*exp(a/(x-1))*(2*x-2+a)/(x-1)^4
    if (Type == 19)
        if (alpha == 0) cType = "L"
        else f2 = a*exp(a/x)*(2*x+a)/x^4
    if (Type == 20)
        if (alpha == 0) cType = "Pi"
        else f2 = a*exp(x^(-a))*(x^(-a-2)*a+x^(-a-2)+x^(-2*a-2)*a)
    if (Type == 21)
        if (alpha == 1) cType = "W"
        else f2 = -(1-(1-x)^a)^(-(-1+2*a)/a)*(1-x)^(2*a-2) +
            (1-(1-x)^a)^(-(-1+a)/a)*(1-x)^(a-2)*a -
            (1-(1-x)^a)^(-(-1+a)/a)*(1-x)^(a-2) +
            (1-(1-x)^a)^(-(-1+2*a)/a)*(1-x)^(2*a-2)*a
    if (Type == 22)
        if (alpha == 0) cType = "Pi"
        else f2 = -a/x^2*(a*x^(2*a)-2*x^(2*a)+x^(3*a))/(2*x^a-x^(2*a))^(3/2)

    if (cType == "Pi") f2 = 1/x^2
    if (cType == "W") f2 = rep(0, times = N)
    if (cType == "L") f2 = 2/x^3

    # Return Value:
    f2
}


# ------------------------------------------------------------------------------


.invPhi <- 
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes inverse of Archimedean generator.

    # FUNCTION:

    # Type:
    type <- match.arg(type)
    Type <- as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check <- archmCheck(alpha, type)

    # Inverse Generator:
    N = length(x)
    cType = "NA"
    if (Type == 1)
        if (alpha == -1) cType = "W"
        else if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else finv = exp(-log(1 + alpha*x)/alpha)
    if (Type == 2)
        if (alpha == 1) cType = "W"
        else finv = 1 - x^(1/alpha)
    if (Type == 3)
        if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else finv = (1-alpha) / (exp(x)-alpha)
    if (Type == 4)
        if (alpha == 1) cType = "Pi"
        else finv = exp(-x^(1/alpha))
    if (Type == 5)
        if (alpha == 0) cType = "Pi"
        else finv = -log(1+exp(-x)*( exp(-alpha)-1 ) ) / alpha
    if (Type == 6)
        if (alpha == 1) cType = "Pi"
        else finv = 1 - (1 - exp(-x))^(1/alpha)
    if (Type == 7)
        if (alpha == 0) cType = "W"
        else if (alpha == 1) Type = "Pi"
        else finv = (1-exp(x)+alpha*exp(x))/alpha/exp(x)
    if (Type == 8)
        if (alpha == 1) cType = "W"
        else finv = (1-x) / ((alpha-1)*x+1)
    if (Type == 9)
        if (alpha == 0) cType = "Pi"
        else  finv = exp((1-exp(x))/alpha)
    if (Type == 10)
        if (alpha == 0) cType = "Pi"
        else finv = ((1+exp(x))/2 )^(-1/alpha)
    if (Type == 11)
        if (alpha == 0) cType = "Pi"
        else finv = (2-exp(x))^(1/alpha)
    if (Type == 12)
        if (alpha == 1) cType = "L"
        else finv = 1/(1+x^(1/alpha))
    if (Type == 13)
        if (alpha == 1) cType = "Pi"
        else finv = exp(1-(1+x)^(1/alpha))
    if (Type == 14)
        if (alpha == 1) cType = "L"
        else finv = (1+x^(1/alpha))^(-alpha)
    if (Type == 15)
        if (alpha == 1) cType = "W"
        else finv = (1-x^(1/alpha))^alpha
    if (Type == 16)
        if (alpha == 0) cType = "W"
        else finv = (1-alpha-x)/2 + sqrt(((1-alpha-x)^2)/4+alpha)
    if (Type == 17)
        if (alpha == -1) cType = "Pi"
        else finv = (exp(-x)*(2^(-alpha)-1)+1)^(-1/alpha) - 1
    if (Type == 18)
        finv = 1+alpha/log(x)
    if (Type == 19)
        if (alpha == 0) cType = "L"
        else finv = alpha / log(x+exp(alpha))
    if (Type == 20)
        if (alpha == 0) cType = "Pi"
        else finv = exp( -log((log(x+exp(1))))/alpha)
    if (Type == 21)
        if (alpha == 1) cType = "W"
        else finv = 1-(1-(1-x)^alpha)^(1/alpha)
    if (Type == 22)
        if (alpha == 0) cType = "Pi"
        else finv = (1-sin(x))^(1/alpha)

    if (cType == "Pi") finv = exp(-x)
    if (cType == "W") finv = 1 - x
    if (cType == "L") finv = 1 / (x+1)

    # Large x Limit:
    finv[which(x >= .Phi0(alpha, type))] = 0

    # Return Value:
    finv
}


# ------------------------------------------------------------------------------


.invPhiFirstDer <- 
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes first Derivative of inverse Archimedean generator.

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # Generator:
    N = length(x)
    cType = "NA"
    a = alpha
    y = x
    ln = log
    if (Type == 1)
        if (alpha == -1) cType = "W"
        else if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else finv1 = -(1+y*a)^(-(a+1)/a)
    if (Type == 2)
        if (alpha == 1) cType = "W"
        else finv1 = -y^(-(a-1)/a)/a
    if (Type == 3)
        if (alpha == 0) cType = "Pi"
        else if (alpha == 1) cType = "L"
        else  finv1 = (a-1)/(exp(y)-1)^2*exp(y)
    if (Type == 4)
        if (alpha == 1) cType = "Pi"
        else finv1 = -y^(-(a-1)/a)/a*exp(-y^(1/a))
    if (Type ==  5)
        if (alpha == 0) cType = "Pi"
        else finv1 = (-1+exp(a))/(-1+exp(a)-exp(y+a))/a
    if (Type == 6)
        if (alpha == 1) cType = "Pi"
        else finv1 = -exp(-(-ln(exp(y)-1)+y)/a)/(exp(y)-1)/a
    if (Type == 7)
        if (alpha == 0) cType = "W"
        else if (alpha == 1) Type = "Pi"
        else finv1 = (-exp(y)+a*exp(y))/a/exp(y)-(1-exp(y)+a*exp(y))/a/exp(y)
    if (Type == 8)
        if (alpha == 1) cType = "W"
        else finv1 = -a/(1+y*a-y)^2
    if (Type == 9)
        if (alpha == 0) cType = "Pi"
        else finv1 = -1/a*exp((y*a-exp(y)+1)/a)
    if (Type == 10)
        if (alpha == 0) cType = "Pi"
        else finv1 = -1/(exp(y)+1)/a*exp((y*a+ln(2)-ln(exp(y)+1))/a)
    if (Type == 11)
        if (alpha == 0) cType = "Pi"
        else finv1 = -(-exp(y)+2)^(-(a-1)/a)/a*exp(y)
    if (Type == 12)
        if (alpha == 1) cType = "L"
        else finv1 = -1/(y^(1/a)+1)^2*y^(-(a-1)/a)/a
    if (Type == 13)
        if (alpha == 1) cType = "Pi"
        else finv1 = -(1+y)^(-(a-1)/a)/a*exp(-(1+y)^(1/a)+1)
    if (Type == 14)
        if (alpha == 1) cType = "L"
        else finv1 = -(y^(1/a)+1)^(-a-1)*y^(-(a-1)/a)
    if (Type == 15)
        if (alpha == 1) cType = "L"
        else finv1 = -(-y^(1/a)+1)^(a-1)*y^(-(a-1)/a)
    if (Type == 16)
        if (alpha == 0) cType = "W"
        else finv1 = -1/2+1/4/(a^2+2*a+2*a*y+1-2*y+y^2)^(1/2)*(2*a-2+2*y)
    if (Type == 17)
        if (alpha == -1) cType = "Pi"
        else finv1 = -(2^(-a)-1+exp(y))^(-1/a)*exp(1/a*y) *
            (-1+2^a)/a/(1-2^a+exp(y)*2^a)
    if (Type == 18)
        finv1 = -a/ln(y)^2/y
    if (Type == 19)
        if (alpha == 0) cType = "L"
        else finv1 = -a/ln(exp(a)+y)^2/(exp(a)+y)
    if (Type == 20)
        if (alpha == 0) cType = "Pi"
        else finv1 = -ln(exp(1)+y)^(-(a+1)/a)/a/(exp(1)+y)
    if (Type == 21)
        if (alpha == 1) cType = "W"
        else finv1 = -exp((log(1-y)*a^2+log(-(1-y)^a+1))/a)/(-1+y)/((1-y)^a-1)
    if (Type == 22)
        if (alpha == 0) cType = "Pi"
        else finv1 = -cos(y)*(1-sin(y))^(-(-1+a)/a)/a

    if (cType == "Pi") finv1 = -exp(-x)
    if (cType == "W") finv1 = rep(-1, times = N)
    if (cType == "L") finv1 = -1 / (x+1)^2

    # Large x Limit:
    finv1[which(x >= .Phi0(a, type))] = 0

    # Return Value:
    finv1
}


# ------------------------------------------------------------------------------


.invPhiSecondDer <- 
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes first Derivative of inverse Archimedean generator.

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # Generator:
    N = length(x)
    cType = "NA"
    a = alpha
    y = x
    ln = log
    if (Type == 1) if (alpha == 0) finv2 = exp(-y) else finv2 =
        finv2 = (1+y*a)^(-(2*a+1)/a)*(a+1)
    if (Type == 2)
        if (alpha == 1) cType = "W"
        else finv2 = y^(-(2*a-1)/a)*(a-1)/a^2
    if (Type == 3)
        if (alpha == 0) cType = "Pi"
        else if (alpha == 1) Type = "L"
        else  finv2 = -(a-1)*exp(y)*(exp(y)+1)/(exp(y)-1)^3
    if (Type == 4)
        if (alpha == 1) cType = "Pi"
        else finv2 = exp(-y^(1/a))*(y^(-(2*a-1)/a)*a-y^(-(2*a-1)/a) +
            y^(-2*(a-1)/a))/a^2
    if (Type == 5)
        if (alpha == 0) cType = "Pi"
        else finv2 = (-1+exp(a))/(-1+exp(a)-exp(y+a))^2/a*exp(y+a)
    if (Type == 6)
        if (alpha == 1) cType = "Pi"
        else finv2 = (-exp(-(-ln(exp(y)-1)+y)/a) +
            exp((ln(exp(y)-1)-y+y*a)/a)*a) / (exp(y)-1)^2/a^2
    if (Type == 7)
        if (alpha == 0) cType = "W"
        else if (alpha == 1) Type = "Pi"
        else finv2 = -(-exp(y)+a*exp(y))/a/exp(y)+(1-exp(y)+a*exp(y))/a/exp(y)
    if (Type == 8)
        if (alpha == 1) cType = "W"
        else finv2 = 2*a/(1+y*a-y)^3*(a-1)
    if (Type == 9)
        if (alpha == 0) cType = "Pi"
        else finv2 = -1/a^2*(a-exp(y))*exp((y*a-exp(y)+1)/a)
    if (Type == 10)
        if (alpha == 0) cType = "Pi"
        else finv2 = -(exp((y*a+ln(2)-ln(exp(y)+1))/a)*a-exp((2*y*a+ln(2) -
            ln(exp(y)+1))/a))/(exp(y)+1)^2/a^2
    if (Type == 11)
        if (alpha == 0) cType = "Pi"
        else finv2 = -exp(y)*((-exp(y)+2)^(-(2*a-1)/a)*exp(y)*a -
            (-exp(y)+2)^(-(2*a-1)/a)*exp(y)+(-exp(y)+2)^(-(a-1)/a)*a)/a^2
    if (Type == 12)
        if (alpha == 1) cType = "L"
        else finv2 = (y^(-2*(a-1)/a)+y^(-2*(a-1)/a)*a+y^(-(2*a-1)/a)*a -
            y^(-(2*a-1)/a))/(y^(1/a)+1)^3/a^2
    if (Type == 13)
        if (alpha == 1) cType = "Pi"
        else finv2 = exp(-(1+y)^(1/a)+1)*((1+y)^(1/a)*a-(1+y)^(1/a) +
            (1+y)^(-2*(a-1)/a)+2*(1+y)^(-2*(a-1)/a)*y +
            (1+y)^(-2*(a-1)/a)*y^2)/a^2/(1+2*y+y^2)
    if (Type == 14)
        if (alpha == 1) cType = "L"
        else finv2 = ((y^(1/a)+1)^(-a-2)*y^(-2*(a-1)/a)*a +
            (y^(1/a)+1)^(-a-2)*y^(-2*(a-1)/a)+(y^(1/a)+1)^(-a-1) *
            y^(-(2*a-1)/a)*a-(y^(1/a)+1)^(-a-1)*y^(-(2*a-1)/a))/a
    if (Type == 15)
        if (alpha == 1) cType = "L"
        else finv2 = (a-1)*((-y^(1/a)+1)^(a-2)*y^(-2*(a-1)/a) +
            (-y^(1/a)+1)^(a-1)*y^(-(2*a-1)/a))/a
    if (Type == 16)
        if (alpha == 0) cType = "W"
        else finv2 = 2*a/(a^2+2*a+2*a*y+1-2*y+y^2)^(3/2)
    if (Type == 17)
        if (alpha == -1) cType = "Pi"
        else finv2 = (2^(-a)-1+exp(y))^(-1/a)*(exp(y*(a+1)/a) -
            2^(a+1)*exp(y*(a+1)/a)+exp(y*(a+1)/a)*4^a +
            exp(1/a*y)*2^(-a)-3*exp(1/a*y)+3*exp(1/a*y)*2^a -
            exp(1/a*y)*4^a-exp(y*(a+1)/a)*a+2^(a+1) *
            exp(y*(a+1)/a)*a- exp(y*(2*a+1)/a)*a*2^a -
            exp(y*(a+1)/a)*a*4^a+exp(y*(2*a+1)/a)*a*4^a)/a^2/(2^(-a)-1 +
            exp(y))/(1-2^a+exp(y)*2^a)^2
    if (Type == 18)
        finv2 = a*(2+ln(y))/ln(y)^3/y^2
    if (Type == 19)
        if (alpha == 0) cType = "L"
        else finv2 = a*(2+ln(exp(a)+y))/ln(exp(a)+y)^3/(exp(a)+y)^2
    if (type == 20)
        if (alpha == 0) cType = "Pi"
        else finv2 = (ln(exp(1)+y)^(-(2*a+1)/a)*a +
            ln(exp(1)+y)^(-(2*a+1)/a) +
            ln(exp(1)+y)^(-(a+1)/a)*a)/a^2/(exp(1)+y)^2
    if (Type == 21)
        if (alpha == 1) cType = "W"
        else finv2 = -(-(1-y)^a+1)^(1/a)*((1-y)^(2*a)-(1-y)^a -
            a*(1-y)^(2*a)+a*(1-y)^a+(1-y)^(2*a-2)*a -
            2*(1-y)^(2*a-2)*a*y+(1-y)^(2*a-2)*a*y^2 -(
            1-y)^(2*a-2)+2*(1-y)^(2*a-2)*y-(1-y)^(2*a-2)*y^2) /
            (-1+y)^2/(-(1-y)^(2*a)+2*(1-y)^a-1)
    if (Type == 22)
        if (alpha == 0) cType = "Pi"
        else finv2 = -(1-sin(y))^(1/a)*(cos(y)^2 +
            a*sin(y)-2*sin(y)+a-2)/cos(y)^2/a^2

    if (cType == "Pi") finv2 = exp(-x)
    if (cType == "W") finv2 = rep(0, times = N)
    if (cType == "L") finv2 = 2 / (x+1)^3

    # Large x Limit:
    finv2[which(x>=.Phi0(a, type))] = 0

    # Return Value:
    finv2
}


################################################################################
# FUNCTION:                  ARCHIMEDEAN DENSITY K GENERATOR:
#  Kfunc                      Computes Archimedean Density Kc and its Inverse
#  KfuncSlider                Displays interactively the density and concordance
#  .Kfunc                     Computes Density for Archimedean Copulae
#  .invK                      Computes Inverse of Density
#  .invK2                     Utility Function
#  .ALPHA                     Utility Function
#  .TAU                       Utility Function
#  .RHO                       Utility Function


Kfunc <-
    function(x, alpha = NULL, type = archmList(), inv = FALSE, lower = 1.0e-8)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes density and its inverse for Archimedean Copulae

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Default alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # Density or its inverse:
    if (!inv) {
        ans = .Kfunc(x, alpha, type)
        names(ans)<-"Kfunc"
    } else {
        ans = .invK(x, alpha, type, lower)
        names(ans)<-"invK"
    }

    # Add Control Attribute:
    attr(ans, "control")<-cbind.data.frame(alpha = alpha, type = type,
        inv = inv, lower = lower, row.names = "")

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


KfuncSlider <- 
    function(B = 5)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays interactively the density and concordance

    # FUNCTION:

    # Graphic Frame:
    par(mfcol = c(2, 2), cex = 0.7)

    # Internal Function:
    refresh.code = function(...)
    {
        # Startup Counter:
        .counter <- getRmetricsOptions(".counter") + 1
        setRmetricsOptions(.counter = .counter)
        if (.counter < 10) return ()

        # Sliders:
        Copula = as.integer(.sliderMenu(no = 1))
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)

        # Skip:
        if (Copula == 13 & alpha == 0) return(invisible())

        # Do we have a strict Copula?
        strict = c(
            "Yes","No","Yes","Yes","Yes","Yes","No","No","Yes","Yes",
            "No","Yes","Yes","Yes","No","Yes","Yes","No","Yes","Yes",
            "No","Yes")[Copula]
        if (alpha < 0 & Copula == 1) strict[1] = "No"
        if (alpha == 0 & Copula == 16) strict[16] = "No"

        # What is the Range?
        RANGE = c(
            "-1|Inf", "1|Inf", "-1|1", "-Inf|inf", "0|1", "0|0.5",
            "0|Inf", "2|Inf")[No]

        # Which one is the Limit Copula?
        limitTitle = rep("NA", times = 22)
        if (alpha == -1)
            limitTitle = c(
                "W ", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "NA", "NA", "NA", "NA", "NA", "NA", "Pi", "NA", "NA", "NA",
                "NA", "NA")
        if (alpha == 0)
            limitTitle = c(
                "Pi", "NA", "Pi", "NA", "Pi", "NA", "W ", "NA", "Pi", "Pi",
                "Pi", "NA", "NA", "NA", "NA", "W ", "NA", "NA", "L ", "Pi",
                "NA", "Pi")
        if (alpha == 1)
            limitTitle = c(
                "L ", "W ", "L ", "Pi", "NA", "Pi", "Pi", "W ", "NA", "NA",
                "NA", "L ", "Pi", "L ", "W ", "NA", "NA", "NA", "NA", "NA",
                "W ", "NA")
        limitTitle = limitTitle[Copula]
        if (limitTitle == "NA") {
            limitTitle = " "
        } else {
            limitTitle = paste("  Copula = ", limitTitle[1])
        }

        # Plot 1 - Kfunc:
        x = (0:N)/N
        y = .Kfunc(x = x, alpha = alpha, type = as.character(Copula))
        plot(x = x, y = y, ylim = c(0, 1), type = "l", xlab = "t", ylab = "K")
        title(main = paste("K - Archimedean Copula No:", as.character(Copula),
            "\nalpha = ", as.character(alpha), "  Strict = ", strict,
            limitTitle))
        if (N < 100) points(x = x, y = y, pch = 19, cex = 0.5)
        y10 = .Kfunc(x = (0:10)/10, alpha = alpha, type = as.character(Copula))
        invK10 = .invK2(y10, alpha = alpha, type = as.character(Copula))
        points(invK10, y10, col = "red")
        text(x = 0.8, y = 0.075, labels = "Test: invK[invK]", col = "red")

        # Plot 2 - archmTau:
        tau = .archmTau(alpha = alpha, type = as.character(Copula))
        rho = approx(.ALPHA[, Copula], .RHO[, Copula], xout = alpha)$y
        plot(x = .ALPHA[, Copula], y = .TAU[, Copula], ylim = c(-1, 1),
            type = "l", col = "red",
            xlab = paste("alpha: ", RANGE, sep = ""), ylab = "Tau")
        # points(x = .ALPHA[, Copula], y = .TAU[, Copula], pch = 19, cex = 0.5)
        lines(x = .ALPHA[, Copula], y = .RHO[, Copula], col = "blue")
        # points(x = .ALPHA[, Copula], y = .RHO[, Copula], pch = 19, cex = 0.5)

        points(x = alpha, y = tau, pch = 19, col = "red")
        abline(h = .archmTauRange(type = as.character(Copula))[1], lty =3,
            col = "steelblue")
        abline(h = .archmTauRange(type = as.character(Copula))[2], lty =3,
            col = "steelblue")
        points(x = alpha, y = rho, col = "blue", pch = 19)
        mtext("rho          ", 4, col = "blue", cex = 0.75)
        mtext("          tau", 4, col = "red ", cex = 0.75)
        title(main = paste("Concordance Measures",
            "\ntau = ", as.character(round(tau, digits = 2)),
            "rho = ", as.character(round(rho, digits = 2)) ) )


        plot(x = y, y = x, xlim = c(0, 1), type = "l", xlab = "K", ylab = "t")
        title(main = "Inverse K")

        # Plot 3 - lambda U:
        # xTail = 1 - (1/2)^(1:20)
        # Tail = .archmTail(alpha = alpha, type = as.character(Copula))
        # plot(x = xTail, y = Tail$lambdaU.Cuv, col = "blue",
        #     xlim = c(0, 1), ylim = c(0, 1), main = "Tail Dependence")
        # points(x = xTail, y = Tail$lambdaU.Phi, col = "red", pch = 3)

        # Rho:
        # Rho = NULL
        # for ( a in Alpha)
        #   Rho = c(Rho, archmRho(alpha = a, type = as.character(Copula)))
        # lines(x = Alpha, y = Rho, type = "l", col = "blue")
        # rho =  archmRho(alpha = alpha, type = as.character(Copula))
        # points(x = alpha, y = rho, col = "red", pty = 19)

        # plot(rnorm(100))
        # plot(rnorm(100))

        # Reset Frame:
        par(mfcol = c(2, 2), cex = 0.7)
    }

    # Open Slider Menu:
    setRmetricsOptions(.counter = 0)
    C1 = "1: [-1,Inf]"
    C2 = "2-4-6-8-12-14-15-21: [1,Inf)"
    C3 = "3: [-1,1)"
    C4 = "5-17: (-Inf,Inf)|{0}"
    C5 = "7-9-10-22: (0,1]"
    C6 = "11: (0, 1/2]"
    C7 = "13-16-19-20: (0,Inf)"
    C8 = "18: [2, Inf)"
    C = c(   C1,  C2,   C3,  C4,   C5,   C6,   C7,  C8 )
    L = c(   -1,   1,   -1,  -B,    0,    0,    0,   2 )
    U = c(    B, 5*B,    1, 5*B,    1,  0.5,    B,   B )
    A = c(  0.5,   2,  0.5,   1,  0.5,  0.2,    1,   3 )
    V = rep(0.01, 20)
    .sliderMenu(refresh.code,
        names       = c("Copula",  "N", C),
        minima      = c(       1,   10, L),
        maxima      = c(      22, 1000, U),
        resolutions = c(       1,   10, V),
        starts      = c(       1,  100, A))
}


# ------------------------------------------------------------------------------


.Kfunc <- 
    function(x, alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Density for Archimedean Copulae

    # Arguments:
    #   x - a numeric vector

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Missing x:
    if (missing(x)) x = (0:10)/10

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # Density:
    Kfunc = x - .Phi(x, alpha, type) / .PhiFirstDer(x, alpha, type)

    # Take care from divergencies:
    Kfunc[is.na(Kfunc)] = 0
    Kfunc[x == 1] = 1

    # Return Value:
    Kfunc
}


# ------------------------------------------------------------------------------


.invK <- 
    function(x, alpha = NULL, type = archmList(), lower = 1.0e-8)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Inverse of Density for Archimedean Copulae

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param

    # Check alpha:
    check = archmCheck(alpha, type)

    # Compute Inverse:
    .fKC = function(x, p, alpha, type) { .Kfunc (x, alpha, type) - p }
    p = x
    z = NULL
    for (P in p) {
        if (P > 1 - lower/2) {
            res = 1
        } else if (P < .Kfunc(0, alpha, type) + lower/2 ) {
            res = 0
        } else {
            res = uniroot(.fKC, c(lower, 1),
                p = P, alpha = alpha, type = type)$root
        }
        z = c(z, res)
    }

    # Return Value:
    z
}


# ------------------------------------------------------------------------------


.invK2 <- 
    function(x, alpha, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes from tabulated values

    # FUNCTION:

    # Type:
    type = match.arg(type)
    Type = as.integer(type)

    # Tabulated Values:
    iK = NULL
    for (i in 1:length(x)) {
        Ord = order(abs(.Kfunc((0:1000)/1000, alpha, type)-x[i]))[1]/1000
        iK = c(iK, Ord)
    }

    # Return Value:
    iK
}


# ------------------------------------------------------------------------------


.makeConcordanceTable <- 
    function(B = 5, dump = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Make Table:
    Counter <- c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
    L = c( -1,  +1, -1, -5*B, 0,   0, 0, 2 )
    U = c(  B, 5*B,  1,  5*B, 1, 0.5, B, B )
    Tau = Alpha = Rho = NULL
    for (i in 1:22) {
        print(i)
        No = Counter[i]
        lower = L[No]
        upper = U[No]
        alpha = seq(lower, upper, length = 25)
        Alpha = cbind(Alpha, alpha)
        tau = archmTau(alpha = alpha, type = i)
        rho = archmRho(alpha = alpha, type = i)
        Tau = cbind(Tau, tau)
        Rho = cbind(Rho, rho)
    }
    .ALPHA = data.frame(Alpha)
    .TAU = data.frame(Tau)
    .RHO = data.frame(Rho)
    colnames(.ALPHA) = colnames(.TAU) = colnames(.RHO) = as.character(1:22)

    # Dump:
    if (dump) {
         dump(".ALPHA", "alpha.R")
         dump(".TAU", "tau.R")
         dump(".RHO", "rho.R")
    }

    # Return Value:
    list(ALPHA = .ALPHA, TAU = .TAU, RHO = .RHO)
}


# ------------------------------------------------------------------------------


".ALPHA" <-
    structure(list(

    "1" =
    c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75,
    2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5),

    "2" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25),

    "3" =
    c(-1, -0.916666666666667, -0.833333333333333, -0.75,
    -0.666666666666667, -0.583333333333333, -0.5, -0.416666666666667,
    -0.333333333333333, -0.25, -0.166666666666667, -0.0833333333333334,
    0, 0.0833333333333333, 0.166666666666667, 0.25, 0.333333333333333,
    0.416666666666667, 0.5, 0.583333333333333, 0.666666666666667,
    0.75, 0.833333333333333, 0.916666666666667, 1),

    "4" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25),

    "5" =
    c(-25, -22.9166666666667, -20.8333333333333, -18.75, -16.6666666666667,
    -14.5833333333333, -12.5, -10.4166666666667, -8.33333333333333,
    -6.25, -4.16666666666666, -2.08333333333333, 0, 2.08333333333334,
    4.16666666666667, 6.25, 8.33333333333334, 10.4166666666667, 12.5,
    14.5833333333333, 16.6666666666667, 18.75, 20.8333333333333,
    22.9166666666667, 25),

    "6" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25),

    "7" =
    c(0, 0.0416666666666667, 0.0833333333333333, 0.125, 0.166666666666667,
    0.208333333333333, 0.25, 0.291666666666667, 0.333333333333333, 0.375,
    0.416666666666667, 0.458333333333333, 0.5, 0.541666666666667,
    0.583333333333333, 0.625, 0.666666666666667, 0.708333333333333,
    0.75, 0.791666666666667, 0.833333333333333, 0.875, 0.916666666666667,
    0.958333333333333, 1),

    "8" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25),

    "9" =
    c(0, 0.0416666666666667, 0.0833333333333333, 0.125, 0.166666666666667,
    0.208333333333333, 0.25, 0.291666666666667, 0.333333333333333, 0.375,
    0.416666666666667, 0.458333333333333, 0.5, 0.541666666666667,
    0.583333333333333, 0.625, 0.666666666666667, 0.708333333333333,
    0.75, 0.791666666666667, 0.833333333333333, 0.875, 0.916666666666667,
    0.958333333333333, 1),

    "10" =
    c(0, 0.0416666666666667, 0.0833333333333333, 0.125, 0.166666666666667,
    0.208333333333333, 0.25, 0.291666666666667, 0.333333333333333, 0.375,
    0.416666666666667, 0.458333333333333, 0.5, 0.541666666666667,
    0.583333333333333, 0.625, 0.666666666666667, 0.708333333333333, 0.75,
    0.791666666666667, 0.833333333333333, 0.875, 0.916666666666667,
    0.958333333333333, 1),

    "11" =
    c(0, 0.0208333333333333, 0.0416666666666667, 0.0625, 0.0833333333333333,
    0.104166666666667, 0.125, 0.145833333333333, 0.166666666666667,
    0.1875, 0.208333333333333, 0.229166666666667, 0.25, 0.270833333333333,
    0.291666666666667, 0.3125, 0.333333333333333, 0.354166666666667,
    0.375, 0.395833333333333, 0.416666666666667, 0.4375, 0.458333333333333,
    0.479166666666667, 0.5),

    "12" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25),

    "13" =
    c(0, 0.208333333333333, 0.416666666666667, 0.625, 0.833333333333333,
    1.04166666666667, 1.25, 1.45833333333333, 1.66666666666667, 1.875,
    2.08333333333333, 2.29166666666667, 2.5, 2.70833333333333,
    2.91666666666667, 3.125, 3.33333333333333, 3.54166666666667, 3.75,
    3.95833333333333, 4.16666666666667, 4.375, 4.58333333333333,
    4.79166666666667, 5),

    "14" = c(1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25),

    "15" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25),

    "16" =
    c(0, 0.208333333333333, 0.416666666666667, 0.625, 0.833333333333333,
    1.04166666666667, 1.25, 1.45833333333333, 1.66666666666667, 1.875,
    2.08333333333333, 2.29166666666667, 2.5, 2.70833333333333,
    2.91666666666667, 3.125, 3.33333333333333, 3.54166666666667, 3.75,
    3.95833333333333, 4.16666666666667, 4.375, 4.58333333333333,
    4.79166666666667, 5),

    "17" =
    c(-25, -22.9166666666667, -20.8333333333333, -18.75, -16.6666666666667,
    -14.5833333333333, -12.5, -10.4166666666667, -8.33333333333333,
    -6.25, -4.16666666666666, -2.08333333333333, 0, 2.08333333333334,
    4.16666666666667, 6.25, 8.33333333333334, 10.4166666666667,
    12.5, 14.5833333333333, 16.6666666666667, 18.75, 20.8333333333333,
    22.9166666666667, 25),

    "18" =
    c(2, 2.125, 2.25, 2.375, 2.5, 2.625, 2.75, 2.875, 3, 3.125, 3.25,
    3.375, 3.5, 3.625, 3.75, 3.875, 4, 4.125, 4.25, 4.375, 4.5, 4.625,
    4.75, 4.875, 5),

    "19" =
    c(0, 0.208333333333333, 0.416666666666667, 0.625,
    0.833333333333333, 1.04166666666667, 1.25, 1.45833333333333,
    1.66666666666667, 1.875, 2.08333333333333, 2.29166666666667,
    2.5, 2.70833333333333, 2.91666666666667, 3.125, 3.33333333333333,
    3.54166666666667, 3.75, 3.95833333333333, 4.16666666666667,
    4.375, 4.58333333333333, 4.79166666666667, 5),

    "20" =
    c(0, 0.208333333333333, 0.416666666666667, 0.625, 0.833333333333333,
    1.04166666666667, 1.25, 1.45833333333333, 1.66666666666667,
    1.875, 2.08333333333333, 2.29166666666667, 2.5, 2.70833333333333,
    2.91666666666667, 3.125, 3.33333333333333, 3.54166666666667,
    3.75, 3.95833333333333, 4.16666666666667, 4.375, 4.58333333333333,
    4.79166666666667, 5),

    "21" =
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 21, 22, 23, 24, 25),

    "22" =
    c(0, 0.0416666666666667, 0.0833333333333333, 0.125, 0.166666666666667,
    0.208333333333333, 0.25, 0.291666666666667, 0.333333333333333, 0.375,
    0.416666666666667, 0.458333333333333, 0.5, 0.541666666666667,
    0.583333333333333, 0.625, 0.666666666666667, 0.708333333333333, 0.75,
    0.791666666666667, 0.833333333333333, 0.875, 0.916666666666667,
    0.958333333333333, 1)),

    .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"),

    row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "23", "24", "25"),

    class = "data.frame")


# ------------------------------------------------------------------------------


".TAU" <-
structure(list(

    "1" =
    c(-1, -0.6, -0.333333333333333, -0.142857142857143,
    0, 0.111111111111111, 0.2, 0.272727272727273, 0.333333333333333,
    0.384615384615385, 0.428571428571429, 0.466666666666667, 0.5,
    0.529411764705882, 0.555555555555556, 0.578947368421053, 0.6,
    0.619047619047619, 0.636363636363636, 0.652173913043478,
    0.666666666666667, 0.68, 0.692307692307692, 0.703703703703704,
    0.714285714285714
    ),

    "2" =
    c(-1, 0, 0.333333333333333, 0.5, 0.6, 0.666666666666667,
    0.714285714285714, 0.75, 0.777777777777778, 0.8, 0.818181818181818,
    0.833333333333333, 0.846153846153846, 0.857142857142857,
    0.866666666666667, 0.875, 0.88235294117647, 0.888888888888889,
    0.894736842105263, 0.9, 0.904761904761905, 0.909090909090909,
    0.91304347826087, 0.916666666666667, 0.92),

    "3" =
    c(-0.181725814826518, -0.168930151452714, -0.155798192853549,
    -0.142309156210049, -0.128440099024957, -0.114165590552606,
    -0.0994573153156502, -0.0842835904937131, -0.068608772818993,
    -0.0523925219034918, -0.0355888743571007, -0.0181450645517658,
    0, 0.0189177438301371, 0.0386926132325796, 0.0594257680440222,
    0.0812402882884418, 0.104288760957381, 0.128764787039966,
    0.154921339236023, 0.183102048111355, 0.21379958230518,
    0.247780252512751, 0.286418218456134, 0.333333333333333),

    "4" =
    c(0, 0.5, 0.666666666666667, 0.75, 0.8, 0.833333333333333,
    0.857142857142857, 0.875, 0.888888888888889, 0.9, 0.909090909090909,
    0.916666666666667, 0.923076923076923, 0.928571428571429,
    0.933333333333333, 0.9375, 0.941176470588235, 0.944444444444444,
    0.947368421052632, 0.95, 0.952380952380952, 0.954545454545455,
    0.956521739130435, 0.958333333333333, 0.96),

    "5" =
    c(-0.85052757802554, -0.837983233335134, -0.823159712179848,
    -0.805382359321779, -0.78368703586404, -0.756652338202137,
    -0.722109024177686, -0.676626253020113, -0.61461896491917,
    -0.527006789744252, -0.400406496234527, -0.222118698154441,
    0, 0.222118698154449, 0.400406496234539, 0.527006789744276,
    0.614618964919029, 0.676626253020132, 0.722109024177453,
    0.756652338200781, 0.783687035871101, 0.805382359356256,
    0.823159712267863, 0.837983231749698, 0.850527554271354),

    "6" =
    c(0, 0.355065933151777, 0.517962498229816, 0.613705638974404,
    0.677220914237255, 0.722592092430507, 0.756685017415291,
    0.783274098241282, 0.80461673005689, 0.822148933158253,
    0.836832638206725, 0.84932812611196, 0.860110789048376,
    0.869526200860125, 0.877832575748863, 0.885224248904,
    0.891855111133839, 0.897842192832803, 0.903279485909824,
    0.90824351995753, 0.912795085448852, 0.91698501728904,
    0.920858299365945, 0.924445190119985, 0.927779794217425),

    "7" =
    c(1, 0.971927944913947, 0.943246768509585, 0.913923522796783,
    0.88392216030227, 0.853203097878133, 0.821722695867944,
    0.789432631089395, 0.756279135134686, 0.722202059745913,
    0.687133717127867, 0.65099742284623, 0.613705638880109,
    0.575157568479307, 0.535235982291939, 0.493802937831557,
    0.450693855665945, 0.40570906309108, 0.358601253084469,
    0.309055967047944, 0.256659242461756, 0.200839120747762,
    0.140745344631603, 0.0749411953484011, 0),

    "8" =
    c(-1, -0.333333333333333, -0.111111111111111, 0, 0.0666666666666667,
    0.111111111111111, 0.142857142857143, 0.166666666666667,
    0.185185185185185, 0.2, 0.212121212121212, 0.222222222222222,
    0.230769230769231, 0.238095238095238, 0.244444444444444,
    0.25, 0.254901960784314, 0.259259259259259, 0.263157894736842,
    0.266666666666667, 0.26984126984127, 0.272727272727273,
    0.275362318840580, 0.277777777777778, 0.28),

    "9" =
    c(0, -0.0204163452169608, -0.0400596555238257, -0.0590081036085306,
    -0.0773261331388824, -0.0950679058715638, -0.112279639253442,
    -0.129001262402105, -0.145267629233813, -0.161109431296128,
    -0.176553899922191, -0.191625356399165, -0.206345649900960,
    -0.220734510872628, -0.234809839618285, -0.24858794447726,
    -0.262083740255211, -0.275310914946985, -0.288282070892981,
    -0.301008845122611, -0.313502012606631, -0.325771575362601,
    -0.337826839765094, -0.349676483955214, -0.361328616888101),

    "10" =
    c(0, -0.0196066744396921,
    -0.0370221472721557, -0.0525703541709941, -0.0665070136005856,
    -0.0790393052012712, -0.0903383149316412, -0.100547418134766,
    -0.109788190716896, -0.118164725554978, -0.125766872846962,
    -0.132672728007605, -0.138950577833938, -0.144660447061511,
    -0.149855344023087, -0.154582275706178, -0.158883083359664,
    -0.162795136572255, -0.166351914405464, -0.16958349544035,
    -0.172516973673804, -0.175176813539516, -0.177585154568996,
    -0.179762074101597, -0.181725814826518),

    "11" =
    c(0, -0.0208398943709387,
    -0.0417175967562695, -0.0626672066307008, -0.083719725295789,
    -0.104903822366714, -0.126246401393688, -0.147773031839223,
    -0.169508288800322, -0.191476027240807, -0.213699608737601,
    -0.236202093256811, -0.259006404912340, -0.282135478276679,
    -0.305612390180847, -0.329460480799121, -0.353703467003579,
    -0.378365550391214, -0.403471521964720, -0.429046865142757,
    -0.455117858555391, -0.481711679925032, -0.50885651222745,
    -0.536581653261374, -0.564917629721708),

    "12" =
    c(0.333333333333333, 0.666666666666667, 0.777777777777778,
    0.833333333333333, 0.866666666666667, 0.888888888888889,
    0.904761904761905, 0.916666666666667, 0.925925925925926,
    0.933333333333333, 0.93939393939394, 0.944444444444444,
    0.948717948717949, 0.952380952380952, 0.955555555555556,
    0.958333333333333, 0.96078431372549, 0.962962962962963,
    0.964912280701754, 0.966666666666667, 0.968253968253968,
    0.96969696969697, 0.971014492753623, 0.972222222222222,
    0.973333333333333),

    "13" =
    c(-0.3613289, -0.269528030161219, -0.187585190523704,
    -0.114164377378166, -0.048139718340646, 0.0114414518639374,
    0.0653882965033201, 0.114390646561491, 0.159038737349337,
    0.199839382405940, 0.237229274320303, 0.271585960700895,
    0.303236932556452, 0.332467174993497, 0.359525461142416,
    0.384629615561405, 0.407970929923072, 0.42971787915222,
    0.450019258484065, 0.469006839692804, 0.486797626862941,
    0.503495777648588, 0.519194244288641, 0.533976179166497,
    0.547916141985897),

    "14" =
    c(0.333333333333333, 0.6, 0.714285714285714, 0.777777777777778,
    0.818181818181818, 0.846153846153846, 0.866666666666667,
    0.88235294117647, 0.894736842105263, 0.904761904761905,
    0.91304347826087, 0.92, 0.925925925925926, 0.93103448275862,
    0.935483870967742, 0.93939393939394, 0.942857142857143,
    0.945945945945946, 0.948717948717949, 0.951219512195122,
    0.953488372093023, 0.955555555555556, 0.957446808510638,
    0.959183673469388, 0.96078431372549),

    "15" =
    c(-1, 0.333333333333333, 0.6, 0.714285714285714,
    0.777777777777778, 0.818181818181818, 0.846153846153846,
    0.866666666666667, 0.88235294117647, 0.894736842105263,
    0.904761904761905, 0.91304347826087, 0.92, 0.925925925925926,
    0.93103448275862, 0.935483870967742, 0.93939393939394,
    0.942857142857143, 0.945945945945946, 0.948717948717949,
    0.951219512195122, 0.953488372093023, 0.955555555555556,
    0.957446808510638, 0.959183673469388),

    "16" =
    c(-1, 0.0199469096156091,
    0.129575836560517, 0.180662881950351, 0.210821233719316,
    0.230868290892863, 0.245206353296857, 0.255989120788036,
    0.264401763304115, 0.271152717969063, 0.276692429117510,
    0.281321400174802, 0.285247984614676, 0.288621292351708,
    0.291550902234041, 0.294119177987248, 0.296389240808115,
    0.298410289410793, 0.300221244635643, 0.301853304387547,
    0.303331771457812, 0.304677385016165, 0.305907306382662,
    0.307035859574841, 0.308075095038758),

    "17" =
    c(-0.505322479883461,
    -0.495828713697966, -0.48454639378008, -0.470935203584076,
    -0.454226630515362, -0.433303941761343, -0.406516652894873,
    -0.371413076160088, -0.324429660012302, -0.26078346047423,
    -0.175313467887867, -0.0654880362471264, 3, 0.198425450290705,
    0.322606327311886, 0.426990238062425, 0.510371695749375,
    0.575835676725875, 0.627386508234144, 0.668494514582698,
    0.701798032118806, 0.729213683422655, 0.752120712677402,
    0.771518248976997, 0.78813985427463),

    "18" =
    c(0.333333333333333, 0.372549019607843, 0.407407407407407,
    0.43859649122807, 0.466666666666667, 0.492063492063492,
    0.515151515151515, 0.536231884057971, 0.555555555555556,
    0.573333333333333, 0.58974358974359, 0.604938271604938,
    0.619047619047619, 0.632183908045977, 0.644444444444444,
    0.655913978494624, 0.666666666666667, 0.676767676767677,
    0.686274509803922, 0.695238095238095, 0.703703703703704,
    0.711711711711712, 0.719298245614035, 0.726495726495726,
    0.733333333333333),

    "19" =
    c(0, 0.429836470415013, 0.492561142661991,
    0.539699842175544, 0.577243238619945, 0.608200347675897,
    0.634340746150618, 0.656827513885057, 0.676426663266239,
    0.693706733577166, 0.709084640531406, 0.722877931761809,
    0.735333844361156, 0.746648476889695, 0.756979782729611,
    0.766456689701118, 0.775185755847842, 0.783255888758933,
    0.790741999213658, 0.797707696363886, 0.80420742554753,
    0.810288035370101, 0.815990188618369, 0.821349301247486,
    0.826396408325272),

    "20" =
    c(0.333333333333333, 0.187581702849446, 0.336923464258114,
    0.453621153661734, 0.544347004922251, 0.615306486593428,
    0.671462346739915, 0.716591196240161, 0.753556906539465,
    0.784548776018356, 0.811245341187264, 0.834925230324345,
    0.856550987713366, 0.87682562788058, 0.896293707460173,
    0.915293066771852, 0.934084281673088, 0.952832982671359,
    0.97164182521641, 0.99056793780414, 1.00963608457645,
    1.02884861875661, 1.04819318064612, 1.06764664223826,
    1.0871818991606),

    "21" =
    c(-0.9999999996, 0.227411277761033, 0.475707247837903,
    0.594420704044238, 0.666780283574186, 0.716296479239256,
    0.752597708588034, 0.780474107458171, 0.80263551556664,
    0.820709018127606, 0.835799583394556, 0.848581734688507,
    0.859631891238008, 0.869228393597159, 0.877745364348898,
    0.885267593021253, 0.8921031172938, 0.89816522529609,
    0.903767761397344, 0.908803416914886, 0.913503328379058,
    0.917746442778645, 0.921749156268669, 0.925374342563027,
    0.92882540945077),

    "22" =
    c(8.88178419700125e-16, -0.0204403642205317, -0.0402325966459149,
    -0.0595398315924127, -0.0784852878388032, -0.09716610551268,
    -0.115661428244375, -0.134037489701166, -0.152350993933844,
    -0.170651459100329, -0.18898289917049, -0.207385066034675,
    -0.225894390250444, -0.244544709719205, -0.263367845858735,
    -0.282394068127122, -0.301652475612564, -0.321171316316645,
    -0.340978259249309, -0.361100630626295, -0.381565622756077,
    -0.402400482266088, -0.423632682913976, -0.445290087203517,
    -0.467401100271068)),

    .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"),

    row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "23", "24", "25"),

    class = "data.frame")


# ------------------------------------------------------------------------------


".RHO" <-
structure(list(
    
    "1" =
        c(-1.00148148148148, -0.738747613322986, -0.466622048681987,
          -0.211687707079451, 0, 0.165652020595619, 0.294857841987463,
          0.396806275669875, 0.478390117460797, 0.544587799346395,
          0.598994031361846, 0.644231561135091, 0.682240753560612,
          0.714478294671625, 0.742053406185035, 0.765822053649082,
          0.786452951364207, 0.80447447147138, 0.820308460074897,
          0.83429494974631, 0.846710452420334, 0.85778166343474,
          0.867695842924806, 0.876608762194626, 0.884650845307784),
    
    "2" =
        c(-1.00148148148148, 0.141567825309872, 0.533448886027939,
          0.708244460527527, 0.800738510405266, 0.855438151990562,
          0.89038551287057, 0.913977650156251, 0.930840822487669,
          0.943140720229384, 0.952367340916965, 0.959636846016786,
          0.965337805120064, 0.969884351339425, 0.973637657583298,
          0.976727483348096, 0.979298582357153, 0.98146084822125,
          0.983296560480555, 0.984868310083549, 0.986245587058835,
          0.987459252328268, 0.988521705438795, 0.989457189777394,
          0.990285276373145),
    
    "3" =
        c(-0.271064557642169, -0.252157028402899, -0.232714526887962,
          -0.212705548849901, -0.19209555265805, -0.170846533523001,
          -0.148916516723668, -0.12625894982673, -0.102821967656877,
          -0.078547495149129, -0.0533701410712804, -0.0272158181938728,
          1.38263885937115e-18, 0.0283745141601889, 0.0580205123349336,
          0.0890702417932126, 0.121680672987631, 0.156040852661039,
          0.192382519137300, 0.230996069414554, 0.272255880785214,
          0.316663434022027, 0.36492869505212, 0.418151365908671,
          0.478390117460797),
    
    "4" =
        c(1.38263885937115e-18, 0.682189978639204, 0.848820654347913,
          0.912515206140176, 0.943205695413621, 0.960245155358946,
          0.970657019434397, 0.977474770619668, 0.982178751195054,
          0.985559960479453, 0.98807177913133, 0.989989026203176,
          0.99148606893177, 0.99267782385585, 0.993642497101149,
          0.994434794399345, 0.99509390407114, 0.995648488592752,
          0.996119899908012, 0.996524305455528, 0.99687412690816,
          0.997179034157432, 0.99744664490479, 0.997683025382027,
          0.997893054234396),
    
    "5" =
        c(-0.972111584358926, -0.967209491068637, -0.960903114161494,
          -0.952607108183223, -0.941402673140619, -0.925789337989117,
          -0.903206982408158, -0.869086491837814, -0.814968529644537,
          -0.725140930480804, -0.573066256464247, -0.328659597722,
          1.38263885937115e-18, 0.328659597722, 0.573066256464246,
          0.725140930480804, 0.814968529644528, 0.86908649183787,
          0.903206982408021, 0.92578933798917, 0.941402673149305,
          0.95260710822185, 0.960903114203512, 0.967209491448638,
          0.972111584081945),
    
    "6" =
        c(1.38263885937115e-18, 0.504193253656214, 0.700093384009142,
          0.798178467968907, 0.854636457142996, 0.890208596051013,
          0.914104060707307, 0.930945457506372, 0.943268127871861,
          0.952561541134428, 0.95968871558049, 0.966247582643763,
          0.970782114931405, 0.976498903835584, 0.979560788569006,
          0.98205552457824, 0.993861933959956, 0.99558219645679,
          1.00597796784080, 1.00765676774247, 1.03633645236535,
          1.03732862637663, 1.03805548919961, 1.05932394767550,
          1.06079424267615),
    
    "7" =
        c(-1.00148148148148, -0.979072702331951, -0.956663923182448,
          -0.934255144032922, -0.910476817558297, -0.885980795610422,
          -0.86087901234568, -0.834111385459534, -0.806575582990396,
          -0.777481481481482, -0.747489711934157, -0.715582990397804,
          -0.682232098765431, -0.646978326474623, -0.60966803840878,
          -0.570271604938272, -0.528215089163238, -0.482955281207133,
          -0.434469135802469, -0.381630727023320, -0.323950617283951,
          -0.259980246913580, -0.187865020576132, -0.104024142661180,
          1.38263885937115e-18),
    
    "8" =
        c(-1.00148148148148, -0.382286405036925, -0.114601585715482,
          0.0310996982933278, 0.121547108159812, 0.182711347208507,
          0.226564943368020, 0.259498634290694, 0.28500277910381,
          0.305383082451978, 0.321951896975906, 0.335684240420725,
          0.347297689031008, 0.35716340065043, 0.365682214112401,
          0.373072635634355, 0.379610561045101, 0.385393549400538,
          0.390535446981844, 0.395136786179741, 0.399301204243857,
          0.403072768984719, 0.406501170036244, 0.409631031685417,
          0.412499583816874),
    
    "9" =
        c(0, -0.0306009955864517, -0.0600171715566712, -0.0883435022143309,
          -0.115662851324462, -0.14204808424762, -0.167563722488530,
          -0.192267256145817, -0.216210197252217, -0.239438934387163,
          -0.261995433039874, -0.283917814885177, -0.305240840989756,
          -0.325996318037027, -0.346213442292882, -0.365919092784125,
          -0.385138082715462, -0.403893376291347, -0.422206276681155,
          -0.440096589759673, -0.457582767389748, -0.474682033331841,
          -0.4914104943232, -0.507783238435509, -0.523814422470025),
    
    "10" =
        c(0, -0.029390926055108, -0.0554838410175942, -0.0787567884244303,
          -0.0995948147001642, -0.118312849424028, -0.135171965480742,
          -0.150391128987132, -0.164155854619207, -0.176624701885803,
          -0.187934236009538, -0.198202876118343, -0.207533922658104,
          -0.216017969401489, -0.223734847121171, -0.230755205988202,
          -0.237141815823774, -0.24295064350982, -0.248231752581214,
          -0.253030059585575, -0.257385974070141, -0.261335943265486,
          -0.264912918148678, -0.268146754209015, -0.271064557642169),
    
    "11" =
        c(0, -0.0312337307463261, -0.0624885635341056, -0.0937774571924054,
          -0.125109119631511, -0.156488211337586, -0.187915529826987,
          -0.219388201730151, -0.250899897834131, -0.282441141124533,
          -0.313999597626781, -0.345560474009285, -0.377106883941072,
          -0.408620285268396, -0.440080745221534, -0.471467680629385,
          -0.502759810642414, -0.533935835408547, -0.564973896022694,
          -0.595854502689236, -0.626555877382724, -0.657058974751787,
          -0.687347683323, -0.717399152001416, -0.747202834006994),
    
    "12" =
        c(0.478390117460797, 0.847457484412861, 0.929514118116192,
          0.959770189940577, 0.974091289816988, 0.98196494553155,
          0.98674930920016, 0.989872618006014, 0.992024786623663,
          0.99357184100672, 0.994722476491513, 0.995602679183581,
          0.99629212083331, 0.996843164325584, 0.99729136172679,
          0.99766153224032, 0.997971424587865, 0.998234001133303,
          0.998458905801764, 0.998653432689402, 0.998823180355318,
          0.998972503148868, 0.99910482845359, 0.99922288351134,
          0.999328860123071),
    
    "13" =
        c(8.60444444444445, -0.396927340010433, -0.279041044592502,
          -0.170790656731972, -0.0721889092723266, 0.0171559644885065,
          0.0978612488308676, 0.170645482690321, 0.236250585187531,
          0.295396788789730, 0.34875879445432, 0.396954947814539,
          0.440544088483541, 0.480026728592929, 0.515848521579894,
          0.548404801961588, 0.578045482558421, 0.605079904444685,
          0.629781421965416, 0.652391617354311, 0.673124105746233,
          0.692167929492244, 0.709690561612045, 0.725840548723069,
          0.740749828053321),
    
    "14" =
        c(0.478390117460797, 0.78697038487367, 0.88669651583035,
          0.930158273783357, 0.952781313104883, 0.966001646844378,
          0.97438037870724, 0.98001921761715, 0.983993663292466,
          0.98689971728531, 0.989089008338032, 0.990779746538743,
          0.992113107410875, 0.99318367420424, 0.99405672549771,
          0.994778481699574, 0.995382405409204, 0.995893200807138,
          0.996329424259493, 0.996705230895822, 0.9970315690086,
          0.997317013020675, 0.99756835474172, 0.997791029818825,
          0.997989429828517),
    
    "15" =
        c(-1.00148148148148, 0.483330421553014, 0.788592827249555,
          0.887213491121673, 0.930351139068034, 0.952863855394297,
          0.966040941166183, 0.974400661352683, 0.980030344652669,
          0.984000050201677, 0.98690350529195, 0.989091305340107,
          0.990781156959334, 0.992113975575947, 0.993184203332592,
          0.99405703905789, 0.994778656516628, 0.995382490232112,
          0.99589322707115, 0.996329412502773, 0.996705194678632,
          0.9970315173572, 0.997316951978555, 0.997568288369395,
          0.997790960846676),
    
    "16" =
        c(-1.00148148148148,
          0.0399918107572716, 0.196922447742505, 0.269217972641498,
          0.311479910536255, 0.339361834072764, 0.359187127162052,
          0.374027186155748, 0.385561510047341, 0.394788531585647,
          0.402340078994373, 0.408635928417183, 0.41396608392397,
          0.418537431930041, 0.422501583532575, 0.425972210130218,
          0.429036230512499, 0.431761254762484, 0.434200670669304,
          0.436397200953502, 0.438385442597004, 0.440193712932073,
          0.441845413865275, 0.443360054973381, 0.444754031079707),
    
    "17" =
        c(-0.644683937053512,
          -0.638242233699325, -0.63002248705626, -0.619327391403488,
          -0.605102498262373, -0.585717699732245, -0.558611529498344,
          -0.519748499835421, -0.462935375218991, -0.379451494373527,
          -0.259494676315209, -0.0980201065759068, -2.86814814814815,
          0.294555695013013, 0.469764827999549, 0.606508141811551,
          0.705924049580753, 0.776198676767287, 0.825867405611722,
          0.86149149285161, 0.887580315101809, 0.907115201815617,
          0.922057497657996, 0.933711847666972, 0.942961892529542),
    
    "18" =
        c(0.579165030761005, 0.612486124687689, 0.641987920637447,
          0.668204867084776, 0.691471975626513, 0.712672417787934,
          0.731485428240074, 0.748541425177017, 0.763998843765579,
          0.7781439158238, 0.79107000000924, 0.802825120690881,
          0.81347401988651, 0.823349095595444, 0.832637848646758,
          0.841114775574677, 0.848922542682376, 0.85620758378877,
          0.862952434000152, 0.869158594153996, 0.87492608635104,
          0.880415774029804, 0.885630004240057,
          0.890455816493047, 0.894963682434538),
    
    "19" =
        c(0.478390117460797, 0.593310952911897, 0.663998561233565,
          0.71419945722106, 0.752174890526936, 0.782049472476609,
          0.806205285404551, 0.826145083754451, 0.842875473608911,
          0.857101027367348, 0.869331955302181, 0.87994792140256,
          0.889237849471334, 0.897425790062242, 0.904688304851229,
          0.911166484840908, 0.916974461375408, 0.922205560262668,
          0.926936833441981, 0.931232450014748, 0.935146270324233,
          0.938723825147773, 0.942003855222286, 0.9450195214372,
          0.947799365319983),
    
    "20" =
        c(0, 0.276598221109225, 0.480091816919932, 0.62207813892774,
          0.720105357077801, 0.788404490421694, 0.83684648649206,
          0.871903512170082, 0.897785391284217, 0.917254967437719,
          0.932155415374907, 0.94373844195542, 0.952871658197302,
          0.960166062946691, 0.96606101324546, 0.970874376612312,
          0.974843755028167, 0.978145796062885, 0.98091616915103,
          0.983256097985314, 0.985247325406535, 0.986952155066644,
          0.988421752891467, 0.98969341361155, 0.990800873367999),
    
    "21" =
        c(-1.00148148148148, 0.347129116118547, 0.65564069685479,
          0.780803624460825, 0.846443027861398, 0.885838814250034,
          0.911560842980429, 0.929365274007805, 0.942234840854385,
          0.951833301777067, 0.96008770231504, 0.965881283889428,
          0.97252949040583, 0.97634531326702, 0.979070799707118,
          0.99154650629721, 0.993599496406914, 1.00432643751157,
          1.00555767252712, 1.03505473383290, 1.03618596470961,
          1.05759086565101, 1.05873221403231, 1.11481166665537,
          1.11507652990229),
    
    "22" =
        c(0, -0.0306367173028787, -0.0602746511903028,
          -0.0891338848924438, -0.117379804930934, -0.145140408516369,
          -0.172516725392443, -0.199589475575909, -0.226423523248459,
          -0.253070965533168, -0.279573384360126, -0.30596355638433,
          -0.332266726350362, -0.358501732266778, -0.384681876941774,
          -0.410815779550173, -0.436907762333632, -0.462959276965944,
          -0.488968041367581, -0.514929666432483, -0.540839175128809,
          -0.566687136331552, -0.592466001518449, -0.618164814543992,
          -0.643774208533738)),
    
    .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"),
    
    row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
        "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
        "23", "24", "25"),
    
    class = "data.frame")


################################################################################

