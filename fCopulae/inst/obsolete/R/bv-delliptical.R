
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             ELLIPTICAL BIVARIATE DISTRIBUTIONS:
#  delliptical2d         Computes density for elliptical distributions
#  .gfunc2d              Generator Function for elliptical distributions
#  .delliptical2dSlider  Slider for bivariate densities
################################################################################


delliptical2d =
function(x, y = x, rho = 0, param = NULL, type = c("norm", "cauchy", "t", 
"logistic", "laplace", "kotz", "epower"), output = c("vector", "list"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Density function for bivariate elliptical distributions
    
    # Arguments:
    #   x, y -  two numeric vectors of the same length.
    #   rho -  a anumeric value specifying the correlation.
    #   param - NULL, a numeric value, or a numeric vector adding
    #       additional parameters to the generator function.
    #   type -  a character string denoting the type of distribution.
    #       This may be either 
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or 
    #       "epower" for the exponential power distribution
        
    # FUNCTION:
    
    # Type:
    type = type[1]
    
    # Settings:
    if (is.list(x)) {
        y = x$y
        x = x$x
    }
    if (is.matrix(x)) {
        y = x[, 2]
        x = x[, 2]
    }

    # Add Default Parameters:
    if (is.null(param)) {
        if (type == "t") param = c(nu = 4)
        if (type == "kotz") param = c(r = sqrt(2))
        if (type == "epower") param = c(r = sqrt(2), s = 1/2)
    }
    
    # Density:
    xoy = ( x^2 - 2*rho*x*y + y^2 ) / (1-rho^2)
    lambda = .gfunc2d(param = param, type = type)[[1]]
    density = lambda * .gfunc2d(x = xoy, param = param, type = type) /
        sqrt(1 - rho^2)
        
    # Add attributes:
    if (is.null(param)) {
        attr(density, "control") = unlist(list(type = type, rho = rho))
    } else {
        attr(density, "control") = unlist(list(type = type, rho = rho, 
            param = param))
    }
    
    # As List ?
    if (output[1] == "list") {
        N = sqrt(length(x))
        x = x[1:N]
        y = matrix(y, ncol = N)[1, ]
        density = list(x = x, y = y,  z = matrix(density, ncol = N))
    }
    
    # Return Value:
    density
}


# ------------------------------------------------------------------------------


.gfunc2d = 
function(x, param = NULL, type = c("norm", "cauchy", "t", "logistic", 
"laplace", "kotz", "epower"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Generator function for elliptical distributions
    
    # Note:
    #   A copy from fExtremes 'gfunc'
    
    # Arguments:
    #   x -  a numeric vector
    #   param - NULL, a numeric value, or a numeric vector adding.
    #       additional parameters to the generator function.
    #   type -  a character string denoting the type of distribution.
    #       This may be either 
    #       "norm" for the normal distribution, or
    #       "cauchy" for the Cauchy distribution, or
    #       "t" for the Student-t distribution, or
    #       "logistic" for the logistic distribution, or
    #       "laplace" for the distribution, or
    #       "kotz" for the original Kotz distribution, or 
    #       "epower" for the exponential power distribution
    
    # Value:
    #   Returns a numeric vector "g(x)" for the generator computed at
    #   the x values taken from the input vector. If x is missing,
    #   the normalizing constant "lambda" will be returned.
    
    # FUNCTION:
    
    # Handle Missing x:
    if (missing(x)) {
        x = NA
        output = "lambda"
    } else {
        output = "g"
    }
    
    # Get Type:
    type = type[1]
    
    # Get Parameters:   
    # if (is.null(param)) param = .ellipticalParam$param
    
    # Create Generator:
    if (type == "norm") {
        g = exp(-x/2)
        lambda = 1 / (2*pi)
        param = NULL
    }
    if (type == "cauchy") {
        g = ( 1 + x )^ (-3/2 )
        lambda = 1 / (2*pi)
        param = NULL
    }
    if (type == "t") {
        if (is.null(param)) {
            nu = 4
        } else {
            nu = param[[1]]
        }
        g = ( 1 + x/nu )^ ( -(nu+2)/2 )
        lambda = 1/(2*pi)
        param = c(nu = nu)
    }
    if (type == "logistic"){
        g = exp(-x/2)/(1+exp(-x/2))^2
        # lambda:
        # integrate(function(x) { exp(-x)/(1+exp(-x))^2}, 0, Inf, 
        #   subdivision = 10000, rel.tol = .Machine$double.eps^0.8)
        # 0.5 with absolute error < 2.0e-13
        lambda = 1 / pi 
        param = NULL
    }
    if (type == "laplace") { # or "double exponential"
        # epower:
        r = sqrt(2)
        s = 1/2
        g = exp(-r*(x/2)^s)
        lambda = s * r^(1/s) / ( 2 * pi * gamma(1/s) )
        param = NULL
    }
    if (type == "kotz") {
        # epower: s = 1
        if (is.null(param)) {
            r = sqrt(2)
        } else {
            r = param
        }
        g = exp(-r*(x/2))
        lambda = r/(2*pi)
        param = c(r = r)
    }
    if (type == "epower") {
        if (is.null(param)) {
            r = sqrt(2) 
            s = 1/2
        } else {
            r = param[[1]]
            s = param[[2]]
        }
        g = exp(-r*(x/2)^s)
        lambda = s * r^(1/s) / ( 2 * pi * gamma(1/s) )
        param = c(r = r, s = s)
    }
    
    # Output:
    output = output[1]
    if (output == "g") {
        ans = g
    } else if (output == "lambda") {
        ans = lambda
    }
    
    # Add Control:
    if (output == "g") {
        attr(ans, "control") = c(type = type, lambda = as.character(lambda))
    } else if (output == "lambda") {
        if (is.null(param)) {
            attr(ans, "control") = unlist(list(type = type))
        } else {
            attr(ans, "control") = unlist(list(type = type, param = param))
        }
    }
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.delliptical2dSlider =
function(B = 10, eps = 1.e-3)
{   # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively perspective plots of density
    
    #FUNCTION:
    
    # Graphic Frame:
    par(mfrow = c(1, 1), cex = 0.7)
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        Distribution = .sliderMenu(no = 1)
        N = .sliderMenu(no = 2)
        rho = .sliderMenu(no = 3)
        nu = .sliderMenu(no = 4)
        r = .sliderMenu(no = 5)
        s = .sliderMenu(no = 6)
        nlev = .sliderMenu(no = 7)
        ncol = .sliderMenu(no = 8)
        if (rho == +1) rho = rho - eps
        if (rho == -1) rho = rho + eps
        
        # Title:
        Names = c("- Normal", "- Cauchy", "- Student t", "- Logistic",
            "- Laplace", "- Kotz", "- Exponential Power")      
        Title = paste("Elliptical Density No:", as.character(Distribution), 
            Names[Distribution], "\nrho = ", as.character(rho)) 
        if (Distribution == 3) Title = paste(Title, "nu =", as.character(nu))
        if (Distribution >= 6) Title = paste(Title, "r =", as.character(r))
        if (Distribution >= 7) Title = paste(Title, "s =", as.character(s))
        
        # Plot: 
        xy= grid2d(x = seq(-5, 5, length = N))
        Type = c("norm", "cauchy", "t", "logistic", "laplace", "kotz", "epower")
        param = NULL
        if (Distribution == 3) param = nu
        if (Distribution == 6) param = r
        if (Distribution == 7) param = c(r, s)
        D = delliptical2d(x = xy, rho = rho, param = param, 
            type = Type[Distribution], output = "list")
        image(D, col = heat.colors(ncol), xlab = "x", ylab = "y" )
        contour(D, nlevels = nlev, add = TRUE)
        title(main = Title)
        
        # Reset Frame:
        par(mfrow = c(1, 1), cex = 0.7)
    }
  
    # Open Slider Menu:
    plot.names = c("Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names = c("Distribution", "N", "rho", "t: nu", "r", "s", plot.names),
        minima      = c(       1,  10,    -1,       1,   0,   0,   10,   12),
        maxima      = c(       7, 100,    +1,       B,   B,   B,  100,  256),
        resolutions = c(       1,  10,   0.1,     0.1, 0.1, 0.1,   10,    1),
        starts      = c(       1,  10,     0,       4,   1,   1,   10,   12)) 
}


################################################################################

