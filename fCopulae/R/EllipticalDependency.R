
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
# FUNCTION:                  ELLIPTICAL COPULAE DEPENDENCE MASURES:
#  ellipticalTau              Computes Kendall's tau for elliptical copulae
#  ellipticalRho              Computes Spearman's rho for elliptical copulae
# FUNCTION:                  ELLIPTICAL COPULAE TAIL COEFFICIENT:
#  ellipticalTailCoeff        Computes tail dependence for elliptical copulae
#  ellipticalTailPlot         Plots tail dependence function
################################################################################


################################################################################
# FUNCTION:                  ELLIPTICAL COPULAE DEPENDENCE MASURES:
#  ellipticalTau              Computes Kendall's tau for elliptical copulae
#  ellipticalRho              Computes Spearman's rho for elliptical copulae


ellipticalTau <- 
    function(rho)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Kendall's tau for elliptical copulae
    
    # Arguments:
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.

    # FUNCTION:
    
    # Compute Kendall's Tau:
    ans = 2 * asin(rho) / pi
    if (length(rho) == 1) {
        names(ans) = "Tau"
    } else {
        names(ans) = paste("Tau", 1:length(rho), sep = "")
    }
    
    # Add Control Attribute: 
    attr(ans, "control") = c(rho = rho)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.ellipticalRho <- 
    function(rho, param = NULL, type = ellipticalList(), subdivisions = 500)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Spearman's rho for elliptical copulae
    
    # Arguments:
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.

    # FUNCTION:
    
    # Settings:
    Type = c("Normal Copula", "Cauchy Copula", "Student-t Copula", 
        "Logistic Copula", "Laplace Copula", "Kotz Copula", 
        "Exponential Power Copula")
    names(Type) = c("norm", "cauchy", "t", "logistic", "laplace", 
        "kotz", "epower")
    type = type[1]
    Type = Type[type]
    
    # Compute Spearman's Rho:
    ans.norm = round(6 * asin(rho/2) / pi, 2)
    
    # Spearman's Rho:
    N = subdivisions
    Pi = pfrechetCopula(u = grid2d((1:(N-1))/N), type = "pi", output = "list")
    D = .dellipticalCopulaGrid(N = N, rho = rho, param = param, 
        type = type, border = FALSE)
    ans = round(12*mean(Pi$z*D$z)-3, 2)
    names(ans) = NULL
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


ellipticalRho <- 
    function(rho, param = NULL, type = ellipticalList(), subdivisions = 500)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Spearman's rho for elliptical copulae
    
    # Arguments:
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.

    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    
    # For all Values of rho:
    ans = NULL
    for (i in 1:length(rho)) {
        ans = c(ans, .ellipticalRho(rho[i], param, type, subdivisions))
    }
   
    # Add Control Attribute:
    control = c(
        rho = rho, 
        param = param, 
        type = type, 
        tau = round(2*asin(rho)/pi, 4))
    attr(ans, "control")<-unlist(control)
    if (length(rho) == 1) {
        names(ans) = "Rho"
    } else {
        names(ans) = paste("Rho", 1:length(rho), sep = "")
    }
    
    # Return Value:
    ans
}


################################################################################
# FUNCTION:                  ELLIPTICAL COPULAE TAIL COEFFICIENT:
#  ellipticalTailCoeff        Computes tail dependence for elliptical copulae
#  ellipticalTailPlot         Plots tail dependence function


ellipticalTailCoeff <- 
    function(rho, param = NULL, type = c("norm", "cauchy", "t"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes tail dependence for elliptical copulae
    
    # Arguments:
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.
    
    # Note:
    #   type = c("logistic", "laplace", "kotz", "epower")
    #   not yet implemented

    # FUNCTION:
    
    # Check:
    stopifnot(length(rho) == 1)
    
    # Match Arguments:
    type = match.arg(type)
    
    # Compute Tail Dependence:
    if (type == "norm") {
        lambda = 0
        param = NULL
    }
    if (type == "cauchy") {
        nu = 1 
        arg = sqrt(nu+1) * sqrt(1-rho) / sqrt(1+rho)
        lambda = 2 * (1 - pt(arg, df = nu+1))
        param = NULL
    }
    if (type == "t") {
        nu = param
        if (is.null(nu)) nu = 4
        arg = sqrt(nu+1) * sqrt(1-rho) / sqrt(1+rho)
        lambda = 2 * (1 - pt(arg, df = nu+1))
        param = c(nu = nu)
    }
    if (type == "logistic") {
        lambda = NA
        param = NULL
    }
    if (type == "laplace") {
        lambda = NA
        param = NULL
    }
    if (type == "kotz") {
        lambda = NA
        param = NULL
    }
    if (type == "epower") {
        lambda = NA
        param = NULL
    }
    
    # Result:
    ans = c(lambda = lambda)
    attr(ans, "control") = c(rho = rho, type = type, param = param)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


ellipticalTailPlot <- 
    function(param = NULL, type = c("norm", "cauchy", "t"), 
    tail = c("Lower", "Upper"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots tail dependence for elliptical copulae
    
    # Arguments:
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.
    
    # Note:
    #   type = c("logistic", "laplace", "kotz", "epower")
    #   not yet implemented

    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    tail = match.arg(tail)
    
    # Settings:
    Title = c("Normal", "Cauchy", "Student-t", "Logistic", "Laplace",
        "Kotz", "Exponential Power")
    Title = paste(Title, "Copula")
    names(Title) = c("norm", "cauchy", "t", "logistic", "laplace", 
        "kotz", "epower")
    Title = Title[type]
    tail = tail[1]
    N = 1000; Points = 20 # don't change these values!
    u = (0:N)/N
    SHOW = N+1
    
    # Parameters:
    if (type == "t" & is.null(param)) {
        param = c(nu = 4)
    }
    if (type == "kotz" & is.null(param)) {
        param = c(r = 1)
    }
    if (type == "epower" & is.null(param)) {
        param = c(r = 1, s = 1)
    }
    
    # Plot Frame:
    if (type == "t") 
        Title = paste(Title, "| nu =", as.character(param))
    if (type == "t") 
        Title = paste(Title, "| r =", as.character(param))
    if (type == "epower") 
        Title = paste(Title, "| r =", as.character(param[1]), 
            " s =", as.character(param[2]))
    plot(c(0,1), c(0,1), type = "n", main = Title, xlab = "u", 
        ylab = paste(tail, "Tail Dependence"))
        
    # Cauchy Tail dependence:
    if (type == "cauchy") {
        type = "t"
        param = c(nu = 1)
    }
    
    # Iterate rho:
    Rho = c(-0.99, seq(-0.9, 0.9, by = 0.3), 0.99)
    for (rho in Rho) {
        
        # Compute Tail Coefficient:
        lambda = ellipticalTailCoeff(rho = rho, param = param, type = type)
        
        # Compute Copula Cross Section C(u,u)"
        if (type == "norm") 
            C.uu = pellipticalCopula(u, rho = rho, type = type)
        if (type == "t") 
            C.uu = .ptCopula(u = u, v = u, rho = rho, nu = param)
        if (type == "logistic" | type == "laplace" | type == "kotz" | 
            type == "epower")
            C.uu = .pellipticalCopulaDiag(N, rho = rho, param = param, 
                type = type)$y
        
        # Compute Copula Tail dependence lambda:
        if (tail == "Upper") {
            lambdaTail = (1-2*u+C.uu)/(1-u)
        } else if (tail == "Lower") {
            lambdaTail = C.uu/u
        }
        
        # Define Plot Elements:
        if (abs(rho) < 0.05) {
            color = "black"
            linetype = 1
        } else if (abs(rho) > 0.95) {
            color = "blue" 
            linetype = 1
        } else {
            color = "black"
            linetype = 3
        }
        
        # Normal Tail Dependence:
        if (type == "norm") { 
            lines(u, lambdaTail, lty = linetype, col = color) 
        }
        
        # Cauchy and Student-t Tail Dependence:
        if (type == "t") {
            if (tail == "Upper") 
                lines(u[u < 0.99], lambdaTail[u < 0.99], lty = linetype, 
                    col = color)
            if (tail == "Lower") 
                lines(u[u > 0.01], lambdaTail[u > 0.01], lty = linetype, 
                    col = color)
        }
        
        # Logistic Tail dependence:
        if (type == "logistic" | type == "laplace" | type == "kotz") {
            if (tail == "Lower") {
                SHOW = which.min(lambdaTail[-1])
                ##
                lines(u[SHOW:(N+1)], lambdaTail[SHOW:(N+1)], type = "l", 
                    lty = linetype, col = color)
            }       
            if (tail == "Upper") {
                SHOW = which.min(lambdaTail[-(N+1)])
                lines(u[1:SHOW], lambdaTail[1:SHOW], type = "l", 
                    lty = linetype, col = color)    
            }
        }
        
        # Add rho Labels
        text(x = 0.5, y = lambdaTail[floor(N/2)]+0.05, col = "red", cex = 0.7,
            labels = as.character(round(rho, 2)))
            
        # Add Points to Curves: 
        if (tail == "Upper") {
            M = min(SHOW, N)
            Index = seq(1, M, by = Points)
            X = 1
        } else if (tail == "Lower") {
            M = max(51, SHOW)
            Index = rev(seq(N+1, M, by = -Points))
            X = 0
        }
        points(u[Index], lambdaTail[Index], pch = 19, cex = 0.7)
        
        # Add Tail Coefficient:
        points(x = X, y = lambda[1], pch = 19, col = "red")
        
    }
    points(1, 0, pch = 19, col = "red")
    abline(h = 0, lty = 3, col = "grey")
    abline(v = X, lty = 3, col = "grey")
    
    # Return Value:
    invisible()
}


################################################################################

