
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
# FUNCTION:                  KENDALL'S TAU AND SPEARMAN'S RHO:
#  archmTau                   Returns Kendall's tau for Archemedean copulae
#  archmRho                   Returns Spearman's rho for Archemedean copulae
#  .archmTauRange              Returns range for Kendall's tau
#  .archm2Tau                  Alternative way to compute Kendall's tau
#  ### .archmGamma                Returns Gini's gamma for Archimedean copulae
#  .archmTail                  Utility Function
# FUNCTION:                  ARCHIMEDEAN COPULAE TAIL COEFFICIENT:
#  archmTailCoeff             Computes tail dependence for Archimedean copulae
#  archmTailPlot              Plots Archimedean tail dependence function
# REQUIREMENT:               DESCRIPTION:
#  adapt                      Contributed R package adapt
################################################################################


################################################################################
# FUNCTION                   KENDALL'S TAU AND SPEARMAN'S RHO:
#  archmTau                   Returns Kendall's tau for Archemedean copulae
#  archmRho                   Returns Spearman's rho for Archemedean copulae
#  .archmTauRange              Returns range for Kendall's tau
#  .archm2Tau                  Alternative way to compute Kendall's tau
#  .archmGamma                Returns Gini's gamma for Archimedean copulae
#  .archmTail                  Utility Function


archmTau <- 
    function(alpha = NULL, type = archmList(), lower = 1.0e-10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Kendall's tau by integration for Archimedean copulae
    
    # FUNCTION:
    
    # Settings:
    type = match.arg(type)
    Type = as.integer(type)
    
    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Check alpha:
    check = archmCheck(alpha, type)

    # Compute tau: 
    if (length(alpha) == 1) {
        ans = .archmTau(alpha, type, lower)
        names(ans) = "Tau"
        names(alpha) = "alpha"
    } else {
        ans = NULL
        for ( i in 1:length(alpha) )
            ans = c(ans, .archmTau(alpha[i], type, lower)[1]) 
        names(ans) = paste("Tau", 1:length(alpha), sep = "")  
        names(alpha) = paste("alpha", 1:length(alpha), sep = "") 
    }
    
    # Add Control Attribute:
    attr(ans, "control")<-cbind.data.frame(
        t(alpha), type = type, lower = lower, row.names= "")
      
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------
    

.archmTau <-  
    function(alpha = NULL, type = archmList(), lower = 1.0e-10)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Kendall's tau by integration for Archimedean copulae
    
    # FUNCTION:
    
    # Type:
    type <- match.arg(type)
    Type <- as.integer(type)
    
    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Check alpha:
    check = archmCheck(alpha, type)
    
    # Select Type:
    if (Type == 1) {
        if (alpha == -1) return(-1)
        if (alpha == 0) return(0)
        tau = alpha/(alpha+2)
        return(tau)
    } else if (Type == 2) {
        if (alpha == 1) return(-1)
        tau = 1 - 2/alpha
        return(tau)
    } else if (Type == 3 & alpha == 0) {
        return(0)
        # tau numeric
    } else if (Type == 3 & alpha == 1) {
        return(1/3)
        # tau numeric
    } else if (Type == 4) {
        if (alpha == 1) return(0)
        tau = 1 - 1/alpha
        return(tau)
    } else if (Type == 5 & alpha == 0) {
        return(0)
        # tau numeric
    } else if (Type == 6 & alpha == 1) {
        return(0)
    } else if (Type == 7) {
        if (alpha == 0) return(1)
        if (alpha == 1) return(0)
        tau = 2*(1-alpha)*(alpha+log(1-alpha)-alpha*log(1-alpha))/alpha^2
        return(tau)
    } else if (Type == 8) {
        if (alpha == 1) return(-1) 
        tau = (-4+alpha)/(3*alpha)
        return(tau)
    } else if (Type == 9 & alpha == 0) { 
        return(0)
        # tau numeric
    } else if (Type == 10 & alpha == 0) { 
        return(0)
        # tau numeric
    } else if (Type == 11 & alpha == 0) { 
        return(0)
    } else if (Type == 12) {
        tau = 1 - 2/(3*alpha)
        return(tau)
    } else if (Type == 13 & alpha == 1) { 
        return(0)
        # tau numeric
    } else if (Type == 13 & alpha == 0) {
        return(-0.3613289) # 1e-8 value
    } else if (Type == 14) {
        tau = 1 - 4/(2+4*alpha) 
        return(tau)
    } else if (Type == 15) {
        if (alpha == 1) return(-1)
        tau = 1 + 4/(2-4*alpha)
        return(tau)
    } else if (Type == 16 & alpha == 0) { 
        return(-1)
        # tau numeric
    } else if (Type == 17 & alpha == -1) { 
        return(0)
        # tau numeric
    } else if (Type == 18) {
        tau = 1 - 4/(3*alpha)
        return(tau)
    } else if (Type == 19 & alpha == 0) { 
        return(0)
        # tau numeric
    } else if (Type == 20 & alpha == 0) { 
        return(1/3)
        # tau numeric
    # } else if (Type == 21) {
        # tau numeric
    # } else if (Type == 22) {
        # tau numeric
    } else {  
        # Integrate:
        ans = integrate(
            f = .Kfunc, lower = lower, upper = 1, alpha = alpha, type = type,
            stop.on.error = FALSE, rel.tol = .Machine$double.eps^0.5)
        tau = 3 - 4 * ans[[1]]
        attr(tau, "control")<-unlist(c(alpha, type = type, ans[2:4]))
        return(tau)
    }
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.archmTauRange <- 
    function(type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Type:
    type = match.arg(type)
    Type = as.integer(type)
        
    # Range:
    range = matrix( c(
         1,     -1,     1,
         2,     -1,     1,
         3, -0.182,   1/3,
         4,      0,     1,
         5,     -1,     1,
         6,      0,     1,
         7,      0,     1,
         8,     -1,   1/3,
         9,      0, 0.361,
        10,      0, 0.182,
        11,      0,-0.565,
        12,    1/3,     1,
        13,  0.361,    NA,
        14,    1/3,     1,
        15,     -1,     1,
        16,     NA,   1/3,
        17,     -1,     1,
        18,    1/3,     1,
        19,    1/3,     1,
        20,      0,     1,
        21,     NA,    NA,
        22,     NA,    NA ), byrow = TRUE, ncol = 3 )
        
    # Result:
    ans <- range[Type, ][-1]
    names(ans) <- c("tau.lower", "tau.upper")
    attr(ans, "control") <- c(type = type)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.archm2Tau <- 
    function (alpha = NULL, type = archmList(), lower = 1e-6) 
{   
    # A function implemented by Diethelm Wuertz

    # Joe's [1997] alternative expression:
    
    # Type:
    type = match.arg(type)
    Type = as.integer(type)
    
    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Check alpha:
    check = archmCheck(alpha, type)
    
    # Integrate:
    K2func = function(x, alpha, type) {
        x * .invPhiFirstDer(x, alpha, type)^2 } 
    upper = .Phi(0, alpha, type) - lower
    ans = integrate(f = K2func, lower = lower, upper = upper, 
        alpha = alpha, type = type)
    tau = 1 - 4 * ans[[1]]
    attr(tau, "control") <- unlist(c(alpha, type = type, ans[2:4]))
    
    # Return Value:
    tau
}


# ------------------------------------------------------------------------------


archmRho <- 
    function(alpha = NULL, type = archmList(), 
    method = c("integrate2d", "adapt"), error = 1.0e-5)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Spearman's Rho by integration for Archimedean copulae
    
    # FUNCTION:
    
    # Match Arguments:
    method = match.arg(method)
    
    # Type:
    type = match.arg(type)
    Type = as.integer(type)
    
    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Check alpha:
    check = archmCheck(alpha, type)
    
    # Compute Rho:
    if (length(alpha) == 1) {
        ans = .archmRho(alpha, type, method, error)
        names(ans) = "Rho"
        names(alpha) = "alpha"
    } else {
        ans = NULL
        for ( i in 1:length(alpha) )
            ans = c(ans, .archmRho(alpha[i], type, method, error)[1])  
        names(ans) = paste("Rho", 1:length(alpha), sep = "")  
        names(alpha) = paste("alpha", 1:length(alpha), sep = "") 
    }
    
    # Add Control Attribute:
    attr(ans, "control")<-cbind.data.frame(
        t(alpha), type = type, method = method, error = error, row.names= "")
    
    # Return Value:
    ans
}
       

# ------------------------------------------------------------------------------


.archmRho <-
    function(alpha = NULL, type = archmList(), 
    method = c("integrate2d", "adapt"), error = 1.0e-5)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Spearman's rho by integration for Archimedean copulae
    
    # Requirements:
    #   Note, method="adapt" requires R-Package adapt
    
    # FUNCTION:
    
    # Match Arguments:
    method <- match.arg(method)
    
    # Type:
    type <- match.arg(type)
    Type <- as.integer(type)
    
    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Check alpha:
    check = archmCheck(alpha, type)
    
    # Global Parameters:
    ## alpha <<- alpha
    ## type <<- type

    # 2D Integration:
    if (method == "integrate2d" ) {
        # Internal Function :
        fun.integrate2d = 
        function(x, y, alpha, type ) 
        { 
            12 * (.parchm1Copula(x, y, alpha = alpha, type = type) - x*y )
        }
        ans = integrate2d(fun.integrate2d, error = error, alpha = alpha, type = type)
    } else if (method == "adapt") {
        # Requires contributed package adapt ...
        fun.adapt = 
        function(z, alpha, type) 
        { 
            x = z[1]
            y = z[2]
            12 * (.parchm1Copula(x, y, alpha = alpha, type = type) - x*y)
        }
        ans = adapt(ndim = 2, lower = c(0, 0), upper = c(1, 1), 
            minpts = 100, maxpts = NULL, functn = fun.adapt, eps = 0.01,
            alpha = alpha, type = type)
    }
    rho = ans$value
    
    # Result:
    control = list(alpha = alpha[[1]])
    attr(rho, "control") <- unlist(control)
    
    # Return Value:
    rho
}


# ------------------------------------------------------------------------------


# .archmGamma <- 
# function(alpha = 0.5, type = archmList())
# {   # A function implemented by Diethelm Wuertz
# 
#     # Description:
#     #   Gini's gamma by integration for Archimedean copulae
#     
#     # FUNCTION:
#     
#     # Type:
#     type = match.arg(type)
#     Type = as.integer(type)
#     
#     # Check alpha:
#     check = archmCheck(alpha, type)
#     
#     # Specification:
#     spec = copulaSpec("archm", model = list(alpha = alpha, type = type))
#     
#     # Internal Function: 
#     fun = function(x, spec) {
#         f = NULL
#         for ( y in x )
#             f = c( f, 4*(pcopula(y, y, spec) + pcopula(y, 1-y, spec) - y) )
#         f }
#         
#     # Integration:
#     ans = integrate(fun, c(0, 0), c(1, 1), spec = spec)
#     
#     # Result:
#     gamma = ans$value
#     attr(gamma, "control") <- unlist(ans[-1])
#     
#     # Return Value:
#     gamma
# }


################################################################################
# FUNCTION:                  ARCHIMEDEAN COPULAE TAIL COEFFICIENT:
#  archmTailCoeff             Computes tail dependence for Archimedean copulae
#  archmTailPlot              Plots Archimedean tail dependence function


archmTailCoeff <-
    function(alpha = NULL, type = archmList())
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Tail Dependence for Archimedean copulae
    
    # FUNCTION:
    
    # Type:
    type = match.arg(type)
    Type = as.integer(type)
    
    # Alpha:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Check alpha:
    check = archmCheck(alpha, type)

    # Tail Coefficient:
    N = 20
    x = 1 - (1/2)^(1:N)
    lambdaU.Cuv = ( 1 - 2*x + 
        parchmCopula(u = x, v = x, alpha = alpha, type = type) ) / (1-x)
    lambdaU.Phi = 2 - 2 * .invPhiFirstDer(2*x, alpha = alpha, type = type) / 
        .invPhiFirstDer(x, alpha = alpha, type = type)
       
    # Return Value:
    list(lambdaU.Cuv = lambdaU.Cuv, lambdaU.Phi = lambdaU.Phi)
}


# ------------------------------------------------------------------------------


archmTailPlot <-
    function(alpha = NULL, type = archmList(), tail = c("Upper", "Lower"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots tail dependence for elliptical copulae
    
    # Arguments:
    #   rho - a numeric value setting the coorelation strength, ranging
    #       between minus one and one.

    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    Type = as.integer(type)
    tail = match.arg(tail)
    
    # Settings:
    Title = paste("Archimedean Copula No.", 1:22)
    names(Title) = paste(1:22)
    Title = Title[type]
    N = 1000; Points = 20 # don't change these values!
    u = (0:N)/N
    
    # Plot Frame:
    plot(c(0, 1), c(0, 1), type = "n", main = Title, xlab = "u", 
        ylab = paste(tail, "Tail Dependence"))
    
    # Iterate alpha:
    B = 10
    lower = max(archmRange(type)[1], -B)
    upper = min(archmRange(type)[2],  B)
    
    # Select alpha:
    if (is.null(alpha)) {
        # from range:
        Alpha = seq(lower, upper, length = 5)
    } else {
        # from arguments:
        Alpha = alpha
    }
    
    # Do for all alpha:
    for (alpha in Alpha) {
        # Compute Copula Tail dependence lambda:
        C.uu = parchmCopula(u, alpha = alpha, type = type)
        if (tail == "Upper") {
            lambdaTail = (1-2*u+C.uu)/(1-u)
        } else if (tail == "Lower") {
            lambdaTail = C.uu/u
        }
        # Add Parameter Labels:
        text(x = 0.52, y = lambdaTail[floor(N/2)]+0.025, col = "red", 
            cex = 0.7, labels = as.character(round(alpha, 2)))
        # Add Lines:
        lines(u, lambdaTail, lty = 3, col = "black")     
        # Add Points to Curves: 
        if (tail == "Upper") {
            Index = round(seq(1, N-1, length = Points)) 
            X = 1
        } else if (tail == "Lower") {
            Index = round(seq(1, N-1, length = Points)) + 1
            X = 0
        }
        points(u[Index], lambdaTail[Index], col = "steelblue",
            pch = 19, cex = 0.7) 
    }
    abline(h = 0, lty = 3, col = "grey")
    abline(v = X, lty = 3, col = "grey")
    
    # Return Value:
    invisible()
}


################################################################################

