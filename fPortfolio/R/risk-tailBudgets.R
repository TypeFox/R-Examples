
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
#  .tailDependenceCoeffs      Returns Lower and Upper Tail Dependence Coeffs 
# FUNCTION:                  MIXED GUMBEL-SURVIVALGUMBEL-NORMAL COPULA:
#  .rgsgnormCopula            Generates G-SG-NORM copula random variates
#  .dgsgnormCopula            Computes G-SG-NORM copula density
# FUNCTION:                  MIXED G-SG-NORM COPULA FIT:
#  .gsgnormCopulaFit          Estimates the parameters of the G-SG-NORM copula
# FUNCTION:                  NON-PARAMETRIC TAIL DEPENDECY ESTIMATOR:
#  .cfgTDE                    Estimates non-parametrically tail dependence
# FUNCTION:                  COPULA FIT WITH NORM, NIG OR  GHT MARGINALS:
#  .empiricalDependencyFit    Estimates tail dependence with empirical marginals
#  .normDependencyFit         Estimates tail dependence with normal marginals
#  .nigDependencyFit          Estimates tail dependence with NIG marginals  
#  .ghtDependencyFit          Estimates tail dependence with GHT marginals     
################################################################################

 
.tailDependenceCoeffs <-
    function(x, method = c("nig", "norm", "ght"), 
    trace = FALSE, doplot = TRUE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns a list with lower and upper bivariate tail depenence matrixes
    
    # Example:
    #   x <- 100 * LPP2005.RET[, 1:6]; .tailDependenceCoeffs(x)
    
    # Notes:
    #   Tested only for NIG marginal distributions.
    
    # FUNCTION:
    
    # Check Settings:
    method <- match.arg(method)
    
    # Compute Coeffs with desired Marginals:
    fun <- paste(".", method, "DependencyFit", sep = "")
    funFit <- match.fun(fun)
    coeffs <- funFit(x, doplot = doplot, trace = trace) 
    
    # Return Value:
    coeffs
}


################################################################################


.rgsgnormCopula <- 
    function(n = 1000, alpha = c(2, 2), rho = 0, weights = c(1/3, 1/3))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes RVs from a mixed Gumbel-SurvivalGumbel-Normal Copula
    
    # Arguments:
    #   n - an integer value, the number of random variates to be
    #       generated.
    #   alpha - a numeric vector with two entries. The first denotes 
    #       the parameter value of alpha for the Gumbel copula, and
    #       the second for the Survival Gumbel Copula.
    #   rho - a numeric value denoting the correlation parameter for
    #       the normal copula.
    #   weights - a numeric vector with two entries. The first denotes 
    #       the weight of the Gumbel copula, and the second the weight
    #       of the Survival Gumbel Copula. The weight for the normal
    #       copula is evaluated by 1 - sum(weights).
    
    # Example:
    #   .rgsgnormCopula(20)
    
    # FUNCTION:
    
    # Checking:
    stopifnot(any(weights >= 0))
    stopifnot(sum(weights) <= 1)
       
    # Upper Gumbel = 1 , Lower Gumbel = 2, t = 3:
    weights <- c(weights, 1-sum(weights))
    N <- round(n*weights[1:2])
    N <- c(N, n-sum(N))
       
    # Random Variates:
    r <- rbind(
        if (N[1] > 0) fCopulae::rgumbelCopula(N[1], alpha[1]),
        if (N[2] > 0) 1-fCopulae::rgumbelCopula(N[2], alpha[2]),
        if (N[3] > 0) fCopulae::rellipticalCopula(N[3], rho, type = "norm") )
    index <- sample(1:n)
    ans <- r[index, ]
    
    N <- c(n, N)
    names(N) <- c("n", "n1", "n2", "n3")
    attr(ans, "control") < -N
    
    # Return Value:
    ans
} 


# ------------------------------------------------------------------------------


.dgsgnormCopula <- 
    function(u = 0.5, v = u, alpha = c(2, 2), rho = 0, weights = c(1/3, 1/3),
    output = c("vector", "list"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes mixed Gumbel-SurvivalGumbel-Normal Copula density
    
    # Example:
    #   .perspPlot(.dgsgnormCopula(u=grid2d()$x, v=grid2d()$y, output = "list"))
   
    # FUNCTION:
    
    # Settings:
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 1]
        u = u[, 2]
    }
    
    # Match Arguments:
    output = match.arg(output)
    
    # Mixed Copula:
    weights <- c(weights, 1-sum(weights))
    dCopula1 <- fCopulae::dgumbelCopula(u, v, alpha[1])
    dCopula2 <- fCopulae::dgumbelCopula(1-u, 1-v, alpha[2])
    dCopula3 <- fCopulae::dellipticalCopula(u, v, rho, type = "norm")
    c.uv <- weights[1]*dCopula1 + weights[2]*dCopula2 + weights[3]*dCopula3
    
    attr(c.uv, "control") <- c(alpha = alpha, rho = rho, weights = weights)
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        c.uv = list(x = x, y = y, z = matrix(c.uv, ncol = N))
    }
    
    # Return Value:
    c.uv
}


# ------------------------------------------------------------------------------


.gsgnormCopulaFit <- 
    function(u, v =  NULL, trace = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters for a mixed GSG copula
    
    # FUNCTION:
    
    # Settings:
    U <- u
    V <- v
    if (is.list(u)) {
        U = u[[1]]
        V = u[[2]]
    }
    if (is.matrix(u)) {
        U = u[, 1]
        V = u[, 2]
    }
   
    # From weights to gamma ...
    #   G-SG-NORM: W1*Gumbel + W2*SurvivalGumbul + W3*NORM
    #   with W3 = 1-W2-W3
    # Transformation: Use Weights W1 and W3
    #   gamma = c(gamma1, gamma2)
    #       0 < gamma1 = W1/(1-W3) < 1
    #       0 < gamma2 = 1-W3 < 1
    #       W3 = 1-gamma2
    #       W1 = gamma1*gamma2
    #       W2 = gamma2*(1-gamma1)
    # Note, this transformation has the advantage, that gamma1
    #   and gamma2 are independent from each other.
     
    # Estimate Copula:
    start = c(alpha1 = 1.5, alpha2 = 1.5, rho = 0, gamma1 = 1/2, gamma2 = 1/2)
    fun = function(x, U, V, trace) 
    {
        alpha = x[1:2]
        rho = x[3]
        gamma = x[4:5]
        weights = c(
            weights1 = gamma[[1]]*gamma[[2]], 
            weights2 = gamma[[2]]*(1-gamma[[1]]))
        
        density = .dgsgnormCopula(u = U, v = V, alpha, rho, weights = weights)
            
        density = density[!is.na(density)]
        f = -mean( log(density) )
        if (trace) {
            params = round(c(x[1:3], 
                weights[1], weights[2], 1-weights[1]-weights[2]), 4)
            names(params) = c("alpha1", "alpha2", "rho", 
                "gumbel", "survival", "norm")
            cat("\n Objective Function Value:  ", -f)
            cat("\n Parameter Estimates:       ", params, "\n") 
                
        }
        f
    }

    # Fit:
    fit = nlminb(
        start = start, objective = fun, 
        lower = c(  1,   1, -0.999, 0, 0), 
        upper = c(Inf, Inf,  0.999, 1, 1), U = U, V = V, trace = trace)
     
    # Fitted Parameters:   
    param = fit$par
    
    # Named Parameters:
    alpha1 = param[1]
    alpha2 = param[2]
    rho = param[3]
    gamma1 = param[4]
    gamma2 = param[5]
    
    # Weights:
    weights3 =  1-gamma2
    weights1 = gamma1*gamma2
    weights2 = gamma2*(1-gamma1)
    
    # Tail Coefficients
    upperLambda = (weights1*(2 - 2^(1/alpha1)))[[1]]
    lowerLambda = (weights2*(2 - 2^(1/alpha2)))[[1]]
    params = c(param[1:3], param[5]*param[4], param[5]*(1-param[4]), 1-param[5])
    names(params) = c("alpha1", "alpha2", "rho", "gumbel", "survival", "norm")
    Lambda = c(lower = lowerLambda, upper = upperLambda)
      
    # Return Value:
    list(param = params, lambda = Lambda, fitted = fit$par)
}


################################################################################


.cfgTDE <- 
    function(x, y)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Estimates non-parametrically tail dependency coefficient
    
    # FUNCTION:
    
    # Upper Tail:
    lambda = NULL
    n = length(x)
    for(i in 1:n){
        lambda = c(lambda,
            log(sqrt(log(1/x[i])*log(1/y[i]))/log(1/max(x[i],y[i])^2)))
    }
    upper <- 2 - 2*exp(sum(lambda/n))
    
    # Lower Tail:
    x = 1-x
    y = 1-y
    lambda = NULL
    n = length(x)
    for(i in 1:n){
        lambda = c(lambda,
            log(sqrt(log(1/x[i])*log(1/y[i]))/log(1/max(x[i],y[i])^2)))
    }
    lower <- 2 - 2*exp(sum(lambda/n))
    
    # Return Value:
    c(lower = lower, upper = upper)
}


################################################################################


.empiricalDependencyFit <- 
    function(x, doplot = TRUE, trace = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Estimates tail dependency coefficients with Normal marginals
    
    # Arguments:
    #   x - a multivariate 'timeSeries' object
    
    # FUNCTION:
    
    # Settings: 
    N = ncol(x)
    lowerLambda = upperLambda = 0*diag(N)
    assetsNames = colnames(x)
    P = NULL
    
    for (i in 1:(N-1)) {
        # First asset:
        r1 = as.vector(x[, i])
        fit1 = nFit(r1)
        estim1 = fit1@fit$estimate
        p1 = pnorm(r1, estim1[1], estim1[2]) 
        Main1 = assetsNames[i]
        P = cbind(P, p1)
        for (j in (i+1):N) 
        {  
            # Second asset:
            r2 = as.vector(x[, j])
            fit2 = nFit(r2) 
            estim2 = fit2@fit$estimate      
            p2 = pnorm(r2, estim2[1], estim2[2]) 
            Main2 = assetsNames[j]
            # Optional Plot:
            if (doplot) 
            {
                # Plot Distribution:
                MainR = paste("Distribution:", Main1, "-", Main2)
                plot(r1, r2, pch = 19, main = MainR)
                grid()
                
                # Plot Copula:
                MainP = paste("Copula:", Main1, "-", Main2)
                plot(p1, p2, pch = 19, main = MainP)
                grid()
            }
            
            # Fit GSG copula parameters:
            fit = .gsgnormCopulaFit(u = p1, v = p2, trace = trace)
            if (trace)
                cat(assetsNames[c(i,j)], round(fit$lambda, 3), "\n")  
            
                # Compose lambda Matrix:
            lowerLambda[i, j] = lowerLambda[j, i] = fit$lambda[1]
            upperLambda[i, j] = upperLambda[j, i] = fit$lambda[2]
        }
    }
    
    # Result:
    colnames(lowerLambda) = rownames(lowerLambda) = assetsNames
    colnames(upperLambda) = rownames(upperLambda) = assetsNames
    ans = list(lower = lowerLambda, upper = upperLambda)
     
    # Return Value:
    ans 
}


# ------------------------------------------------------------------------------


.normDependencyFit <- 
    function(x, doplot = TRUE, trace = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Estimates tail dependency coefficients with Normal marginals
    
    # Arguments:
    #   x - a multivariate 'timeSeries' object
    
    # FUNCTION:
    
    # Settings: 
    N = ncol(x)
    lowerLambda = upperLambda = 0*diag(N)
    assetsNames = colnames(x)
    P = NULL
    
    for (i in 1:(N-1)) {
        # First asset:
        r1 = as.vector(x[, i])
        fit1 = nFit(r1)
        estim1 = fit1@fit$estimate
        p1 = pnorm(r1, estim1[1], estim1[2]) 
        Main1 = assetsNames[i]
        P = cbind(P, p1)
        for (j in (i+1):N) 
        {  
            # Second asset:
            r2 = as.vector(x[, j])
            fit2 = nFit(r2) 
            estim2 = fit2@fit$estimate      
            p2 = pnorm(r2, estim2[1], estim2[2]) 
            Main2 = assetsNames[j]
            # Optional Plot:
            if (doplot) 
            {
                # Plot Distribution:
                MainR = paste("Distribution:", Main1, "-", Main2)
                plot(r1, r2, pch = 19, main = MainR)
                grid()
                
                # Plot Copula:
                MainP = paste("Copula:", Main1, "-", Main2)
                plot(p1, p2, pch = 19, main = MainP)
                grid()
            }
            
            # Fit GSG copula parameters:
            fit = .gsgnormCopulaFit(u = p1, v = p2, trace = trace)
            if (trace)
                cat(assetsNames[c(i,j)], round(fit$lambda, 3), "\n")  
            # Compose lambda Matrix:
            lowerLambda[i, j] = lowerLambda[j, i] = fit$lambda[1]
            upperLambda[i, j] = upperLambda[j, i] = fit$lambda[2]
        }
    }
    
    # Result:
    colnames(lowerLambda) = rownames(lowerLambda) = assetsNames
    colnames(upperLambda) = rownames(upperLambda) = assetsNames
    ans = list(lower = lowerLambda, upper = upperLambda)
     
    # Return Value:
    ans 
}


# ------------------------------------------------------------------------------


.nigDependencyFit <- 
    function(x, doplot = TRUE, trace = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Estimates tail dependency coefficients with NIG marginals
    
    # Arguments:
    #   x - a multivariate 'timeSeries' object
    
    # FUNCTION:
    
    # Settings:
    N = ncol(x)
    lowerLambda = upperLambda = 0*diag(N)
    assetsNames = colnames(x)
    P = NULL
    
    for (i in 1:(N-1)) {
        # First asset:
        r1 = as.vector(x[, i])
        fit1 = nigFit(r1, doplot = FALSE, trace = trace)
        estim1 = fit1@fit$estimate
        p1 = .pnigC(r1, estim1[1], estim1[2], estim1[3], estim1[4]) 
        Main1 = assetsNames[i]
        P = cbind(P, p1)
        for (j in (i+1):N) {  
            # Second asset:
            r2 = as.vector(x[, j])
            fit2 = nigFit(r2, doplot = FALSE, trace = trace) 
            estim2 = fit2@fit$estimate      
            p2 = .pnigC(r2, estim2[1], estim2[2], estim2[3], estim2[4]) 
            Main2 = assetsNames[j]
            # Optional Plot:
            if (doplot) {
                ## MainR = paste("Distribution:", Main1, "-", Main2)
                ## plot(r1, r2, pch = 19, main = MainR)
                ## grid()
                MainP = paste("Copula:", Main1, "-", Main2)
                plot(p1, p2, pch = 19, main = MainP, xlab = "", ylab = "")
                grid()
            }
            # Fit GSG copula parameters:
            fit = .gsgnormCopulaFit(u = p1, v = p2, trace = trace)
            if (trace)
                cat(assetsNames[c(i,j)], round(fit$lambda, 3), "\n")  
            # Compose lambda Matrix:
            lowerLambda[i, j] = lowerLambda[j, i] = fit$lambda[1]
            upperLambda[i, j] = upperLambda[j, i] = fit$lambda[2]
        }
    }
    
    # Result:
    colnames(lowerLambda) = rownames(lowerLambda) = assetsNames
    colnames(upperLambda) = rownames(upperLambda) = assetsNames
    ans = list(lower = lowerLambda, upper = upperLambda)
     
    # Return Value:
    ans 
}


# ------------------------------------------------------------------------------


.ghtDependencyFit <- 
    function(x, doplot = TRUE, trace = FALSE)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Estimates tail dependency coefficients with GH Student-t marginals
    
    # Arguments:
    #   x - a multivariate 'timeSeries' object
    
    # FUNCTION:
    
    # Settings: 
    N = ncol(x)
    lowerLambda = upperLambda = 0*diag(N)
    assetsNames = colnames(x)
    P = NULL
    
    for (i in 1:(N-1)) {
        # First asset:
        r1 = as.vector(x[, i])
        fit1 = ghtFit(r1, doplot = FALSE, trace = trace)
        estim1 = fit1@fit$estimate
        p1 = pght(r1, estim1[1], estim1[2], estim1[3], estim1[4]) 
        Main1 = assetsNames[i]
        P = cbind(P, p1)
        for (j in (i+1):N) {  
            # Second asset:
            r2 = as.vector(x[, j])
            fit2 = ghtFit(r2, doplot = FALSE, trace = trace) 
            estim2 = fit2@fit$estimate      
            p2 = pght(r2, estim2[1], estim2[2], estim2[3], estim2[4]) 
            Main2 = assetsNames[j]
            # Optional Plot:
            if (doplot) {
                MainR = paste("Distribution:", Main1, "-", Main2)
                plot(r1, r2, pch = 19, main = MainR)
                grid()
                MainP = paste("Copula:", Main1, "-", Main2)
                plot(p1, p2, pch = 19, main = MainP)
                grid()
            }
            # Fit GSG copula parameters:
            fit = .gsgnormCopulaFit(u = p1, v = p2, trace = trace)
            if (trace)
                cat(assetsNames[c(i,j)], round(fit$lambda, 3), "\n")  
            # Compose lambda Matrix:
            lowerLambda[i, j] = lowerLambda[j, i] = fit$lambda[1]
            upperLambda[i, j] = upperLambda[j, i] = fit$lambda[2]
        }
    }
    
    # Result:
    colnames(lowerLambda) = rownames(lowerLambda) = assetsNames
    colnames(upperLambda) = rownames(upperLambda) = assetsNames
    ans = list(lower = lowerLambda, upper = upperLambda)
     
    # Return Value:
    ans 
}


################################################################################

