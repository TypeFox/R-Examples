
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


################################################################################
# FUNCTION:             DESCRIPTION:
#  .filter.RMT           Returns filtered correlation matrix from RMT
#  .mp.density.kernel    Returns kernel density estimate
#  .mp.fit.kernel        Function for fitting the density   
#  .mp.rho               Theoretical density for a set of eigenvalues. 
#  .mp.theory            Calculate and plot the theoretical density distribution
#  .mp.lambdas           Generate eigenvalues for theoretical MP distribution
#  .dmp                  Density in R notation style
################################################################################


# Rmetrics:
#   Note that tawny is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 


# Package: tawny
# Title: Provides various portfolio optimization strategies including
#    random matrix theory and shrinkage estimators
# Version: 1.0
# Date: 2009-03-02
# Author: Brian Lee Yung Rowe
# Maintainer: Brian Lee Yung Rowe <tawny-help@muxspace.com>
# License: GPL-2


# Modifications done by Diethelm Wuertz
# ... works with Rmetrics S4 timeSeries objects
# ... using DEoptim (David Ardia) instead of optim


# ------------------------------------------------------------------------------


.filter.RMT <- 
function(h, trace = TRUE, doplot = TRUE)
{   
    # Description:
    #   Returns filtered correlation matrix from random matrix theory
    
    # Arguments:
    #   h - a multivariate time series object of class timeSeries
    
    # Example:
    #   h = 100 * LPP2005.RET; cor = .filter.RMT(h, FALSE, FALSE)
    
    # FUNCTION:
    
    # Get Data Part:
    h = getDataPart(h)
    
    # .mp.density.kernel()
    # Calculating eigenvalue distribution
    mp.hist <- 
        .mp.density.kernel(h, adjust = 0.2, kernel = 'e', doplot = doplot)
      
    # .mp.fit.kernel()
    # Here we use the DEoptim solver. The reason for this is that the 
    # objective function is not convex, there exist a lot of local minima
    # ... using David Ardia's DEoptim Package
    # DW: To do: modify .DEoptim for a better stop criterion for Q and sigma 
    mp.result <- .DEoptim( 
        FUN = .mp.fit.kernel,
        # Empirically, Q < 0 and sigmas < 0.2 are unrealistic
        lower = c(Q = 0, sigma = 0.2),
        upper = c(10, 10),
        control = list(itermax = 200),
        trace = trace,
        hist = mp.hist)   
    # The solution Q and Sigma: 
    mp.Q <- mp.result$optim$bestmem[1]
    mp.sigma <- mp.result$optim$bestmem[2] 
    if (trace) print(c(mp.Q, mp.sigma))
    
    # Plot:
    if (doplot) rho <- .mp.theory(mp.Q, mp.sigma)
    
    # Cleaning eigenvalues:
    lambda.1 <- mp.hist$values[1]
    sigma.2 <- sqrt(1 - lambda.1/length(mp.hist$values))
    lambda.plus <- sigma.2^2 * (1 + sqrt(1/mp.Q))^2 
    
    # Cleaning correlation matrix: 
    ans = .denoise(mp.hist, lambda.plus, h)
    
    if (trace)  {
        cat("Upper cutoff (lambda.max) is",lambda.plus,"\n")
        cat("Variance is", sigma.2, "\n")
        cat("Greatest eigenvalue is", lambda.1, "\n")
    }
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.mp.density.kernel <- 
function(h, adjust = 0.2, kernel = 'e', doplot = TRUE, ...)
{
    # Description:
    #   Returns kernel density estimate
    
    # Arguments:
    #   h - a multivariate time series object of class timeSeries
    #   adjust, kernel - arguments passed to function density()
    
    # FUNCTION:
    
    # Compute normalized correlation matrix:
    e = cov2cor(cov(h/colSds(h)))
    
    # Calculate eigenvalues
    lambda <- eigen(e, symmetric = TRUE, only.values = FALSE)
    ds <- density(lambda$values, adjust = adjust, kernel = kernel, ...)
    ds$ adjust <- adjust
    ds$kernel <- kernel
    ds$values <- lambda$values
    ds$vectors <- lambda$vectors
    
    # Plot:
    if(doplot) plot(ds, xlim = c(0, max(ds$values)*1.2), 
        main = 'Eigenvalue Distribution')
    
    # Return Value:
    return(ds)
}


# ------------------------------------------------------------------------------


.mp.fit.kernel <- 
function(ps, hist) 
{
    # Description:
    #   Function for fitting the density   
    
    # Arguments:
    #   ps - a numeric vector with two numeric entries, Q and sigma
    #   hist - histogram as returned by the function .mp.density.kernel(h)
    
    # Note:
    #   Calls function .mp.rho()
    
    # FUNCTION:
    
    # Settings:
    BIG <- 1e14
    zeros <- which(hist$y == 0)
    wholes <- which(hist$y > 0)
    after <- head(zeros[zeros > wholes[1]], 1)
    l.plus <- hist$x[after]
    
    Q <- ps[1]
    sigma <- ps[2]     
    rhos <- .mp.rho(Q, sigma, hist$x)
    
    # Just use some very large number to prevent it from being used 
    # as optimal score
    if (max(rhos) == 0) return(BIG)
   
    # Scale densities so that the max values of each are about the same.
    # This is a bit of hand-waving to get the best fit
    scale <- max(rhos) / max(hist$y) + 0.25
    
    # Shift the densities to get a better fit
    whole.idx <- head(rhos[rhos > 0], 1)
    hist$y <- c(
        rep(0, whole.idx-1), 
        tail(hist$y, length(hist$y) - whole.idx+1))
    
    # Normalize based on amount of density below MP upper limit
    # This is basically dividing the distance by the area under 
    # the curve, which gives a bias towards larger areas
    norm.factor <- sum(rhos[hist$x <= l.plus])
    # DW: Check this ...
    hist$y = hist$y[1:length(rhos)]
    dy <- (rhos - (hist$y * scale)) / norm.factor
    
    # Just calculate the distances of densities less than the MP 
    # upper limit
    dist <- as.numeric(dy %*% dy)
    if (is.na(dist)) dist = BIG
    
    # Return Value:
    dist
}


# ------------------------------------------------------------------------------


.mp.rho <- 
function(Q, sigma, e.values)
{
    # Description:
    #   This provides the theoretical density for a set of eigenvalues. 
    #   These are really just points along the x axis for which the 
    #   eigenvalue density is desired.
    
    # Arguments:
    #   Q, sigma - Marcenko-Pastur distribution parameters.
    #   e.values - can be a vector of eigen values or a single eigen value.
    
    # Example:
    #   e.values = seq(-0.5, 4.5, length = 101)
    #   plot(e.values, .mp.rho(2, 1, e.values), type = "h")
    #   points(e.values, .mp.rho(2, 1, e.values), type = "l", col = "red")

    # FUNCTION:
    
    # Get min and max eigenvalues specified by Marcenko-Pastur
    l.min <- sigma^2 * (1 - sqrt(1/Q))^2
    l.max <- sigma^2 * (1 + sqrt(1/Q))^2 

    # Provide theoretical density:
    k <- (Q / 2*pi*sigma^2)
    rho <- k * sqrt(pmax(0, (l.max-e.values)*(e.values-l.min)) ) / e.values
    rho[is.na(rho)] <- 0
    
    # Return Value:
    attr(rho, "e.values") <- e.values
    rho
}


# ------------------------------------------------------------------------------


.mp.theory <- 
function(Q, sigma, e.values = NULL, steps = 200)
{
    # Description:
    #   Calculate and plot the theoretical density distribution
    
    # Arguments:
    #   Q, sigma - Marcenko-Pastur distribution parameters.
    #   e.values - The eigenvalues to plot the density against. 
    #       This can really be any point on the xaxis.
    
    # Note:
    #   calls function .mp.lambdas(), .mp.rho()
    
    # Example:

    # FUNCTION:
    
    # Plot a range of values
    if (is.null(e.values)) { 
        e.values <- .mp.lambdas(Q, sigma, steps) 
    }
    rho <- .mp.rho(Q, sigma, e.values)
    
    if (length(e.values) > 1) {
        l.min <- sigma^2 * (1 - sqrt(1/Q))^2
        l.max <- sigma^2 * (1 + sqrt(1/Q))^2 
        xs <- seq(round(l.min-1), round(l.max+1), (l.max-l.min)/steps)
        main <- paste('Marcenko-Pastur Distribution for Q',Q,'and sigma',sigma)
        plot(xs, rho, xlim = c(0, 6), type = 'l', main = main)
    }
    
    # Return Value:
    rho
}


# ------------------------------------------------------------------------------


.mp.lambdas <- 
function(Q, sigma, steps, trace = FALSE)
{
    # Descrption:
    #   Generate eigenvalues for theoretical Marcenko-Pastur distribution
    
    # Arguments:
    #   Q, sigma - Marcenko-Pastur distribution parameters
    #   steps -
    #   trace - 
    
    # FUNCTION:
    
    # Min and Max Eigenvalues:
    l.min <- sigma^2 * (1 - sqrt(1/Q))^2
    l.max <- sigma^2 * (1 + sqrt(1/Q))^2
    if (trace) {
        cat("min eigenvalue:", l.min, "\n")
        cat("max eigenvalue:", l.max, "\n")}
    
    evs <- seq(round(l.min-1), round(l.max+1), (l.max-l.min)/steps)
    evs[evs < l.min] <- l.min
    evs[evs > l.max] <- l.max
    if (trace) {
        # cat("x labels: ", xs, "\n")
        cat("eigenvalues: ", evs, "\n")
    }
    
    # Return Value:
    evs
}


# ------------------------------------------------------------------------------


.denoise <- 
function(hist, lambda.plus = 1.6, h = NULL)
{
    # Description:
    #   Clean a correlation matrix based on calculated value of lambda.plus 
    #   and the computed eigenvalues.
    #   This takes flattened eigenvalues and returns a new cleaned  
    #   correlation matrix
    
    # Arguments:
    #   e.values - Cleaned eigenvalues
    #   e.vectors - Eigenvectors of correlation matrix of normalized returns
    #   h - non-normalized returns matrix (only used for labels)

    # FUNCTION:
    
    e.values <- hist$values
    avg <- mean(e.values[e.values < lambda.plus])
    e.values[e.values < lambda.plus] <- avg
    
    e.vectors <- hist$vectors
    c.clean <- e.vectors %*% diag(e.values) %*% t(e.vectors)
    diags <- diag(c.clean) %o% rep(1, nrow(c.clean))
    c.clean <- c.clean / sqrt(diags * t(diags))
    
    if (! is.null(h)) {
      rownames(c.clean) <- colnames(h)
      colnames(c.clean) <- colnames(h)
    }
    
    # Return Value:
    c.clean
}


# ------------------------------------------------------------------------------


.dmp = 
function(x, Q = 2, sigma = 1)
{
    # Description:
    #   This provides the theoretical density for a set of eigenvalues. 
    #   These are really just points along the x axis for which the 
    #   eigenvalue density is desired.
    
    # Arguments:
    #   x - 
    #   Q, sigma - Marcenko-Pastur distribution parameters.
    
    # Example:
    #   x = seq(-0.5, 4.5, length = 1001); plot(x, dmp(x, 2, 1), type = "l")

    # FUNCTION:
    
    # Get min and max eigenvalues specified by Marcenko-Pastur
    l.min <- sigma^2 * (1 - sqrt(1/Q))^2
    l.max <- sigma^2 * (1 + sqrt(1/Q))^2 

    # Provide theoretical density:
    k <- (Q / 2*pi*sigma^2)
    rho <- k * sqrt(pmax(0, (l.max-x)*(x-l.min)) ) / x
    rho[is.na(rho)] <- 0 
    
    # Return Value:
    rho 
}


################################################################################

