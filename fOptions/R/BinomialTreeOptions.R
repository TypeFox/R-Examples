
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


###############################################################################
# FUNCTION:                 DESCRIPTION:
#  CRRBinomialTreeOption     Cox-Ross-Rubinstein Binomial Tree Option Model
#  JRBinomialTreeOption      JR Modfication to the Binomial Tree Option
#  TIANBinomialTreeOption    Tian's Modification to the Binomial Tree Option
# FUNCTION:
#  BinomialTreeOption        CRR Binomial Tree Option with Cost of Carry Term
#  BinomialTreePlot          Plots results from the CRR Option Pricing Model
###############################################################################


CRRBinomialTreeOption = 
function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n,
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz           
    
    # Description:
    #   Cox-Ross-Rubinstein Binomial Tree Option Model
    
    # FUNCTION:
    
    # Check Flags:
    TypeFlag = TypeFlag[1]
    z = NA
    if (TypeFlag == "ce" || TypeFlag == "ca") z = +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z = -1
    if (is.na(z)) stop("TypeFlag misspecified: ce|ca|pe|pa")
  
    # Parameters:
    dt = Time/n
    u  = exp(sigma*sqrt(dt))
    d  = 1/u
    p  = (exp(b*dt)-d)/(u-d)
    Df = exp(-r*dt)
    
    # Iteration:
    OptionValue = z*(S*u^(0:n)*d^(n:0) - X)
    OptionValue = (abs(OptionValue) + OptionValue) / 2
    
    # European Option:
    if (TypeFlag == "ce" || TypeFlag == "pe") {
        for ( j in seq(from = n-1, to = 0, by = -1) ) 
            for ( i in 0:j )         
                OptionValue[i+1] = 
                (p*OptionValue[i+2] + (1-p)*OptionValue[i+1]) * Df }
    
    # American Option:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
        for ( j in seq(from = n-1, to = 0, by = -1) )  
            for ( i in 0:j )  
                OptionValue[i+1] = max((z * (S*u^i*d^(abs(i-j)) - X)), 
                    (p*OptionValue[i+2] + (1-p)*OptionValue[i+1]) * Df) }
    
    # Return Value:
    
    
    # Parameters:
    # TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$n = n
    
    # Add title and description:
    if (is.null(title)) title = "CRR Binomial Tree Option"
    if (is.null(description)) description = as.character(date())
    
    # Return Value:
    new("fOPTION", 
        call = match.call(),
        parameters = param,
        price = OptionValue[1], 
        title = title,
        description = description
        )     
}


# ------------------------------------------------------------------------------


JRBinomialTreeOption = 
function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n,
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz           
    
    # Description:
    #   JR Modfication to the Binomial Tree Option
    
    # FUNCTION:
    
    # Check Flags:
    TypeFlag = TypeFlag[1]
    if (TypeFlag == "ce" || TypeFlag == "ca") z = +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z = -1
    
    # Parameters:
    dt = Time/n
    # DW Bug Fix: r -> b
    u = exp( (b-sigma^2/2)*dt+sigma*sqrt(dt) )
    d = exp( (b-sigma^2/2)*dt-sigma*sqrt(dt) )
    # DW End of Bug Fix
    p = 1/2
    Df = exp(-r*dt)
    
    # Iteration:
    OptionValue = z*(S*u^(0:n)*d^(n:0) - X)
    OptionValue = (abs(OptionValue) + OptionValue) / 2
    
    # European Option:
    if (TypeFlag == "ce" || TypeFlag == "pe") {
        for ( j in seq(from = n-1, to = 0, by = -1) ) 
            for ( i in 0:j )         
                OptionValue[i+1] = 
                (p*OptionValue[i+2] + (1-p)*OptionValue[i+1]) * Df }
    
                # American Option:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
        for ( j in seq(from = n-1, to=0, by = -1) )  
            for ( i in 0:j )  
                OptionValue[i+1] = max((z * (S*u^i*d^(abs(i-j)) - X)), 
                    (p*OptionValue[i+2] + (1-p)*OptionValue[i+1]) * Df) }
    
    # Return Value:
    OptionValue[1]
    
    # Parameters:
    # TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$n = n
    
    # Add title and description:
    if (is.null(title)) title = "JR Binomial Tree Option"
    if (is.null(description)) description = as.character(date())
    
    # Return Value:
    new("fOPTION", 
        call = match.call(),
        parameters = param,
        price = OptionValue[1], 
        title = title,
        description = description
        )     
}


# ------------------------------------------------------------------------------


TIANBinomialTreeOption = 
function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n,
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz           
    
    # Description:
    #   Tian's Modification to the Binomial Tree Option
    
    # FUNCTION:
    
    # Check Flags:
    TypeFlag = TypeFlag[1]
    if (TypeFlag == "ce" || TypeFlag == "ca") z = +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z = -1  
    
    # Parameters:
    dt = Time/n 
    M = exp ( b*dt )
    V = exp ( sigma^2 * dt )
    u = (M*V/2) * ( V + 1 + sqrt(V*V + 2*V - 3) )
    d = (M*V/2) * ( V + 1 - sqrt(V*V + 2*V - 3) )
    p = (M-d)/(u-d)
    Df = exp(-r*dt)
    
    # Iteration:
    OptionValue = z*(S*u^(0:n)*d^(n:0) - X)
    OptionValue = (abs(OptionValue) + OptionValue) / 2
    
    # European Option:
    if (TypeFlag == "ce" || TypeFlag == "pe") {
        for ( j in seq(from = n-1, to = 0, by = -1) ) 
            for ( i in 0:j )         
                OptionValue[i+1] = 
                (p*OptionValue[i+2] + (1-p)*OptionValue[i+1]) * Df }
    
    # American Option:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
        for ( j in seq(from = n-1, to = 0, by = -1) )  
            for ( i in 0:j )  
                OptionValue[i+1] = max((z * (S*u^i*d^(abs(i-j)) - X)), 
                    (p*OptionValue[i+2] + (1-p)*OptionValue[i+1]) * Df) }
                    
    # Return Value:
    OptionValue[1]
    
    # Parameters:
    # TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n
    param = list()
    param$TypeFlag = TypeFlag
    param$S = S
    param$X = X
    param$Time = Time
    param$r = r
    param$b = b
    param$sigma = sigma
    param$n = n
    
    # Add title and description:
    if (is.null(title)) title = "TIAN Binomial Tree Option"
    if (is.null(description)) description = as.character(date())
    
    # Return Value:
    new("fOPTION", 
        call = match.call(),
        parameters = param,
        price = OptionValue[1], 
        title = title,
        description = description
        )     
}


# ******************************************************************************


BinomialTreeOption = 
function(TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n,
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz           
    
    # Description:
    #   Calculates option prices from the Cox-Ross-Rubinstein
    #   Binomial tree model.
    
    # Note:
    #   The model described here is a version of the CRR Binomial
    #   Tree model. Including a cost of carry term b, the model can
    #   used to price European and American Options on
    #     b=r       stocks
    #     b=r-q     stocks and stock indexes paying a continuous  
    #               dividend yield q
    #     b=0       futures
    #     b=r-rf    currency options with foreign interst rate rf
    
    # Example:
    #   par(mfrow=c(1,1))
    #   Tree = BinomialTree("pa", 100, 95, 0.5, 0.08, 0.08, 0.3, 5)
    #   print(round(Tree, digits=3))
    #   BinomialTreePlot(Tree, main="American Put Option")
    #
    # Reference:
    #   E.G. Haug, The Complete Guide to Option Pricing Formulas, 
    #   1997, Chapter 3.1.1
    
    # FUNCTION:
    
    # Check Flags:
    TypeFlag = TypeFlag[1]
    if (TypeFlag == "ce" || TypeFlag == "ca") z = +1
    if (TypeFlag == "pe" || TypeFlag == "pa") z = -1    
    
    # Parameters:
    dt = Time / n
    u  = exp(sigma*sqrt(dt))
    d  = 1 / u
    p  = (exp(b*dt) - d) / (u - d)
    Df = exp(-r*dt)
    
    # Algorithm:
    OptionValue = z*(S*u^(0:n)*d^(n:0) - X)
    offset = 1
    Tree = OptionValue = (abs(OptionValue)+OptionValue)/2   
    
    # European Type:
    if (TypeFlag == "ce" || TypeFlag == "pe") {
        for (j in (n-1):0) {
            Tree <-c(Tree, rep(0, times=n-j))
            for (i in 0:j) {         
                OptionValue[i+offset] = 
                    (p*OptionValue[i+1+offset] + 
                (1-p)*OptionValue[i+offset]) * Df 
                Tree = c(Tree, OptionValue[i+offset]) } } }
                
    # American Type:
    if (TypeFlag == "ca" || TypeFlag == "pa") {
        for (j in (n-1):0) { 
            Tree <-c(Tree, rep(0, times=n-j))
            for (i in 0:j) { 
                OptionValue[i+offset] = 
                max((z * (S*u^i*d^(abs(i-j)) - X)), 
                        (p*OptionValue[i+1+offset] + 
                (1-p)*OptionValue[i+offset]) * Df ) 
                Tree = c(Tree, OptionValue[i+offset]) } } } 
                
    # Tree-Matrix of form (here n=4):
    # x x x x
    # . x x x
    # . . x x
    # . . . x
    Tree = matrix(rev(Tree), byrow = FALSE, ncol = n+1)
    
    # Tree Output:
    # if (doprint) print(Tree)
    
    # Parameters:
    # TypeFlag = c("ce", "pe", "ca", "pa"), S, X, Time, r, b, sigma, n
    # param = list()
    # param$TypeFlag = TypeFlag
    # param$S = S
    # param$X = X
    # param$Time = Time
    # param$r = r
    # param$b = b
    # param$sigma = sigma
    # param$n = n
    
    # Add title and description:
    # if (is.null(title)) title = "Binomial Tree Option"
    # if (is.null(description)) description = as.character(date())
    
    # Return Value:
    # new("fOPTION", 
    #    call = match.call(),
    #    parameters = param,
    #    price = Tree[1], 
    #    title = title,
    #    description = description
    #    )
    
    # Return Value:
    invisible(Tree)     
}


# ------------------------------------------------------------------------------


BinomialTreePlot = 
function(BinomialTreeValues, dx = -0.025, dy = 0.4, cex = 1, digits = 2, ...) 
{   # A function implemented by Diethelm Wuertz           
    
    # Description:
    #   Plots the binomial tree of the Cox-Ross-Rubinstein
    #   binomial tree model.
    
    # Example:
    #   par(mfrow=c(1,1))
    #   Tree = BinomialTree("a", "p", 100, 95, 0.5, 0.08, 0.08, 0.3, 5)
    #   print(round(Tree, digits = 3))
    #   BinomialTreePlot(Tree, main = "American Put Option")

    # FUNCTION:
    
    # Tree:
    Tree = round(BinomialTreeValues, digits = digits)
    depth = ncol(Tree)
    plot(x = c(1,depth), y = c(-depth+1, depth-1), type = "n", col = 0, ...)
    points(x = 1, y = 0)
    text(1+dx, 0+dy, deparse(Tree[1, 1]), cex = cex)
    for (i in 1:(depth-1) ) {
        y = seq(from = -i, by = 2, length = i+1)
        x = rep(i, times = length(y))+1
        points(x, y, col = 1) 
        for (j in 1:length(x))
            text(x[j]+dx, y[j]+dy, deparse(Tree[length(x)+1-j,i+1]), cex = cex)   
        y = (-i):i
        x = rep(c(i+1,i), times = 2*i)[1:length(y)]
        lines(x, y, col = 2)
    }
    
    # Return Value:
    invisible()
}


# --- 3.1.2 --------------------------------------------------------------------


# Options on a Stock Paying a Known Dividend Yield
# not yet implemented


# --- 3.1.3 --------------------------------------------------------------------


# BarrierBinomialTree
# not yet implemented


# --- 3.1.4 --------------------------------------------------------------------


# ConvertibleBond
# not yet implemented


# --- 3.2 ----------------------------------------------------------------------


# TrinomialTree
# not yet implemented


# --- 3.3 ----------------------------------------------------------------------


# ThreeDimensionalBinomialTree
# PayoffFunction
# not yet implemented


# --- 3.4.1 --------------------------------------------------------------------


# ImpliedBinomialTree
# not yet implemented


# --- 3.4.2 --------------------------------------------------------------------


# ImpliedTrinomialTree
# not yet implemented


# ******************************************************************************

