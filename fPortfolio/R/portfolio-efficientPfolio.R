
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                     DESCRIPTION:
#  efficientPortfolio            Returns a frontier portfolio
#  maxratioPortfolio             Returns the max return/risk ratio portfolio
#  tangencyPortfolio             Returns the tangency portfolio
#  minriskPortfolio              Returns the minimum risk portfolio
#  minvariancePortfolio          Returns the minimum variance portfolio
#  maxreturnPortfolio            Returns the maximum return portfolio
################################################################################


efficientPortfolio <-
    function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes target risk and weights for an efficient portfolio

    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec = portfolioSpec(); setTargetReturn(spec) <- mean(data)
    #   efficientPortfolio(data, spec)

    # FUNCTION:

    # Match Spec Versus Constraints:
    # .checkSpecVsConstraints(spec, constraints)

    # Optimize Portfolio:
    Solver <- match.fun(getSolver(spec))
    portfolio <- Solver(data, spec, constraints)

    # Set Parameters:
    # Do not use ...
    # setWeights(spec) = portfolio$weights
    # setTargetReturn(spec) = portfolio$targetReturn
    # setTargetRisk(spec) = portfolio$targetRisk
    # to provide overwriting use:
    spec@portfolio$weights <- portfolio$weights
    spec@portfolio$targetReturn <- portfolio$targetReturn
    spec@portfolio$targetRisk <- portfolio$targetRisk

    # Add Status:
    setStatus(spec) <- portfolio$status

    # Add Title:
    Title <- "Efficient Portfolio"

    # Compose Portfolio:
    portfolio <- feasiblePortfolio(data, spec, constraints)
    portfolio@call <- match.call()
    portfolio@title <- Title

    # Return Value:
    portfolio
}


# ------------------------------------------------------------------------------


maxratioPortfolio <-
    function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Capital Market Line

    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL

    # Example:
    #   data <- as.timeSeries(data(LPP2005REC))[, 1:6]
    #   maxratioPortfolio(data)

    # FUNCTION:

    # Match Spec Versus Constraints:
    # .checkSpecVsConstraints(spec, constraints)

    # Transform Data:
    Data <- portfolioData(data, spec)

    # Compute Sharpe ratio to be minimized:
    ratioFun <- function(x, data, spec, constraints)
    {
        # x is the target return ...
        setTargetReturn(spec) <- x[1]
        Solver <- match.fun(getSolver(spec))
        ans <- Solver(data, spec, constraints)
        # 2012-02-21 DW: Return if Solver does not converge
          BIG <- 1e10
          if(ans$status != 0) return(-BIG)
        ratio = (x[1] - getRiskFreeRate(spec)) / ans$objective
        attr(ratio, "weights") <- ans$weights
        attr(ratio, "status") <- ans$status
        return(ratio)
    }

    # Start Solution - Equal Weights Portfolio:
    nAssets <- getNAssets(Data)
    setWeights(spec) <- rep(1/nAssets, times = nAssets)
    fp <- feasiblePortfolio(Data, spec, constraints)
    setTargetReturn(spec) <- getTargetReturn(fp)
    
    ## 2012-03-10 DW: 
    ## tol = 10*.Machine$double.eps - higher tolerance added
    portfolio <- optimize(f = ratioFun, interval = range(getMu(Data)),
        maximum = TRUE, data = Data, spec = spec, constraints = constraints,
        tol = 10*.Machine$double.eps)
        
    ## 2009-04-19 DW:
    ## It may happen, that the maximum ratio portfolio cannot be computed.
    ## One reason is that the portfolio does not exist since the constraints 
    ## are too restrictive. 
    ## Another reason is that the risk free rate is above the highest return
    ## point.
    ## In these cases we stop the computation here, and return an error message.
    
    STATUS = attr(portfolio$objective, "status")
    if (STATUS != 0) {
        # Error Message:
        cat("\nExecution stopped:")   
        cat("\n  The maximum ratio portfolio could not be computed.")
        cat("\nPossible Reason:")
        cat("\n  Your portfolio constraints may be too restrictive.")
        cat("\nStatus Information:")
        cat("\n  status=", STATUS, " from solver ", getSolver(spec), ".", sep = "")
        cat("\n")
        stop(call. = FALSE, show.error.messages = "\n  returned from Rmetrics")
    }
       
    # Continue: Succesfully computed the minimum risk portfolio ...     
    setWeights(spec) <- attr(portfolio$objective, "weights")
    setStatus(spec) <- attr(portfolio$objective, "status")

    # Compose Portfolio:
    portfolio <- feasiblePortfolio(data, spec, constraints)
    portfolio@call <- match.call()
    portfolio@title <- "Max Return/Risk Ratio Portfolio"

    # Return Value:
    portfolio
}


# ------------------------------------------------------------------------------


tangencyPortfolio <-
    function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes Markowitz tangency portfolio
    
    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL
    
    # FUNCTION:
    
    # Portfolio:
    portfolio <- maxratioPortfolio(data, spec, constraints)
    portfolio@title <- "Tangency Portfolio"

    # Return Value:
    portfolio
}


################################################################################


.minriskPortfolio <-
    function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes minimum risk portfolio

    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL

    # Example:
    #   minriskPortfolio(SWX[, 1:3])

    # FUNCTION:

    # Match Spec Versus Constraints:
    #   .checkSpecVsConstraints(spec, constraints)

    # Transform Data:
    Data <- portfolioData(data, spec)

    # Compute target risk to be minimized:
    targetRiskFun <- function(x, data, spec, constraints) {
        # x is the target return ...
        setTargetReturn(spec) = x[1]
        Solver <- match.fun(getSolver(spec))
        ans <- Solver(data, spec, constraints)
        targetRisk <- ans$objective
        attr(targetRisk, "weights") <- ans$weights
        attr(targetRisk, "status") <- ans$status
        return(targetRisk)
    }

    # Minimal Risk:
    portfolio <- optimize(targetRiskFun, interval = range(getMu(Data)),
        data = Data, spec = spec, constraints = constraints, 
        tol = .Machine$double.eps^0.5)
        
    ## 2009-04-19 DW:
    ## It may happen, that the minimum risk protfolio cannot be computed.
    ## One reason is that the portfolio does not exist since the constraints 
    ## are too restrictive. 
    ## In this case we stop the computation here, and return an error message.
    
    STATUS = attr(portfolio$objective, "status")
    if (STATUS != 0) {
        # Error Message:
        cat("\nExecution stopped:")   
        cat("\n  The minimum risk portfolio could not be computed.")
        cat("\nPossible Reason:")
        cat("\n  Your portfolio constraints may be too restrictive.")
        cat("\nStatus Information:")
        cat("\n  status=", STATUS, " from solver ", getSolver(spec), ".", sep = "")
        cat("\n")
        stop(call.= FALSE, show.error.messages = "\n  returned from Rmetrics")
    }
       
    # Continue: Succesfully computed the minimum risk portfolio ...  
    setWeights(spec) <- attr(portfolio$objective, "weights")
    setStatus(spec) <- attr(portfolio$objective, "status")

    # Compose Portfolio:
    portfolio <- feasiblePortfolio(data, spec, constraints)
    portfolio@call <- match.call()
    portfolio@title <- "Minimum Risk Portfolio"

    # Return Value:
    portfolio
}


minriskPortfolio <-
    function(data, spec = portfolioSpec(), constraints = "LongOnly")
{   
    # A function implemented by Diethelm Wuertz
    
    # Note:
    #    NEW VERSION DW
    
    # FUNCTION:
    
    setTargetReturn(spec) <- NULL
    
    portfolio <- efficientPortfolio(data, spec, constraints)
    portfolio@call <- match.call()
    portfolio@title <- "Minimum Risk Portfolio"
    
    # Return Value:
    portfolio
}    


# ------------------------------------------------------------------------------


minvariancePortfolio <-
    function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes global minimum variance portfolio
    
    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL
    
    # FUNCTION:
    
    # Portfolio:
    portfolio <- minriskPortfolio(data, spec, constraints)
    portfolio@title <- "Minimum Variance Portfolio"

    # Return Value:
    portfolio
}


################################################################################


maxreturnPortfolio <-
    function(data, spec=portfolioSpec(), constraints="LongOnly")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes target risk and weights for an efficient portfolio

    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL

    # FUNCTION:

    # Match Spec Versus Constraints:
    # .checkSpecVsConstraints(spec, constraints)

    # Transform Data:
    data = portfolioData(data, spec)

    # Maximize Return:
    if(is.null(getTargetRisk(spec))) {
        stop("Missing target risk for maximum return optimization.")
    } else {
        # Optimize Portfolio:
        Solver = match.fun(getSolver(spec))
        portfolio = Solver(data, spec, constraints)
        setWeights(spec) = portfolio$weights
        setStatus(spec) = portfolio$status
        Title = "Return Maximized Efficient Portfolio"
    }

    # Compose Portfolio:
    portfolio <- feasiblePortfolio(data, spec, constraints)
    portfolio@call <- match.call()
    portfolio@title <- Title

    # Return Value:
    portfolio
}


################################################################################

