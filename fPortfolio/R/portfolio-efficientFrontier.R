

################################################################################
# FUNCTION:                DESCRIPTION:
#  portfolioFrontier        Returns the efficient frontier of a portfolio
# DEPRECATED:          
#  .portfolioFrontier       old/alt Version
################################################################################


portfolioFrontier <-
function(data, spec = portfolioSpec(), constraints = "LongOnly",
    include.mvl = TRUE, title = NULL, description = NULL)
{
    # A function implemented by Rmetrics

    # Description:
    #   Computes the efficient frontier of a portfolio

    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL
    #   include.mvl - a logical flag, should the minimum variance
    #       locus be added to the plot?

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec = portfolioSpec()
    #   constraints = c("minW[3:4]=0.1", "maxW[5:6]=0.8", "minsumW[1:3]=0.2", "maxsumW[c(2,4)]=0.8")
    #   portfolioFrontier(data, spec, constraints)

    # FUNCTION:

    DEBUG <- TRUE

    # Match Spec Versus Constraints?
    .checkSpecVsConstraints(spec, constraints)

    # Transform Data and Constraints:
    Data <- portfolioData(data, spec)

    # Optimize in N Points the Portfolios along the frontier:
    nFrontierPoints = getNFrontierPoints(spec)

    # The Target Return - get problems in the first and last point for
    #   long only portfolios, just move a little bit aside ...
    mu <- getMu(Data)
    targetReturns <- seq(min(mu), max(mu), length = nFrontierPoints)
    eps <- .Machine$double.eps^0.5
    targetReturns[1] = targetReturns[1]*(1+eps)
    targetReturns[nFrontierPoints] = targetReturns[nFrontierPoints]*(1-eps)

    # How to Go Along the Frontier ?
    #   The Idea is to start from the minvariance portfolio and to explore
    #   the efficient frontier and the minimum variance locus starting from
    #   this point ...
    #   Then we stop when the status flag fails ...

    # Compute minvariancePortfolio:
    mvPortfolio <- minvariancePortfolio(Data, spec, constraints)
    mvReturn <- getTargetReturn(mvPortfolio@portfolio)["mean"]
    minIndex <- which.min(abs(mvReturn-targetReturns))

    # Upper Frontier Part:
    Status <- 0
    IDX <- minIndex
    weights <- targetReturn <- targetRisk <- covRiskBudgets <- maxDD <- NULL
    while (Status == 0 & IDX <= nFrontierPoints) {
        # Add Target Return to Specification:
        setTargetReturn(spec) = targetReturns[IDX]
        # Optimize Efficient Portfolio:
        ans = try(efficientPortfolio(Data, spec, constraints), silent = TRUE)
        if (class(ans) == "try-error") {
            Status = 1
        } else {
            portfolio = ans
            Status = getStatus(portfolio)
        }
        if (Status == 0) {
            Weights = getWeights(portfolio)
            weights = rbind(weights, Weights)
            targetReturn = rbind(targetReturn, getTargetReturn(portfolio@portfolio))
            targetRisk = rbind(targetRisk, getTargetRisk(portfolio@portfolio))
            covRiskBudgets = rbind(covRiskBudgets, getCovRiskBudgets(portfolio@portfolio))
            ### maxDD = rbind(maxDD,
            ###    min(drawdowns(pfolioReturn(data/100, as.vector(Weights)))) )
        }
        IDX = IDX + 1
    }

    # Lower Min Variance Locus:
    if (include.mvl) {
        if (minIndex > 1) {
            weights2 = targetReturn2 = targetRisk2 = covRiskBudgets2 = maxDD2 = NULL
            Status = 0
            IDX = minIndex - 1
            while (Status == 0 & IDX > 0) {
                # Add Target Return to Specification:
                setTargetReturn(spec) = targetReturns[IDX]
                # Optimize Efficient Portfolio:
                ans = try(efficientPortfolio(Data, spec, constraints), silent = TRUE)
                if (class(ans) == "try-error") {
                    Status = 1
                } else {
                    portfolio = ans
                    Status = getStatus(portfolio)
                }
                if (Status == 0) {
                    Weights2 = getWeights(portfolio)
                    weights2 = rbind(Weights2, weights2)
                    targetReturn2 =
                        rbind(getTargetReturn(portfolio@portfolio), targetReturn2)
                    targetRisk2 =
                        rbind(getTargetRisk(portfolio@portfolio), targetRisk2)
                    covRiskBudgets2 =
                        rbind(getCovRiskBudgets(portfolio@portfolio), covRiskBudgets2)
                    ### maxDD2 = rbind(maxDD2, min(drawdowns(
                    ###    pfolioReturn(data/100, as.vector(Weights2)))) )
                }
                IDX = IDX - 1
            }
            weights = rbind(weights2, weights)
            targetReturn = rbind(targetReturn2, targetReturn)
            targetRisk = rbind(targetRisk2, targetRisk)
            covRiskBudgets = rbind(covRiskBudgets2, covRiskBudgets)
            ### maxDD = rbind(maxDD2, maxDD)
        }
    }
    colnames(weights) <- names(getMu(Data))
    rownames(weights) <- NULL
    rownames(covRiskBudgets) <- NULL
    rownames(targetReturn) <- NULL
    rownames(targetRisk) <- NULL

    # Check: Did we find points on the frontier?
    if (is.null(weights)) {
        portfolio <- mvPortfolio
        ### portfolio@portfolio$maxDD = min(drawdowns(
        ### pfolioReturn(data/100, as.vector(getWeights(mvPortfolio)))))
        return(portfolio)
    }

    # Reset Target Return:
    setTargetReturn(spec) <- NULL

    # Call:
    portfolio@call <- match.call()

    # Compose Portfolio:
    portfolio@portfolio <- new("fPFOLIOVAL",
        portfolio = list(
            weights = weights,
            covRiskBudgets = covRiskBudgets,
            targetReturn = targetReturn,
            targetRisk = targetRisk,
            targetAlpha = getAlpha(spec),
            minriskPortfolio = mvPortfolio,
            status = 0))

    ### portfolio@portfolio$maxDD = maxDD

    # Update Title
    portfolio@title = "Portfolio Frontier"

    # Return Value:
    portfolio
}


################################################################################
# DEPRECATED


.portfolioFrontier <-
function(data, spec = portfolioSpec(), constraints = "LongOnly",
    include.mvl = TRUE, title = NULL, description = NULL)
{
    # A function implemented by Rmetrics

    # Description:
    #   Computes the efficient frontier of a portfolio

    # Arguments:
    #   data - a rectangular object of assets
    #   spec - an object of class 'fPFOLIOSPEC'
    #   constraints - a character vector or NULL

    # Example:
    #   data = as.timeSeries(data(LPP2005REC))[, 1:6]
    #   spec = portfolioSpec()
    #   constraints = c("minW[3:4]=0.1", "maxW[5:6]=0.8", "minsumW[1:3]=0.2", "maxsumW[c(2,4)]=0.8")
    #   portfolioFrontier(data, spec, constraints)

    # FUNCTION:

    DEBUG = TRUE

    # Match Spec Versus Constraints?
    .checkSpecVsConstraints(spec, constraints)

    # Transform Data and Constraints:
    Data = portfolioData(data, spec)
    data <- getSeries(Data)

    # Optimize in N Points the Portfolios along the frontier:
    nFrontierPoints = getNFrontierPoints(spec)

    # The Target Return - get problems in the first and last point for
    #   long only portfolios, just move a little bit aside ...
    mu = getMu(Data)
    targetReturns <- seq(min(mu), max(mu), length = nFrontierPoints)
    eps = .Machine$double.eps^0.5
    targetReturns[1] = targetReturns[1]*(1+eps)
    targetReturns[nFrontierPoints] = targetReturns[nFrontierPoints]*(1-eps)

    # How to Go Along the Frontier ?
    #   The Idea is to start from the minvariance portfolio and to explore
    #   the efficient frontier and the minimum variance locus starting from
    #   this point ...
    #   Then we stop when the status flag fails ...

    # Compute minvariancePortfolio:
    mvPortfolio = minvariancePortfolio(Data, spec, constraints)
    mvReturn = getTargetReturn(mvPortfolio)[, "mean"]
    minIndex = which.min(abs(mvReturn-targetReturns))

    # Upper Frontier Part:
    Status = 0
    IDX = minIndex
    weights = targetReturn = targetRisk = covRiskBudgets = maxDD = NULL
    while (Status == 0 & IDX <= nFrontierPoints) {
        # Add Target Return to Specification:
        setTargetReturn(spec) = targetReturns[IDX]
        # Optimize Efficient Portfolio:
        ans = try(efficientPortfolio(Data, spec, constraints), silent = TRUE)
        if (class(ans) == "try-error") {
            Status = 1
        } else {
            portfolio = ans
            Status = getStatus(portfolio)
        }
        if (Status == 0) {
            Weights = getWeights(portfolio)
            weights = rbind(weights, Weights)
            targetReturn = rbind(targetReturn, getTargetReturn(portfolio))
            targetRisk = rbind(targetRisk, getTargetRisk(portfolio))
            covRiskBudgets = rbind(covRiskBudgets, getCovRiskBudgets(portfolio))
            maxDD = rbind(maxDD,
                min(drawdowns(pfolioReturn(data/100, as.vector(Weights)))) )
        }
        IDX = IDX + 1
    }

    # Lower Min Variance Locus:
    if (include.mvl) {
        if (minIndex > 1) {
            weights2 = targetReturn2 = targetRisk2 = covRiskBudgets2 = maxDD2 = NULL
            Status = 0
            IDX = minIndex - 1
            while (Status == 0 & IDX > 0) {
                # Add Target Return to Specification:
                setTargetReturn(spec) = targetReturns[IDX]
                # Optimize Efficient Portfolio:
                ans = try(efficientPortfolio(Data, spec, constraints), silent = TRUE)
                if (class(ans) == "try-error") {
                    Status = 1
                } else {
                    portfolio = ans
                    Status = getStatus(portfolio)
                }
                if (Status == 0) {
                    Weights2 = getWeights(portfolio)
                    weights2 = rbind(Weights2, weights2)
                    targetReturn2 =
                        rbind(getTargetReturn(portfolio), targetReturn2)
                    targetRisk2 =
                        rbind(getTargetRisk(portfolio), targetRisk2)
                    covRiskBudgets2 =
                        rbind(getCovRiskBudgets(portfolio), covRiskBudgets2)
                    maxDD2 = rbind(maxDD2, min(drawdowns(
                        pfolioReturn(data/100, as.vector(Weights2)))) )
                }
                IDX = IDX - 1
            }
            weights = rbind(weights2, weights)
            targetReturn = rbind(targetReturn2, targetReturn)
            targetRisk = rbind(targetRisk2, targetRisk)
            covRiskBudgets = rbind(covRiskBudgets2, covRiskBudgets)
            maxDD = rbind(maxDD2, maxDD)
        }
    }

    # Check: Did we find points on the frontier?
    if (is.null(weights)) {
        portfolio = mvPortfolio
        portfolio@portfolio$maxDD = min(drawdowns(
            pfolioReturn(data/100, as.vector(getWeights(mvPortfolio)))))
        return(portfolio)
    }

    # Reset Target Return:
    setTargetReturn(spec) <- NULL

    # Compose Frontier:
    portfolio@call = match.call()
    portfolio@portfolio$weights  = weights
    portfolio@portfolio$targetReturn = targetReturn
    portfolio@portfolio$targetRisk = targetRisk
    portfolio@portfolio$covRiskBudgets = covRiskBudgets
    portfolio@portfolio$maxDD = maxDD
    portfolio@portfolio$status = 0
    portfolio@portfolio$minriskPortfolio = mvPortfolio
    portfolio@title = "Portfolio Frontier"

    # Return Value:
    portfolio
}


###############################################################################

