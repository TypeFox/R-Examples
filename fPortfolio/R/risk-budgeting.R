
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


###############################################################################
# FUNCTION:                DESCRIPTION:
#  pfolioReturn             Returns portfolio returns
# FUNCTION:                DESCRIPTION:
#  sampleCOV                Returns sample covariance risk
#  normalVaR                Returns normal Value at Risk
#  modifiedVaR              Returns modified Cornish Fisher VaR
#  sampleVaR                Returns sammple VaR from historical quantiles
# FUNCTION:                DESCRIPTION:
#  budgetsSampleCOV         Covariance risk contribution and budgets
#  budgetsNormalVAR         Normal VaR risk contribution and budgets
#  budgetsModifiedVAR       Modified VaR risk contribution and budgets
#  budgetsNormalES          Normal ES (CVaR) risk contribution and budgets
#  budgetsModifiedES        Modified ES (CVaR) risk contribution and budgets
# UTILITIES:               DESCRIPTION:
#  .M34.MM                  Internal fast computing of M3 and M4
#  .run                     Returns execution time information
#  .Ipower                  Internal utility function to compute M3
#  .derIpower               Internal utility function to compute M4
#  .myVaR                   ... to do
# DEPRECATED:              DESCRIPTION:
#  .covarRisk               Computes covariance portfolio risk
#  .mcr                     Computes marginal contribution to covariance risk
#  .mcrBeta                 Computes beta, the rescaled mcr to covariance risk
#  .riskContributions       Computes covariance risk contributions
#  .riskBudgets             Computes covariance risk budgets
###############################################################################


sampleCOV <- 
    function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns sample covariance risk
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # FUNCTION:
    
    # Return Value:
    cov(x)
}


# -----------------------------------------------------------------------------


normalVaR <- 
    function(x, alpha=0.05)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns normal Value at Risk
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # FUNCTION:
    
    # Mean and Centered 2nd Moment:
    x.mean <- colMeans(x)
    x.centered <- t(t(x) - x.mean)
    m2 <- colMeans(x.centered^2)
    
    # Gaussian:
    q <- qnorm(alpha)
    
    # Return Value:
    x.mean + q * sqrt(m2)
}


# -----------------------------------------------------------------------------


modifiedVaR <- 
    function(x, alpha=0.05)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns modified Cornish Fisher VaR
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # Details:
    #     Includes Code Borrowed from Peterson and Boudt, GPL
    
    # FUNCTION:
    
    # Mean and Centered Moments:
    x.mean <- colMeans(x)
    x.centered <- t(t(x) - x.mean)
    m2 <- colMeans(x.centered^2)
    m3 <- colMeans(x.centered^3)
    m4 <- colMeans(x.centered^4)
    skew <- m3 / sqrt(m2^3)
    kurt <- m4 / (m2*m2) - 3
    
    # Cornish Fisher:
    z <- qnorm(alpha)
    q <- z + (z*z-1)*skew/6 + z*(z*z-3)*kurt/24 - z*(2*z*z-5)*skew*skew/36
    
    # Return Value:
    x.mean + q * sqrt(m2)
}


# -----------------------------------------------------------------------------


sampleVaR <- 
    function(x, alpha=0.05)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns sammple VaR from historical quantiles
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # FUNCTION:
    
    # Return Value:
    colQuantiles(x, alpha)
}


# -----------------------------------------------------------------------------


budgetsSampleCOV <- 
    function(x, weights, mu=NULL, Sigma=NULL)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # Details:
    #     Includes Code Borrowed from Peterson and Boudt, GPL
    #     Rmetrics Re-Implementation
    
    # FUNCTION:
    
    # Risk:
    if(is.null(mu)) mu <- colMeans(x)
    if(is.null(Sigma)) Sigma <- cov(x)
    risk <- sqrt( t(weights) %*% Sigma %*% weights )[[1]]
    attr(risk, "estimator") <- substitute(FUN)
    
    # Risk Contributions:
    m2.pfolio <- (t(weights) %*% Sigma %*% weights)[[1]]
    dm2.pfolio <- as.vector(Sigma %*% weights)
    contribution <- dm2.pfolio/sqrt(m2.pfolio) * weights
    names(contribution) <- colnames(x)
    
    # Risk Budgets:
    budgets <- contribution/risk
    names(budgets) <- colnames(x)
    attr(budgets, "control") <- sum(budgets)
    
    # Return Value:
    list(riskCOV=risk, contribution=contribution, budgets=budgets)
}


# -----------------------------------------------------------------------------


budgetsNormalVAR <- 
    function(x, weights, alpha=0.05, mu=NULL, Sigma=NULL)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # Details:
    #     Includes Code Borrowed from Peterson and Boudt, GPL
    
    # FUNCTION:
    
    # Risk:
    if(is.null(mu)) mu <- colMeans(x)
    if(is.null(Sigma)) Sigma <- cov(x)
    risk <- -(t(weights) %*% mu + qnorm(alpha) * 
                sqrt( t(weights) %*% Sigma %*% weights))[[1]]
    attr(risk, "estimator") <- substitute(FUN)
    attr(risk, "alpha") <- alpha
    
    # Risk Contributions:
    m2.pfolio <- (t(weights) %*% Sigma %*% weights)
    dm2.pfolio <- as.vector(Sigma %*% weights)
    contribution <- - (mu + qnorm(alpha)* dm2.pfolio/sqrt(m2.pfolio)) * 
      weights
    names(contribution) <- colnames(x)
    
    # Risk Budgets:
    budgets <- contribution/risk
    names(budgets) <- colnames(x)
    attr(budgets, "sumBudgets") <- sum(budgets)
    
    # Return Value:
    list(riskVAR=risk, contributionVAR=contribution, budgetsVAR=budgets)
}


# -----------------------------------------------------------------------------
    

budgetsModifiedVAR  <- 
    function(x, weights, alpha=0.05, mu=NULL, Sigma=NULL, M3=NULL, M4=NULL) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # Details:
    #     Includes code borrowed from Peterson and Boudt, GPL
    
    # FUNCTION:
    
    # Compute Moments:
    if(is.null(mu)) mu <- colMeans(x)
    if(is.null(Sigma)) Sigma <- cov(x)
    if(is.null(M3) || is.null(M4)) {
        MM <- .M34.MM(x, mu=mu) 
        M3 <- MM$M3
        M4 <- MM$M4
    }
    
    # Risk:
    z <- qnorm(alpha)
    location <- t(weights) %*% mu
    pm2 <- t(weights) %*% Sigma %*% weights 
    dpm2 <- as.vector(2 * Sigma %*% weights)
    pm3 <- weights %*% M3 %*% (weights %x% weights)
    dpm3 <- as.vector(3 * M3 %*% (weights %x% weights))
    pm4 <- t(weights) %*% M4 %*% (weights %x% weights %x% weights) 
    dpm4 <- as.vector(4 * M4 %*% (weights %x% weights %x% weights))
    skew <- (pm3/pm2^(3/2))[[1]]
    exkurt <- (pm4/pm2^(2) - 3)[[1]]
    derskew <- (2 * (pm2^(3/2)) * dpm3 - 3 * pm3 * sqrt(pm2) * dpm2) / 
      (2 * pm2^3)
    derexkurt <- ((pm2) * dpm4 - 2 * pm4 * dpm2)/(pm2^3)  
    h <- z + (1/6) * (z^2 - 1) * skew
    h <- h + (1/24) * (z^3 - 3 * z) * exkurt - (1/36) * (2 * z^3 - 5 * z) * 
      skew^2 
    risk <- -(location + h * sqrt(pm2))
    
    # Risk Contribution:
    derGausVaR <- -as.vector(mu) - qnorm(alpha) * 
        (0.5 * as.vector(dpm2))/sqrt(pm2)
    derMVaR <- derGausVaR + (0.5 * dpm2/sqrt(pm2)) * (-(1/6) * 
        (z^2 - 1) * skew - (1/24) * (z^3 - 3 * z) * exkurt + 
        (1/36) * (2 * z^3 - 5 * z) * skew^2)
    derMVaR <- derMVaR + sqrt(pm2) * (-(1/6) * (z^2 - 1) * derskew - 
        (1/24) * (z^3 - 3 * z) * derexkurt + (1/36) * (2 * z^3 - 
        5 * z) * 2 * skew * derskew)
    contribution <- as.vector(weights) * as.vector(derMVaR)
    names(contribution) <- colnames(x)
    
    # Risk Budgets:
    budgets <- contribution/risk
    names(budgets) <- colnames(x)
    budgets
    attr(budgets, "sum(contribution)-risk") <- sum(contribution)-risk
    attr(budgets, "sum(budgets)") <- sum(budgets)
    
    # Return Value:
    list(modifiedVAR=risk, contribution=contribution, budgets=budgets)
}


# -----------------------------------------------------------------------------


budgetsNormalES <- 
    function(x, weights, alpha=0.05, mu=NULL, Sigma=NULL)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    x - a 'timeSeries' object
    
    # Arguments:
    
    # Details:
    #     Includes Code Borrowed from Peterson and Boudt, GPL
    
    # FUNCTION:
    
    # Risk:
    if(is.null(mu)) mu <- colMeans(x)
    if(is.null(Sigma)) Sigma <- cov(x)
    location <- t(weights) %*% mu
    pm2 <- t(weights) %*% Sigma %*% weights 
    dpm2 <- as.vector(2 * Sigma %*% weights)
    risk <- -location + dnorm(qnorm(alpha)) * sqrt(pm2)/alpha
    attr(risk, "estimator") <- substitute(FUN)
    attr(risk, "alpha") <- alpha
    
    # Contribution:
    derES <- -mu + (1/alpha) * dnorm(qnorm(alpha)) * (0.5 * dpm2)/sqrt(pm2)
    contribution <- weights * derES
    names(contribution) <- colnames(x)
    
    # Budgets:
    budgets <- contribution/risk
    names(budgets) <- colnames(x)
    attr(budgets, "sumBudgets") <- sum(budgets)
    
    # Return Value:
    list(normalES=risk, contribution=contribution, budgets=budgets)
}


# -----------------------------------------------------------------------------


budgetsModifiedES <- 
    function(x, weights, alpha=0.05, mu=NULL, Sigma=NULL, M3=NULL, M4=NULL)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # Details:
    #     Includes code borrowed from Peterson and Boudt, GPL
    
    # FUNCTION:
    
    # Settings:
    if(is.null(mu)) mu <- colMeans(x)
    if(is.null(Sigma)) Sigma <- cov(x)
    if(is.null(M3) || is.null(M4)) {
        MM <- .M34.MM(x, mu=mu) 
        M3 <- MM$M3
        M4 <- MM$M4
    }
    
    # Risk:
    z <- qnorm(alpha)
    location <- t(weights) %*% mu
    pm2 <- (t(weights) %*% Sigma %*% weights)[[1]]
    dpm2 <- as.vector(2 * Sigma %*% weights)
    pm3 <- (weights %*% M3 %*% (weights %x% weights))[[1]]
    dpm3 <- as.vector(3 * M3 %*% (weights %x% weights))
    pm4 <- (t(weights) %*% M4 %*% (weights %x% weights %x% weights))[[1]] 
    dpm4 <- as.vector(4 * M4 %*% (weights %x% weights %x% weights))
    skew <- (pm3/pm2^(3/2))[[1]]
    exkurt <- (pm4/pm2^(2) - 3)[[1]]
    derskew <- (2 * (pm2^(3/2)) * dpm3 - 
        3 * pm3 * sqrt(pm2) * dpm2)/(2 * pm2^3)
    derexkurt <- ((pm2) * dpm4 - 2 * pm4 * dpm2)/(pm2^3)
    h <- z + 
        (1/6) * (z^2 - 1) * skew
    h <- h + 
        (1/24) * (z^3 - 3 * z) * exkurt - 
        (1/36) * (2 * z^3 - 5 * z) * skew^2
    derh <- (1/6) * (z^2 - 1) * derskew + 
        (1/24) * (z^3 - 3 * z) * derexkurt - 
        (1/18) * (2 * z^3 - 5 * z) * skew * derskew
    E <- dnorm(h)
    E <- E + (1/24) * (.Ipower(4, h) - 6 * .Ipower(2, h) + 
        3 * dnorm(h)) * exkurt
    E <- E + (1/6) * (.Ipower(3, h) - 3 * .Ipower(1, h)) * skew
    E <- E + (1/72) * (.Ipower(6, h) - 15 * .Ipower(4, h) + 
        45 * .Ipower(2, h) - 15 * dnorm(h)) * (skew^2)
    E <- E/alpha
    risk <- MES <- -location + sqrt(pm2) * E
    
    # Risk Contributions:
    derMES <- -mu + 0.5 * (dpm2/sqrt(pm2)) * E
    derE <- (1/24) * (.Ipower(4, h) - 6 * .Ipower(2, h) + 
        3 * dnorm(h)) * derexkurt
    derE <- derE + (1/6) * (.Ipower(3, h) - 3 * .Ipower(1, h)) * derskew
    derE <- derE + (1/36) * (.Ipower(6, h) - 15 * .Ipower(4, h) + 
        45 * .Ipower(2, h) - 15 * dnorm(h)) * skew * derskew
    X <- -h * dnorm(h) + (1/24) * (.derIpower(4, h) - 6 * .derIpower(2, h) - 
        3 * h * dnorm(h)) * exkurt
    X <- X + (1/6) * (.derIpower(3, h) - 3 * .derIpower(1, h)) * skew
    X <- X + (1/72) * (.derIpower(6, h) - 15 * .derIpower(4, h) + 
        45 * .derIpower(2, h) + 15 * h * dnorm(h)) * skew^2
    derE <- derE + derh * X
    derE <- derE/alpha
    derMES <- derMES + sqrt(pm2) * derE
    contribution <- as.vector(weights) * as.vector(derMES)
    names(contribution) <- colnames(x)
    
    # Risk Budgets:
    budgets <- contribution/risk
    names(budgets) <- colnames(x)
    attr(budgets, "sumBudgets") <- sum(budgets)
    
    # Return Value:
    list(modifedES=risk, contribution=contribution, budgets=budgets)
}


###############################################################################


.M34.MM <- 
    function (x, mu=NULL)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    #    x - a 'timeSeries' object
    
    # Details:
    #     Includes Code Borrowed from Peterson and Boudt, GPL
    
    # Fast Rmetrics Implementation:
    n <- ncol(x)
    m <- nrow(x)
    if(is.null(mu)) mu <- colMeans(x)
    M3 <- matrix(rep(0, n^3), nrow=n, ncol=n^2)
    M4 <- matrix(rep(0, n^4), nrow=n, ncol=n^3)
    centret <- series(x) - matrix(rep(mu, each=m), ncol=n)
    for (i in c(1:m)) {
        cent <- centret[i, ]
        tcent <- t(cent)
        M <- (cent %*% tcent) %x% tcent
        M3 <- M3 + M
        M4 <- M4 + M %x% tcent
    }
    
    # Return Value:
    list(M3=M3/m, M4=M4/m, mu=mu)
}


# -----------------------------------------------------------------------------


.run <- 
    function(FUN, times=10, mult=100, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    
    # Arguments:
    
    # FUNCTION:
    
    # Timing:
    fun <- match.fun(FUN)
    now <- Sys.time()
    for (i in 1:as.integer(times)) ans <- fun(...)
    done <- Sys.time()
    time <- mult * as.numeric(done - now)
    
    # Print Timing Results:
    cat("Timing:\n")
    print(c(sec=round(time,0), times=times, mult=mult, runs=times*mult))
    cat("\nResults:\n\n")
    
    # Return Value:
    ans
}


# -----------------------------------------------------------------------------


.Ipower <- 
    function (power, h) 
{
    # Description:
    
    # Arguments:
    
    # Details:
    #     A function borrowed from PerformanceAnalytics, GPL
    
    fullprod <- 1
    if ((power%%2) == 0) {
        pstar <- power/2
        for (j in c(1:pstar)) fullprod <- fullprod * (2 * j)
        I <- fullprod * dnorm(h)
        for (i in c(1:pstar)) {
            prod <- 1
            for (j in c(1:i)) prod <- prod * (2 * j)
            I <- I + (fullprod/prod) * (h^(2 * i)) * dnorm(h)
        }
    } else {
        pstar <- (power - 1)/2
        for (j in c(0:pstar)) {
            fullprod = fullprod * ((2 * j) + 1)
        }
        I <- -fullprod * pnorm(h)
        for (i in c(0:pstar)) {
            prod = 1
            for (j in c(0:i)) prod = prod * ((2 * j) + 1)
            I <- I + (fullprod/prod) * (h^((2 * i) + 1)) * dnorm(h)
        }
    }
    return(I)
}


# -----------------------------------------------------------------------------


.derIpower <- 
    function (power, h) 
{
    # Description:
    
    # Arguments:
    
    # Details:
    #     A function borrowed from PerformanceAnalytics, GPL
    
    fullprod <- 1
    if ((power%%2) == 0) {
        pstar <- power/2
        for (j in c(1:pstar)) fullprod = fullprod * (2 * j)
        I <- -fullprod * h * dnorm(h)
        for (i in c(1:pstar)) {
            prod = 1
            for (j in c(1:i)) prod = prod * (2 * j)
            I <- I + (fullprod/prod) * (h^(2 * i - 1)) * 
                (2 * i - h^2) * dnorm(h)
        }
    } else {
        pstar = (power - 1)/2
        for (j in c(0:pstar)) fullprod <- fullprod * ((2 * j) + 1)
        I <- -fullprod * dnorm(h)
        for (i in c(0:pstar)) {
            prod = 1
            for (j in c(0:i)) prod <- prod * ((2 * j) + 1)
            I <- I + (fullprod/prod) * (h^(2 * i) * 
                (2 * i + 1 - h^2)) * dnorm(h)
        }
    }
    return(I)
}


# -----------------------------------------------------------------------------


.myVaR <- 
    function(x, alpha=0.05, method=c("normal", "modified", "sample"))
{
    # A function implemented by Diethelm Wuertz
    
    # todo ...
    
    # Funcion Selection:
    fun <- match.fun(paste(match.arg(method), "VaR", sep=""))
    
    # Return Value:
    fun(x, alpha)
}


###############################################################################
# DEPRECATED - DO NOT REMOVE - REQUIRED BY PACKAGE appRmetricsHandbook

.covarRisk <-
    function(data, weights=NULL, FUN="cov", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes covariance portfolio risk
    
    # Arguments:
    #   data - a multivariate timeSeries object of financial returns
    #   weights - numeric vector of portfolio weights
    #   FUN - a covariance estimator, which returns a matrix of 
    #       covariance estimates, by default the sample covariance
    #   ... - Optional arguments passed to the function FUN
    
    # Example:
    #   covarRisk(data)
    
    # FUNCTION:
    
    # Covariance Risk:
    covFun <- match.fun(FUN)
    COV <- covFun(data)
    
    # Portfolio Weights:
    N <- ncol(COV)
    if (is.null(weights)) weights = rep(1/N, N)
    names(weights) <- colnames(COV)
    
    # Covariance Portfolio Risk:
    covarRisk <- sqrt( t(weights) %*% COV %*% weights )[[1, 1]]
    
    # Return Value:
    covarRisk
}


# -----------------------------------------------------------------------------


.mcr <- 
    function(data, weights=NULL, FUN="cov", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description
    #   Computes marginal contribution to covariance risk
    
    # Arguments:
    #   data - a multivariate timeSeries object of financial returns
    #   weights - numeric vector of portfolio weights
    #   FUN - a covariance estimator, which returns a matrix of 
    #       covariance estimates, by default the sample covariance
    #   ... - Optional arguments passed to the function FUN
    
    # Details:
    #   The formula are implemented according to Goldberg et al., 
    #   see also R script assetsPfolio.R
    
    # References:
    #   Lisa Goldberg et al., Extreme Risk Management, 2009
    #   Scherer and Martin, Introduction to modern portfolio Optimimization
    
    # Example:
    #   data <- assetsSim(100, 6); mcr(data)
    
    # FUNCTION:
    
    # Covariance Risk:
    covFun <- match.fun(FUN)
    COV <- covFun(data)
    N <- ncol(data)
    if (is.null(weights)) weights <- rep(1/N, N)
    
    # Marginal Contribution to Risk  
    mcr <- (COV %*% weights)[, 1] / .covarRisk(data, weights, FUN, ...) 
    names(mcr) <- colnames(data) 
    
    # Return Value:
    mcr
}


# -----------------------------------------------------------------------------


.mcrBeta <- 
    function(data, weights=NULL, FUN="cov", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes beta, the rescaled marginal contribution to covariance risk
    
    # Arguments:
    #   data - a multivariate timeSeries object of financial returns
    #   weights - numeric vector of portfolio weights
    #   FUN - a covariance estimator, which returns a matrix of 
    #       covariance estimates, by default the sample covariance
    #   ... - Optional arguments passed to the function FUN
    
    # Example:
    #    .mcrBeta(data)
   
    # FUNCTION:
    
    # Portfolio Beta:
    beta <- .mcr(data, weights, FUN = FUN, ...) / 
        .covarRisk(data, weights, FUN = FUN, ...)
   
    # Return Value:
    beta
}


# -----------------------------------------------------------------------------


.riskContributions <-
    function(data, weights=NULL, FUN="cov", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes covariance risk contributions
    
    # Arguments:
    #   data - a multivariate timeSeries object of financial returns
    #   weights - numeric vector of portfolio weights
    #   FUN - a covariance estimator, which returns a matrix of 
    #       covariance estimates, by default the sample covariance
    #   ... - Optional arguments passed to the function FUN
    
    # Example:
    #    .riskContributions(data)
    
    # FUNCTION:
    
    # Risk Contributions:
    if (is.null(weights)) {
        N <- ncol(data)
        weights <- rep(1/N, times = N)
    }
    riskContributions <- weights * .mcr(data, weights, FUN, ...)
    
    # Return Value:
    riskContributions
}


# -----------------------------------------------------------------------------


.riskBudgets <-
    function(data, weights=NULL, FUN="cov", ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes covariance risk budgets
    
    # Arguments:
    #   data - a multivariate timeSeries object of financial returns
    #   weights - numeric vector of portfolio weights
    #   FUN - a covariance estimator, which returns a matrix of 
    #       covariance estimates, by default the sample covariance
    #   ... - Optional arguments passed to the function FUN
    
    # Example:
    #    data <- 100*LPP2005.RET[, 1:6]; .riskBudgets(data)
    
    # FUNCTION:
    
    # Risk Budgets:
    riskBudgets <- .riskContributions(data, weights, FUN, ...) /
        .covarRisk(data, weights, FUN, ...)
    
    # Return Value:
    riskBudgets
} 


###############################################################################


