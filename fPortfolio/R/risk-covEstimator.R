

################################################################################
# FUNCTION:                 DESCRIPTION:
#  covEstimator              Uses sample covariance estimation
#  mveEstimator              Uses robust estimation "cov.mve" from [MASS]
#  mcdEstimator              Uses robust estimation "cov.mcd" from [MASS]
# FUNCTION:                 DESCRIPTION:
#  lpmEstimator              Returns Lower Partial Moment Estimator
#  slpmEstimator             Returns Symmetric Lower Partial Moment Estimator
# FUNCTION:                 DESCRIPTION:
#  kendallEstimator          Returns Kendall's Covariance Estimator
#  spearmanEstimator         Returns Spearman's Covariance Estimator
# FUNCTION:                 DESCRIPTION:
#  covMcdEstimator           Requires "covMcd" from [robustbase]
#  covOGKEstimator           Requires "covOGK" from [robustbase]
#  shrinkEstimator           Requires "cov.shrink" from [corpcor]
#  nnveEstimator             Requires "cov.nnve" from [covRobust]
# FUNCTION:                 DESCRIPTION:
#  .studentEstimator         uses "cov.trob" from [MASS]
#  .baggedEstimator          uses builtin from [corpcor]
#  .donostahEstimator        uses builtin from [robust]
#  .bayesSteinEstimator      copy from Alexios Ghalanos
#  .ledoitWolfEstimator      uses builtin from [tawny]
#  .rmtEstimator             uses builtin from [tawny]
# FUNCTION:                 DESCRIPTION:
#  .mveEstimator2            Uses robust estimation "cov.mve" from [MASS]
#  .mcdEstimator2            Uses robust estimation "cov.mcd" from [MASS]
#  .covMcdEstimator2         Requires "covMcd" from [robustbase]
#  .covOGKEstimator2         Requires "covOGK" from [robustbase]
#  .arwEstimator2            Uses robust estimation ".cov.arw"from [mvoutlier]
################################################################################


covEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses sample covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covEstimator(x)

    # FUNCTION:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = cov(x.mat)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


mveEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses robust estimation "cov.mve" from [MASS]

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; mveEstimator(x)

    # FUNCTION:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate: {
    mu = colMeans(x.mat)
    Sigma = MASS::cov.rob(x = x.mat, method = "mve")$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


mcdEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses robust estimation "cov.mve" from [MASS]

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; mcdEstimator(x)

    # FUNCTION:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = MASS::cov.rob(x = x.mat, method = "mcd")$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


lpmEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns lower partial moment estimator

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; lpmEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu <- colMeans(x)
    if (is.null(spec)) {
        FUN = colMeans
        a = 2
    } else {
        FUN = match.fun(spec@model$param$tau)
        a = spec@model$param$a
    }
    Sigma <- assetsLPM(x, tau = FUN(x), a = a)$Sigma
    colnames(Sigma) <- rownames(Sigma) <- names(mu)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


slpmEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns symmetric lower partial moment estimator

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; slpmEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu <- colMeans(x)
    if (is.null(spec)) {
        FUN = colMeans
        a = 2
    } else {
        FUN = match.fun(spec@model$param$tau)
        a = spec@model$param$a
    }
    Sigma <- assetsSLPM(x, tau = FUN(x), a = a)$Sigma
    colnames(Sigma) <- rownames(Sigma) <- names(mu)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


kendallEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses Kendall's rank covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = cov(x.mat, method = "kendall")

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


spearmanEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses Spearman's rank covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = cov(x.mat, method = "spearman")

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


covMcdEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covMcdEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = robustbase::covMcd(x.mat, alpha = 1/2, ...)$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


covOGKEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covOGKEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = robustbase::covOGK(x.mat, sigmamu = robustbase::scaleTau2, ...)$cov
    colnames(Sigma) <- rownames(Sigma) <- names(mu)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


shrinkEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; shrinkEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = fAssets:::.cov.shrink(x = x.mat, verbose = FALSE, ...)
    attr(Sigma, "lambda.var") <- NULL
    attr(Sigma, "lambda.var.estimated") <- NULL

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


nnveEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Eample:
    #   x  = as.timeSeries(data(LPP2005REC))[, 1:6]; nnveEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    mu = colMeans(x.mat)
    Sigma = fAssets:::.cov.nnve(datamat = x.mat, ...)$cov
    colnames(Sigma) <- rownames(Sigma) <- names(mu)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


.studentEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses mean/student-d covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Note:
    #   Source from package MASS

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; .studentEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    robust = fAssets:::.studentMeanCov(x.mat, ...)
    mu = robust$center
    Sigma = robust$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.baggedEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses bagged mean/covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Note:
    #   Source from package corpcor

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; .baggedEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    robust = fAssets:::.baggedMeanCov(x, ...)
    mu = robust$center
    Sigma = robust$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.donostahEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Description:
    #   Uses Donostah's mean/covariance estimation

    # Note:
    #   Source from package robust

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; .donostahEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    robust = fAssets:::.donostahMeanCov(x, ...)
    mu = robust$center
    Sigma = robust$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.bayesSteinEstimator <-
    function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses Bayes Stein mean/covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Note:
    #   Source from Alexios Ghalanos

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; .bayesSteinEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    robust = fAssets:::.bayesSteinMeanCov(x, ...)
    mu = robust$center
    Sigma = robust$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.ledoitWolfEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses Ledoit-Wolf mean/covariance estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Note:
    #   Source from package tawny

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; .ledoitWolfEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    robust = fAssets:::.ledoitWolfMeanCov(x, ...)
    mu = robust$center
    Sigma = robust$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.rmtEstimator <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses Random Matrix Theory correlation estimation

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Note:
    #   Source from package tawny

    # Eample:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; .rmtEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    robust = fAssets:::.rmtMeanCov(x, ...)
    mu = robust$center
    Sigma = robust$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


.mveEstimator2 <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses robust estimation "cov.mve" from [MASS]

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; mveEstimator(x)

    # FUNCTION:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate: 
    ans = MASS::cov.rob(x = x.mat, method = "mve")
    mu = ans$center
    Sigma = ans$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.mcdEstimator2 <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Uses robust estimation "cov.mve" from [MASS]

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; mcdEstimator(x)

    # FUNCTION:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate: 
    ans = MASS::cov.rob(x = x.mat, method = "mcd")
    mu = ans$center
    Sigma = ans$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


################################################################################


.covMcdEstimator2 <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covMcdEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    ans = robustbase::covMcd(x.mat, alpha = 1/2, ...)
    mu = ans$center
    Sigma = ans$cov

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.covOGKEstimator2 <-
function(x, spec = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:
    #   x - an object of class timeSeries
    #   spec - a portfolio specification of class fPFOLIOSPEC

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))[, 1:6]; covOGKEstimator(x)

    # FUNCTION:

    # Check Arguments:
    stopifnot(inherits(x, "timeSeries"))

    # Extract Matrix:
    x.mat = getDataPart(x)

    # Estimate:
    ans = robustbase::covOGK(x.mat, sigmamu = robustbase::scaleTau2, ...)
    mu = ans$center
    Sigma = ans$cov
    colnames(Sigma) <- rownames(Sigma) <- names(mu)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.arwEstimator2 <-
function (x, spec = NULL, ...) 
{
    x.mat <- as.matrix(x)
    N <- ncol(x)
    assetNames <- colnames(x)
    fit <- fAssets:::.cov.arw(x = x.mat, center = colMeans(x.mat), cov = cov(x) , ...)
    
    # Estimate:
    mu <- fit$center
    Sigma <- fit$cov
    names(mu) <- assetNames
    rownames(Sigma) <- colnames(Sigma) <- assetNames
    
    # Return Value:
    list(mu = fit$center, Sigma = fit$cov)
}


################################################################################

