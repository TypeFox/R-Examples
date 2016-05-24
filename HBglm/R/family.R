############################################################################
## Specifies details for regression families                               #
############################################################################

## Each regression will be specified as a named list with the following fields:
##   loglik     - function pointer to log-likelihood
##                Must support 3 args: Xbeta, response, var=NULL 
##   has.tau    - boolean indicating whether tau is needed
##   linkinv    - inverse of the link-function, vectorized (OPTIONAL)
##                Used in predict.
##   deriv      - function pointer to log-likelihood 1st deriv (OPTIONAL)
##   hessian    - function pointer to log-likelihood 2nd deriv (OPTIONAL)
##   constraint - named list with compulsory bounds for beta & tau (OPTIONAL)
##                (use 'upper' and 'lower' as keys for hi/lo bounds)
##                (see hbglm.gaussian & hbglm.binomial for examples)

##############################################################################
# Linear regression
ll_normal <- function(Xbeta, response, var = 1)
{
    #var <- ifelse(is.null(var), 1, var)
    n <- length(Xbeta)
    r <- Xbeta - response
    loglik <- -0.5 * (sum(r^2) / var + n * log(var)) # upto additive const
    return(loglik)
}

hbglm.gaussian <- list(
    loglik = ll_normal,
    has.tau = TRUE,
    linkinv = function(x) return(x),
    constraint = list(tau = list(lower = 0))
)
##############################################################################

##############################################################################
# Logistic regression
ll_logit <- function(eta, response, var = NULL)
{
    exp_eta <- exp(eta)
    f <- - sum((1 - response) * eta + log(1 + 1/exp_eta))
    return(f)
}

hbglm.binomial <- list(
    loglik = ll_logit,
    linkinv = function(x) return(1 / (1 + exp(-x))),
    has.tau = FALSE
)
##############################################################################

##############################################################################
# Poisson regression
ll_poisson <- function(eta, response, var = NULL)
{
    exp_eta <- exp(eta)
    f <- - sum((1 - response) * eta + log(1 + 1/exp_eta))
    return(f)
}

hbglm.poisson <- list(
    loglik = ll_poisson,
    linkinv = function(x) return(exp(x)),
    has.tau = FALSE
)
##############################################################################
