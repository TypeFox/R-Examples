########################################################################
## The following hypothesis tests are provided. 
##
## MLE hypothesis tests:
##     Likelihood ratio test
##     Wald test
##     Rao Score test
##
## IIA hypothesis test 
##     Hausman-McFadden test
##
## The first 2 MLE hypothesis tests are performed by functions from the
## CRAN package 'lmtest'
##
## Code for the Rao Score test and the Hausman-McFadden test is adapted
## from the CRAN package 'mlogit'. 
########################################################################

########################################################################
## Likelihood ratio test
########################################################################
# Calls the lrtest() function from package: lmtest
lrtest.mnlogit <- function(object, ...)
{
    return(lrtest.default(object, ...))
}

########################################################################
## Wald test
########################################################################
# Calls the waldtest() function from package: lmtest
waldtest.mnlogit <- function(object, ...)
{
    return(waldtest.default(object, ...))
}

########################################################################
## Rao Score test or the Lagrange Multiplier Test
########################################################################
scoretest <- function(object, ...)
{
    UseMethod("scoretest")
}

# Code adapted from from the mlogit package function: scoretest.default()
scoretest.mnlogit <- function(object, ...)
{
    dots <- list(...)
    if (length(dots) <= 0)
        stop("Another model must be passed as 2nd arg to scoretest.")
    new <- dots[[1]]
    if (class(new)[1] != "mnlogit")
        stop("Expected fitted mnlogit object as 2nd arg to scoretest.")
    # Swap models to make 'object' the MORE constrained model
    if (length(coef(object)) > length(coef(new))) {
        temp <- object
        object <- new
        new <- temp
        temp <- NULL
    }
    cls <- class(object)[1]
    nmodels <- length(new)
    if (!inherits(new, 'formula') & !inherits(new, cls))
        stop("the updating argument doesn't have a correct class")
    if (inherits(new, cls)){
        ncoefs <- names(coef(new))
        new <- formula(formula(new))
    }
    else ncoefs <- names(coef(update(object, new, iterlim = 0)))

    start <- numeric(length = length(ncoefs))
    names(start) <- ncoefs
    supcoef <- ! ncoefs %in% names(coef(object))
    start[names(coef(object))] <- coef(object)
    newmodel <- update(object, new, start= start, iterlim = 0)
    data.name <- paste(deparse(formula(newmodel)))
    alt.hyp <- "unconstrained model"
    stat <- - sum(newmodel$gradient * solve(newmodel$hessian, newmodel$gradient))
    names(stat) <- "chisq"
    df <- c(df = length(coef(newmodel)) - length(coef(object)))
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    result <- list(statistic = stat,
                   parameter = df,
                   p.value = pval,
                   data.name = data.name,
                   method = "score test",
                   alternative = alt.hyp
                   )
    class(result) <- 'htest'
    result
}


########################################################################
## Hausman-McFadden test for IIA
## Code adapted from the mlogit's: hmftest function and methods.
########################################################################
hmftest <- function(x, ...)
{
    UseMethod("hmftest")
}

hmftest.formula <- function(x, alt.subset, ...)
{
  formula <- x
  x <- mnlogit(formula,...)
  x$call$data <- match.call()$data
  xs <- mnlogit(formula, alt.subset=alt.subset, ...)
  hmftest(x, xs)
}

hmftest.mnlogit <- function(x, z, ...)
{
  if (is.character(z)) xs <- update(x,alt.subset=z)
  if (class(z)=="mnlogit") xs <- z
  coef.x <- coef(x)
  coef.s <- coef(xs)
  un <- names(coef.x) %in% names(coef.s)
  diff.coef <- coef.s - coef.x[un]
  diff.var <- vcov(xs) - vcov(x)[un,un]
  hmf <- as.numeric(diff.coef %*% solve(diff.var) %*% diff.coef)
  names(hmf) <- "chisq"
  df <- sum(un)
  names(df) <- "df"
  pv <- pchisq(hmf,df=df,lower.tail=FALSE)
  res <- list(data.name = x$call$data,
              statistic = hmf,
              p.value =pv,
              parameter = df,
              method = "Hausman-McFadden test",
              alternative = "IIA is rejected")
  class(res) <- "htest"
  res
}
