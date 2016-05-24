##' A package for estimating strategic statistical models.
##' 
##' @name games-package
##' @docType package
##' @section Acknowledgements: We thank the Wallis Institute of Political
##' Economy for financial support.
##' @references
##' Brenton Kenkel and Curtis S. Signorino.  2014.  "Estimating Extensive Form
##' Games in R."  \emph{Journal of Statistical Software} 56(8):1--27.
##' 
##' Curtis S. Signorino.  2003.  "Structure and Uncertainty
##' in Discrete Choice Models."  \emph{Political Analysis} 11:316--344.
##' @import stringr
##' @import Formula
##' @import maxLik
##' @import MASS
##' @export predict.egame12
##' @export predict.egame122
##' @export predict.egame123
##' @export predict.ultimatum
##' @export clarke
NULL

##' Print a strategic model object
##' 
##' The default method for printing a \code{game} object.
##'
##' Prints the call and coefficients of a fitted strategic model.
##' @param x a fitted model of class \code{game}
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
print.game <- function(x, ...)
{
    oldx <- x
    cat("\nA fitted strategic model\n\nCALL:\n\n")
    print(x$call)
    cat("\nCOEFFICIENTS:\n")

    ## print a table of coefficients for each utility (or variance) equation
    for (i in seq_along(x$equations)) {
        eq <- x$equations[i]
        hc <- attr(x$equations, "hasColon")[i]

        ## finds which coefficients contain the name of the relevant equation.
        ## this is a hack; in the future, fitted models should contain a matrix
        ## where each row is a regressor name and each column is an equation
        ## name in order to avoid this
        cf <- grep(eq, names(x$coefficients), fixed = TRUE)

        cat("\n  ", prefixToString(eq), "\n", sep = "")
        if (length(cf) > 0) {
            isFixed <- all(x$fixed[cf])
            cf <- x$coefficients[cf]

            if (hc) {
                ## this strips out the equation prefix in each term; e.g.,
                ## "u1(war):x1" becomes "x1"
                names(cf) <- sapply(strsplit(names(cf),
                                             paste(eq, ":", sep = ""),
                                             fixed = TRUE), "[", -1)
            } else {  ## i.e., the term is estimated itself, without regressors
                names(cf) <- if (isFixed) "fixed to" else "estimated as"
            }
            
            names(cf) <- paste("     ", names(cf), sep = "")
            cf <- data.frame(as.matrix(cf))
            names(cf) <- " "
            print(cf)
        } else {
            ## this is for cases when there is a utility equation with no terms
            ## estimated
            cat("\n     fixed to 0\n")
        }
    }

    cc <- convergenceCriterion(x$convergence$method)
    if (!(x$convergence$code %in% cc)) {
        cat("\nWarning: Model fitting did not converge\nCode:",
            x$convergence$code, "\nMessage:", x$convergence$message)
    }

    if (!x$localID)
        warning("Hessian is not negative definite; estimate may not be a strict local maximum")

    cat("\n")
    invisible(oldx)
}

##' Summarize a strategic model object
##' 
##' The default method for summarizing a \code{game} object.
##'
##' Forms the standard regression results table from a fitted strategic model.
##' Normally used interactively, in conjunction with
##' \code{\link{print.summary.game}}.
##' @param object a fitted model of class \code{game}
##' @param useboot logical: use bootstrap estimates (if present) to construct
##' standard error estimates?
##' @param ... other arguments, currently ignored
##' @return an object of class \code{summary.game}, containing the coefficient
##' matrix and other information needed for printing
##' @seealso \code{\link{print.summary.game}}
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
summary.game <- function(object, useboot = TRUE, ...)
{
    useboot <- useboot && !is.null(object$boot.matrix)
    if (useboot)
        object$vcov <- var(object$boot.matrix)

    ## makes the standard "regressor table" (as in summary.lm)
    cf <- object$coefficients[!object$fixed]
    se <- sqrt(diag(object$vcov[!object$fixed, !object$fixed, drop = FALSE]))
    zval <- cf / se
    pval <- 2 * pnorm(-abs(zval))

    ans <- list()
    ans$coefficients <- cbind(cf, se, zval, pval)
    colnames(ans$coefficients) <- c("Estimate", "Std. Error", "z value",
                                    "Pr(>|z|)")
    ans$call <- object$call
    ans$log.likelihood <- sum(object$log.likelihood)
    ans$nobs <- nrow(object$model)
    ans$fixed.terms <- object$coefficients[object$fixed]
    ans$convergence <- object$convergence
    ans$useboot <- useboot
    ans$localID <- object$localID
    class(ans) <- "summary.game"

    return(ans)
}

##' Print strategic model summary
##' 
##' Print output from \code{summary.game}
##'
##' Prints the standard regression results table from a fitted strategic model,
##' along with the log-likelihood, AIC, and number of observations.
##' @method print summary.game
##' @param x an object of class \code{summary.game}, typically produced by
##'running \code{summary} on a fitted model of class \code{game}
##' @param ... other arguments, currently ignored
##' @return \code{x}, invisibly.
##' @seealso \code{\link{summary.game}}
##' @export
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
print.summary.game <- function(x, ...)
{
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients)
    if (x$useboot) {
        cat("\nStandard errors estimated from bootstrap results\n")
    } else {
        cat("\nStandard errors estimated from inverse Hessian\n")
    }
    if (length(x$fixed.terms)) {
        cat("\nFixed terms:\n")
        print(x$fixed)
    }
    cat("\nLog-likelihood:", x$log.likelihood)
    cat("\nAIC:", AIC(x))
    cat("\nNo. observations:", x$nobs, "\n\n")
    cc <- convergenceCriterion(x$convergence$method)
    if (!(x$convergence$code %in% cc)) {
        cat("\nWarning: Model fitting did not converge\nCode:",
            x$convergence$code, "\nMessage:", x$convergence$message)
    }
    if (!x$localID)
        warning("Hessian is not negative definite; coefficients may not be a strict local maximum")
    invisible(x)
}

##' @export
coef.game <- function(object, ...)
{
    object$coefficients
}

##' @export
vcov.game <- function(object, ...)
{
    object$vcov
}

##' @export
logLik.game <- function(object, ...)
{
    ans <- sum(object$log.likelihood)
    attr(ans, "df") <- length(object$coefficients) - sum(object$fixed)
    attr(ans, "nobs") <- nrow(object$model)
    class(ans) <- "logLik"
    return(ans)
}

##' @method logLik summary.game
##' @export
logLik.summary.game <- function(object, ...)
{
    ans <- object$log.likelihood
    attr(ans, "df") <- nrow(object$coefficients)
    attr(ans, "nobs") <- object$nobs
    class(ans) <- "logLik"
    return(ans)
}

##' Predicted probabilities for strategic models
##' 
##' Makes predicted probabilities from a strategic model.
##'
##' This method uses a fitted strategic model to make predictions for a new
##' set of data.  This is useful for cross-validating or for graphical
##' analysis.  For many uses, such as analyzing the marginal effect of a
##' particular independent variable, the function \code{\link{predProbs}} will
##' be more convenient.
##'
##' In the \code{\link{ultimatum}} model, there is not an analytic expression
##' for the expected value of Player 1's offer.  Therefore, predicted values
##' are instead generating via simulation by drawing errors from a logistic
##' distribution.  The number of draws per observation can be controlled via
##' the \code{n.sim} argument.  For replicability, we recommend seeding the
##' random number generator via \code{\link{set.seed}} before using
##' \code{predict.ultimatum}.
##' @aliases predict.game predict.egame12 predict.egame122 predict.egame123
##' predict.ultimatum
##' @usage
##' \method{predict}{game}(object, ...)
##'
##' \method{predict}{egame12}(object, newdata, type=c("outcome", "action"), na.action = na.pass, ...)
##' \method{predict}{egame122}(object, newdata, type=c("outcome", "action"), na.action = na.pass, ...)
##' \method{predict}{egame123}(object, newdata, type=c("outcome", "action"), na.action = na.pass, ...)
##' \method{predict}{ultimatum}(object, newdata, na.action = na.pass, n.sim = 1000, ...)
##' @param object a fitted model of class \code{game}.
##' @param newdata data frame of values to make the predicted probabilities for.
##' If this is left empty, the original dataset is used.
##' @param type whether to provide probabilities for outcomes (e.g., L, RL, or
##' RR in \code{egame12}) or for actions (e.g., whether 2 moves L or R given
##' that 1 moved R).
##' @param na.action how to deal with \code{NA}s in \code{newdata}
##' @param n.sim number of simulation draws to use per observation for
##' \code{ultimatum} models (see Details).
##' @param ... other arguments, currently ignored.
##' @return A data frame of predicted probabilities.
##' @export
##' @seealso \code{\link{predProbs}} provides a more full-featured and
##' user-friendly wrapper, including plots and confidence bands.
##' @author Brenton Kenkel (\email{brenton.kenkel@@gmail.com})
predict.game <- function(object, ...)
{
    NextMethod("predict", object, ...)
}

##
## INPUT:
## x: character string
##
## RETURN:
## plain English version of 'x' (e.g., "u1(war)" becomes "Player 1's utility for
## war")
## 
prefixToString <- function(x)
{
    first <- substr(x, 1, 1)  # first character
    
    if (first == "u") {  # utility equation
        ## this is for the unlikely event of a 10+ player model, i.e., to ensure
        ## that "u10(war)" isn't translated into "Player 1's utility for (war"
        n <- 3
        while (substr(x, n, n) != "(") n <- n + 1

        player <- substr(x, 2, n - 1)  # player number
        outcome <- substr(x, n + 1, nchar(x) - 1)  # outcome name

        x <- paste("Player ", player, "'s utility for ", outcome, ":", sep = "")
    } else if (first == "l") {  # variance equation
        ## this one is not yet safe in case of a 10+ player model, just takes
        ## the second-last character (e.g., "1" in "log(sigma1)") and checks
        ## whether it's numeric.  if so, this is taken to be the player number;
        ## if not, the scale parameter must be common to all players
        player <- substr(x, nchar(x) - 1, nchar(x) - 1)
        player <- tryCatch(as.numeric(player), warning = identity)
        if (!inherits(player, "warning")) {
            x <- paste("Logged scale parameter for player ", player, ":", sep
                       = "")
        } else {  # a warning means the second-last character couldn't be
                  # coerced to numeric, so it's not a number
            x <- "Logged scale parameter:"
        }
    } else if (first == "R") {  # reservation value in ultimatum model
        x <- paste("Player ", substr(x, 2, nchar(x)), "'s reservation value:",
                   sep = "")
    }

    return(x)
}
