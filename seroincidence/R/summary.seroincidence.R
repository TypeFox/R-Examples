#' @title
#' Summary Method for Seroincidence Object
#'
#' @description
#' Calculate seroincidence from output of the seroincidence calculator \code{\link{estimateSeroincidence}}.
#'
#' @param
#' object A dataframe containing output of function \code{\link{estimateSeroincidence}}.
#' @param
#' ... Additional arguments affecting the summary produced.
#' @param
#' quantiles A vector of length 2 specifying quantiles for lower (first element) and upper (second element) bounds of \code{lambda}. Default \code{c(0.025, 0.975)}.
#' @param
#' showDeviance Logical flag (\code{FALSE}/\code{TRUE}) for reporting deviance (-2*log(likelihood) at estimated seroincidence. Default is \code{TRUE}.
#' @param
#' showConvergence Logical flag (\code{FALSE}/\code{TRUE}) for reporting convergence (see help for \code{\link{optim}} for details). Default is \code{TRUE}.
#'
#' @return
#' A list with the following items:
#' \describe{
#' \item{\code{Results}}{Dataframe with maximum likelihood estimate of \code{lambda} (the seroincidence) (column \code{Lambda}) and
#'                          corresponding lower (\code{Lambda.lwr}) and upper (\code{Lambda.upr} bounds.
#'                          Optionally \code{Deviance} (Negative log likelihood (NLL) at estimated (maximum likelhood) \code{lambda})
#'                          and \code{Covergence} (Convergence indicator returned by \code{\link{optim}}. Value of 0 indicates convergence)
#'                          columns are included.}
#' \item{\code{Antibodies}}{Character vector with names of input antibodies used in \code{\link{estimateSeroincidence}}.}
#' \item{\code{Strata}}{Character with names of strata used in \code{\link{estimateSeroincidence}}.}
#' \item{\code{CensorLimits}}{List of cutoffs for each of the antibodies used in \code{\link{estimateSeroincidence}}.}
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' # estimate seroincidence
#' seroincidence <- estimateSeroincidence(...)
#'
#' # calculate summary statistics for the seroincidence object
#' seroincidenceSummary <- summary(seroincidence)
#' }
#'
#' @export
summary.seroincidence <- function(object, ..., quantiles = c(0.025, 0.975), showDeviance = TRUE, showConvergence = TRUE)
{

    if (length(quantiles) != 2 | any(quantiles < 0) | any(quantiles > 1))
        stop("Incorrectly specified quantiles")

    if (quantiles[1] > quantiles[2])
        stop("Quantile for upper bound of incidence estimate cannot be less than the lower bound")

    fits <- object[["Fits"]]
    results <- as.data.frame(t(sapply(fits, FUN = function(elem)
        {
            with(elem, c(
                Lambda = 365.25 * exp(par + qnorm(0.5) * sqrt(1 / hessian)),
                Lambda.lwr = 365.25 * exp(par + qnorm(quantiles[1]) * sqrt(1 / hessian)),
                Lambda.upr = 365.25 * exp(par + qnorm(quantiles[2]) * sqrt(1 / hessian)),
                Deviance = 2 * value,
                Convergence = convergence))
        })))
    results$Stratum <- rownames(results)
    rownames(results) <- NULL

    if (!showDeviance)
        results$Deviance <- NULL

    if (!showConvergence)
        results$Convergence <- NULL

    output <- list(Results = results,
                   Antibodies = object[["Antibodies"]],
                   Strata = object[["Strata"]],
                   CensorLimits = object[["CensorLimits"]],
                   Quantiles = quantiles)

    class(output) <- c("summary.seroincidence", "list")

    return(output)
}
