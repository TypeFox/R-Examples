#' @title
#' Estimate Seroincidences
#'
#' @description
#' Function to estimate seroincidences based on cross-section serology data and
#' longitudinal response model.
#'
#' @param
#' data Dataframe with cross-sectional serum antibody data. At least one
#' antibody should be specified. Rows represent multiple antibody measurements
#' in the same serum sample. Strata may be specified in additional columns.
#' @param
#' antibodies Character vector with names of antibodies to be included
#' for incidence calculations. At least one antibody name must be specified,
#' names must match with those in \code{data} and in \code{Ak}.
#' @param
#' strata Character vector with names of strata to be used, default ="". Names
#' must match with data.
#' @param
#' Ak List with two dataframes containing longitudinal parameters named
#' A (the peak serum antibody level) and k (the decay rate). Each dataframe
#' contains sets of Monte Carlo samples, one for each antibody that is available.
#' @param
#' censorLimits List of cutoffs for each of the antibodies specified
#' in 'antibodies'. Any cross-sectional observations below cutoff will be treated
#' as censored at the cutoff level.
#' @param
#' showProgress Indicator whether or not to show text progress bar indicating progress of estimation.
#' The setting is effective only if the number of levels per stratum if greater than 1. Default is \code{FALSE}.
#' @param
#' ... Additional arguments passed to function \code{\link{optim}}.
#'
#' @return
#' A list with the following items:
#' \describe{
#' \item{\code{Fits}}{Outputs of \code{\link{optim}} function per stratum.}
#' \item{\code{Antibodies}}{Input parameter \code{antibodies} saved for \code{\link{summary}} function.}
#' \item{\code{Strata}}{Input parameter \code{strata} saved for \code{\link{summary}} function.}
#' \item{\code{CensorLimits}}{Input parameter \code{censorLimits} saved for \code{\link{summary}} function.}
#' }
#'
#' @examples
#'
#' \dontrun{
#' estimateSeroincidence(data = serologyData, antibodies = c("IgG", "IgM", "IgA"),
#'                          strata = c("age"), Ak = responseParams)
#' }
#'
#' @export
estimateSeroincidence <- function(data, antibodies, strata = "", Ak, censorLimits, showProgress = FALSE, ...)
{
    # Check antibodies
    if (!is.character(antibodies))
        stop("Argument \"antibodies\" is not a character vector.\nProvide a character vector with at least one antibody name")

    if (all(antibodies == ""))
        stop("Argument \"antibodies\" is empty.\nProvide a character vector with at least one antibody name")

    # Check data
    if (!is.data.frame(data))
        stop("Argument \"data\" is not a dataframe.\nProvide a dataframe with cross-sectional serology data per antibody")

    if (!all(is.element(antibodies, names(data))))
        stop("Antibody names in argument \"data\" and argument \"antibodies\" do not match")

    # Check strata
    if (!is.character(strata))
        stop("Argument \"strata\" is not a character vector.\nProvide a character vector with strata names")

    if (!all(is.element(strata, c("", names(data)))))
        stop("Strata names in argument \"data\" and argument \"strata\" do not match")

    # Check AkSim
    if (!is.list(Ak) | length(Ak) != 2)
        stop("Argument \"Ak\" is not a list with two dataframes.\nProvide a list with two dataframes named A and k\nEach dataframe contains Monte Carlo simulations per antibody")

    if (!is.data.frame(Ak[[1]]) | !is.data.frame(Ak[[2]]))
        stop("Argument \"Ak\" does not contain two dataframes.\nProvide a list with two dataframes named A and k\nEach dataframe contains Monte Carlo simulations per antibody")

    if (!all(is.element(names(Ak), c("A", "k"))))
        stop("The two dataframes in argument \"Ak\" are not named A and k.\nProvide a list with two dataframes named A and k\nEach dataframe contains Monte Carlo simulations per antibody")

    if (nrow(Ak$A) != nrow(Ak$k))
        stop("The two dataframes A and k in argument \"Ak\" are of different length")

    if (!all(is.element(antibodies, names(Ak$A))))
        stop("Antibody names in argument \"antibodies\" and argument \"Ak\" do not match")

    # Make stratum variable (if needed)
    if (all(strata != ""))
        data <- within(data, stratum <- interaction(data[, strata]))
    else
        data <- within(data, stratum <- factor(1))

    levelsStrata <- levels(data$stratum)
    numLevels <- nlevels(data$stratum)

    # Loop over levelsStrata
    progressBarCreated <- FALSE
    if (showProgress & numLevels > 1)
    {
        pb <- utils::txtProgressBar(min = 0, max = numLevels)
        progressBarCreated <- TRUE
    }

    fits <- list()
    for (i in seq_len(nlevels(data$stratum)))
    {

        # Make subset of data: per stratum, select antibodies
        y <- subset(data, subset = stratum == levelsStrata[i], select = antibodies)

        # Incidence can not be calculated if there are zero observations
        if (nrow(y) == 0) next

        # Estimate log.lambda, starting value = log(1/365.25) day^-1
        fit <- stats::optim(par = log(1 / 365.25), fn = .nll, y = y, m = 0, Ak = Ak, censorLimits = censorLimits, method = "BFGS", hessian = TRUE, ...)

        # Collect results
        fits[[levelsStrata[i]]] <- fit

        if (progressBarCreated)
            utils::setTxtProgressBar(pb, i)

    }
    if (progressBarCreated)
        close(pb)

    incidenceData <- list(Fits = fits, Antibodies = antibodies, Strata = strata, CensorLimits = censorLimits)

    class(incidenceData) <- c("seroincidence", "list")

    # Return seroincidence object
    return(incidenceData)
}
