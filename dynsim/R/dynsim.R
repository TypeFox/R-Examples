#' Dynamic simulations of autoregressive relationships
#'
#' \code{dynsim} dynamic simulations of autoregressive relationships
#'
#' @param obj the output object the estimation model.
#' @param ldv character. Names the lagged dependent variable
#' @param scen data frame or list of data frames. Specifies the values of the
#' variables used to generate the predicted values when \eqn{t = 0}. If only one
#' scenario is desired then \code{scen} should be a data frame. If more than one
#' scenario is desired then the \eqn{t = 0} values should be in data frames
#' contained in a list.
#' @param n numeric. Specifies the number of iterations (or time period) over
#' which the program will generate the predicted value of the dependent
#' variable. The default is 10.
#' @param sig numeric. Specifies the level of statistical significance of the
#' confidence intervals. Any value allowed be greater than 0 and cannot be
#' greater than 1.
#' @param num numeric. Specifies the number of simulations to compute for each
#' value of \code{n}. The default is 1000.
#' @param shocks data frame. Allows the user to choose independent variables,
#' their values, and times to introduce these values. The first column of the
#' data frame must be called \code{times} this will contain the times in
#' \code{n} to use the shock values. The following columns' names must match the
#' names of the variables whose values you wish to alter. You do not need to
#' specify values for variables that you want to remain the same as in
#' \code{scen}. In times \code{n} where shock values are not specified,
#' non-\code{ldv} variable values will revert to those in \code{scen}. If
#' \code{*} is used to create interactions, interaction terms will be fitted
#' appropriately.
#' @param ... arguments to pass to methods.
#'
#' @details A post-estimation technique for producing dynamic simulations of
#' autoregressive models.
#'
#' @return The command returns a \code{data.frame} and \code{dynsim} class
#' object. This can contain up to columns elements:
#' \itemize{
#'  \item{\code{scenNumber}: }{The scenario number.}
#'  \item{\code{time}: }{The time points.}
#'  \item{\code{shock.}: }{Columns containing the values of the shock variables
#'  at each point in \code{time}.}
#'  \item{\code{ldvMean}: }{Mean of the simulation distribution.}
#'  \item{\code{ldvLower}: }{Lower bound of the simulation distribution's
#'  central interval set with \code{sig}.}
#'  \item{\code{ldvUpper}: }{Upper bound of the simulation distribution's
#'  central interval set with \code{sig}.}
#'  \item{\code{ldvLower50}: }{Lower bound of the simulation distribution's
#'  central 50 percent interval.}
#'  \item{\code{ldvUpper50}: }{Upper bound of the simulation distribution's
#'  central 50 percent interval.}
#' }
#' The output object is a data frame class object. Do with it as you like.
#'
#' @examples
#' # Load package
#' library(DataCombine)
#'
#' # Load Grunfeld data
#' data(grunfeld, package = "dynsim")
#'
#' # Create lag invest variable
#' grunfeld <- slide(grunfeld, Var = "invest", GroupVar = "company",
#'                NewVar = "InvestLag")
#'
#' # Convert company to factor for fixed-effects specification
#' grunfeld$company <- as.factor(grunfeld$company)
#'
#' # Estimate basic model
#' M1 <- lm(invest ~ InvestLag + mvalue + kstock + company, data = grunfeld)
#'
#' # Estimate model with interaction between mvalue and kstock
#' M2 <- lm(invest ~ InvestLag + mvalue*kstock + company, data = grunfeld)
#'
#' # Set up scenarios for company 4
#' ## List version ##
#' attach(grunfeld)
#' Scen1 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
#'                     mvalue = quantile(mvalue, 0.05),
#'                     kstock = quantile(kstock, 0.05),
#'                     company4 = 1)
#' Scen2 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
#'                     mvalue = mean(mvalue),
#'                     kstock = mean(kstock),
#'                     company4 = 1)
#' Scen3 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
#'                     mvalue = quantile(mvalue, 0.95),
#'                     kstock = quantile(kstock, 0.95),
#'                     company4 = 1)
#' detach(grunfeld)
#'
#' \dontrun{
#' ## Alternative data frame version of the scenario builder ##
#' attach(grunfeld)
#' ScenComb <- data.frame(InvestLag = rep(mean(InvestLag, na.rm = TRUE), 3),
#'                       mvalue = c(quantile(mvalue, 0.95), mean(mvalue), 
#'                                  quantile(mvalue, 0.05)),
#'                       kstock = c(quantile(kstock, 0.95), mean(kstock),
#'                                  quantile(kstock, 0.05)),
#'                       company4 = rep(1, 3)
#' )
#' detach(grunfeld)
#' }
#' 
#' # Combine into a single list
#' ScenComb <- list(Scen1, Scen2, Scen3)
#'
#' ## Run dynamic simulations without shocks and no interactions
#' Sim1 <- dynsim(obj = M1, ldv = "InvestLag", scen = ScenComb, n = 20)
#'
#' ## Run dynamic simulations without shocks and interactions
#' Sim2 <- dynsim(obj = M2, ldv = "InvestLag", scen = ScenComb, n = 20)
#'
#' ## Run dynamic simulations with shocks
#'
#' # Create data frame of shock values
#' mShocks <- data.frame(times = c(5, 10), kstock = c(100, 1000),
#'                       mvalue = c(58, 5000))
#'
#' # Run simulations without interactions
#' Sim3 <- dynsim(obj = M1, ldv = "InvestLag", scen = ScenComb, n = 20,
#'                shocks = mShocks)
#'
#' # Run simulations with interactions
#' Sim4 <- dynsim(obj = M2, ldv = "InvestLag", scen = ScenComb, n = 20,
#'                shocks = mShocks)
#'
#' @references
#' Williams, L. K., & Whitten, G. D. (2011). Dynamic Simulations of
#' Autoregressive Relationships. The Stata Journal, 11(4), 577-588.
#'
#' Williams, L. K., & Whitten, G. D. (2012). But Wait, There's More! Maximizing
#' Substantive Inferences from TSCS Models. Journal of Politics, 74(03),
#' 685-693.
#' 
#' @importFrom stats coef
#'
#' @export

dynsim <- function(obj, ldv, scen, n = 10, sig = 0.95, num = 1000,
                   shocks = NULL, ...) {
    # Zelig no longer used
    if ('zelig' %in% class(obj)) {
        stop(paste0('dynsim no longer relies on Zelig.\n',
            '---- Please use `lm`. ----'),
            call. = FALSE)
    }

    ModCoefNames <- names(coef(obj))

    # Create mean fitted values if scen is not specified
    if (missing(scen)) {
        stop('\nNo scen provided. Please specify scenario values.',
            call. = FALSE)
    }

    # Make sure that the variables in scen are in the model
    VarMisError <- '\nAt least one variable name in scen was not found in the estimation model.'
    if (is.data.frame(scen)) {
        if (any(!(names(scen) %in% ModCoefNames))) {
            stop(VarMisError, call. = FALSE)
        }
    } else if (is.list(scen)) {
        for (a in seq_along(scen)) {
            TempDF <- scen[[a]]
            if (any(!(names(TempDF) %in% ModCoefNames))) {
                stop(VarMisError, call. = FALSE)
            }
        }
    }

    # Make sure that both shocks is a data frame and the first column of shocks
    # is a variable called "times".
    if (!is.null(shocks)) {
        if (!is.data.frame(shocks)) {
            stop("\nShocks must be a data frame.", call. = FALSE)
        }
        if (names(shocks)[1] != "times") {
            stop("\nThe first variable of shocks must be called 'times' and contain the shock times.",
               call. = FALSE)
        }
    }
    # Error if number of iterations is <= 0.
    if (n <= 0) {
        stop("\nYou must specify at least 1 iteration with the n argument.",
             call. = FALSE)
    }
    # Make sure sig is between 0 and 1.
    if (sig <= 0 | sig > 1) {
        stop("\nsig must be greater than 0 and not greater than 1.",
             call. = FALSE)
    }

    # Determine if 1 or more scenarios are desired and simulate scenarios
    if (is.data.frame(scen)) {
        if (nrow(scen) == 1) {
            SimOut <- OneScen(obj = obj, ldv = ldv, n = n, num = num,
                                scen = scen, sig = sig, shocks = shocks)
        }
        else if (nrow(scen > 1)) {
            results <- list()
            for (u in 1:nrow(scen)) {
                scen_sub <- scen[u, ]
                SimTemp <- OneScen(obj = obj, ldv = ldv, n = n,
                                    scen = scen_sub, sig = sig, num = num,
                                    shocks = shocks)
                results[[u]] <- cbind("scenNumber" = u, SimTemp)
            }
            SimOut <- do.call("rbind", results)
        }
    }
    else if (is.list(scen)) {
        results <- list()
        for (u in seq_along(scen)) {
            SimTemp <- OneScen(obj = obj, ldv = ldv, n = n,
                                scen = scen[[u]], sig = sig, num = num,
                                shocks = shocks)
            results[[u]] <- cbind("scenNumber" = u, SimTemp)
        }
        SimOut <- do.call("rbind", results)
    }

    # Ascribe class and return output
    class(SimOut) <- c("data.frame", "dynsim")
    return(SimOut)
}
