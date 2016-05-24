#' Build a series of linear models using automated variable selection
#' 
#' This function allows building a series of linear models (\code{lm}) using 
#' one or more automated variable selection implemented in function 
#' \code{stepVIF} and \code{stepAIC}.
#' 
#' @param formula A list containing one or several model formulas (a symbolic
#' description of the model to be fitted).
#' @param data Data frame containing the variables in the model formulas.
#' @param vif Logical for performing backward variable selection using the 
#' Variance-Inflation Factor (VIF). Defaults to \code{VIF = FALSE}.
#' @param vif.threshold Numeric value setting the maximum acceptable VIF value.
#' Defaults to \code{vif.threshold = 10}.
#' @param vif.verbose Logical for printing iteration results of backward 
#' variable selection using the VIF. Defaults to \code{vif.verbose = FALSE}.
#' @param aic Logical for performing variable selection using Akaike 
#' Information Criterion (AIC). Defaults to \code{aic = FALSE}.
#' @param aic.direction Character string setting the direction of variable 
#' selection when using AIC. Available options are \code{"both"}, 
#' \code{"forward"}, and \code{"backward"}. Defaults to 
#' \code{aic.direction = "both"}.
#' @param aic.trace Logical for printing iteration results of variable selection
#' using the AIC. Defaults to \code{aic.trace = FALSE}.
#' @param aic.steps Integer value setting the maximum number of steps to be
#' considered for variable selection using the AIC. Defaults to 
#' \code{aic.steps = 5000}.
#' @param ... Further arguments passed to the function \code{stepAIC}.
#' 
#' @details
#' This function was devised to deal with a list of linear model formulas. The 
#' main objective is to bring together several functions commonly used when
#' building linear models, such as automated variable selection. In the current 
#' implementation, variable selection can be done using \code{stepVIF} or 
#' \code{stepAIC} or both. \code{stepVIF} is a backward variable selection
#' procedure, while \code{stepAIC} supports backward, forward, and bidirectional
#' variable selection. For more information about these functions, please visit 
#' their respective help pages.
#' 
#' An important feature of \code{buildMS} is that it records the initial number
#' of candidate predictor variables and observations offered to the model, and 
#' adds this information as an attribute to the final selected model. Such 
#' feature was included because variable selection procedures result biased 
#' linear models (too optimistic), and the effective number of degrees of 
#' freedom is close to the number of candidate predictor variables initially
#' offered to the model (Harrell, 2001). With the initial number of candidate
#' predictor variables and observations offered to the model, one can calculate
#' penalized or adjusted measures of model performance. For models built using
#' \code{builtMS}, this can be done using \code{statsMS}.
#' 
#' Some important details should be clear when using \code{buildMS}:
#' 
#' \enumerate{
#' \item this function was originally devised to deal with a list of formulas, 
#' but can also be used with a single formula;
#' \item in the current implementation, \code{stepVIF} runs before 
#' \code{stepAIC};
#' \item function arguments imported from \code{stepAIC} and \code{stepVIF} 
#' were named as in the original functions, and received a prefix (\code{aic} 
#' or \code{vif}) to help the user identifying which function is affected by a 
#' given argument without having to go check the documentation.
#' }
#' 
#' @return A list containing the fitted linear models.
#' 
#' @references
#' Harrell, F. E. (2001) \emph{Regression modelling strategies: with 
#' applications to linear models, logistic regression, and survival analysis.}
#' First edition. New York: Springer.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern applied statistics 
#' with S.} Fourth edition. New York: Springer.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' 
#' @section TODO: Add option to set the order in which \code{stepAIC} and
#' \code{stepVIF} are run.
#' 
#' @seealso \code{\link[MASS]{stepAIC}}, \code{\link[pedometrics]{stepVIF}},
#' \code{\link[pedometrics]{statsMS}}.
#' @export
#' @examples
#' \dontrun{
#' # based on the second example of function stepAIC
#' require(MASS)
#' cpus1 <- cpus
#' for(v in names(cpus)[2:7])
#'   cpus1[[v]] <- cut(cpus[[v]], unique(stats::quantile(cpus[[v]])),
#'                     include.lowest = TRUE)
#' cpus0 <- cpus1[, 2:8]  # excludes names, authors' predictions
#' cpus.samp <- sample(1:209, 100)
#' cpus.form <- list(formula(log10(perf) ~ syct + mmin + mmax + cach + chmin +
#'                   chmax + perf),
#'                   formula(log10(perf) ~ syct + mmin + cach + chmin + chmax),
#'                   formula(log10(perf) ~ mmax + cach + chmin + chmax + perf))
#' data <- cpus1[cpus.samp,2:8]
#' cpus.ms <- buildMS(cpus.form, data, vif = TRUE, aic = TRUE)
#' }
#' @keywords iteration models
#' 
# FUNCTION #####################################################################
buildMS <- 
  function (formula, data,
            vif = FALSE, vif.threshold = 10, vif.verbose = FALSE,
            aic = FALSE, aic.direction = "both", aic.trace = FALSE,
            aic.steps = 5000, ...) {
    
    # Check if suggested packages are installed
    pkg <- c("MASS", "pbapply")
    id <- !sapply(pkg, requireNamespace, quietly = TRUE)
    if (any(id)) {
      pkg <- paste(pkg[which(id)], collapse = " ")
      stop(paste("Package(s) needed for this function to work but not",
                 "installed: ", pkg, sep = ""), call. = FALSE)
    }
    
    # check arguments ##########################################################
    if (missing(formula)) {
      stop("<formula> is a mandatory argument")
    }
    if (class(formula) != "list") {
      formula <- list(formula)
    }
    if (missing(data)) {
      stop("<data> is a mandatory argument")
    }
    if (class(data) != "data.frame") {
      data <- as.data.frame(data)
    }
    
    # lm()
    print("fitting linear model using ols")
    model <- pbapply::pblapply(formula, function (X) stats::lm(X, data))
    # get the initial number of candidate predictors and observations
    p <- sapply(model, function (X) dim(stats::model.matrix(X))[2])
    n <- sapply(model, function (X) dim(stats::model.matrix(X))[1])
    
    # stepVIF()
    if (vif) {
      print("backward variable selection using VIF")
      model <- pbapply::pblapply(model, function (X)
        stepVIF(X, threshold = vif.threshold, verbose = vif.verbose))
    }
    
    # stepAIC()
    if (aic) {
      print(paste(aic.direction, " variable selection using AIC", sep = ""))
      model <- pbapply::pblapply(model, function (X) 
        MASS::stepAIC(X, direction = aic.direction, steps = aic.steps, 
                      trace = aic.trace, ...))
    }
    
    # Prepare output - add attributes to the final model
    a <- lapply(model, attributes)
    for (i in 1:length(a)) {
      a[[i]]$p <- p[i]
    }
    for (i in 1:length(a)) {
      a[[i]]$n <- n[i]
    }
    for (i in 1:length(model)) {
      attributes(model[[i]]) <- a[[i]]
    }
    return(model)
  }
# End!
