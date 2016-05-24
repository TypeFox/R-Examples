#' Calculate cross-validated decision curves
#'
#' This is a wrapper for 'decision_curve' that computes k-fold cross-validated estimates of sensitivity, specificity, and net benefit so that cross-validated net benefit curves can be plotted.
#'
#' @param formula an object of class 'formula' of the form outcome ~ predictors, giving the prediction model to be fitted using glm. The outcome must be a binary variable that equals '1' for cases and '0' for controls.
#' @param data data.frame containing outcome and predictors. Missing data on any of the predictors will cause the entire observation to be removed.
#' @param family a description of the error distribution and link function to pass to 'glm' used for model fitting. Defaults to binomial(link = "logit") for logistic regression.
#' @param thresholds Numeric vector of high risk thresholds to use when plotting and calculating net benefit values.
#' @param folds Number of folds for k-fold cross-validation.
#' @param study.design Either 'cohort' (default) or 'case-control' describing the study design used to obtain data. See details for more information.
#' @param population.prevalence  Outcome prevalence rate in the population used to calculate decision curves when study.design = 'case-control'.
#' @return List with components
#' \itemize{
#'   \item derived.data: derived.data: A data frame in long form showing the following for each predictor and each 'threshold', 'FPR':false positive rate, 'TPR': true positive rate, 'NB': net benefit, 'sNB': standardized net benefit, 'rho': outcome prevalence, 'prob.high.risk': percent of the population considered high risk. DP': detection probability = TPR*rho, 'model': name of prediction model or 'all' or 'none', and cost.benefit.ratio's.
#'   \item folds: number of folds used for cross-validation.
#'   \item call: matched function call.
#' }
#'
#' @seealso \code{\link{summary.decision_curve}}, \code{\link{decision_curve}},  \code{\link{Add_CostBenefit_Axis}}
#' @examples
#'
#' full.model_cv <- cv_decision_curve(Cancer~Age + Female + Smokes + Marker1 + Marker2,
#'                                   data = dcaData,
#'                                   folds = 5,
#'                                   thresholds = seq(0, .4, by = .01))
#'
#'full.model_apparent <- decision_curve(Cancer~Age + Female + Smokes + Marker1 + Marker2,
#'                                      data = dcaData,
#'                                      thresholds = seq(0, .4, by = .01),
#'                                      confidence.intervals = 'none')
#'
#'plot_decision_curve( list(full.model_apparent, full.model_cv),
#'                     curve.names = c("Apparent curve", "Cross-validated curve"),
#'                     col = c("red", "blue"),
#'                     lty = c(2,1),
#'                     lwd = c(3,2, 2, 1),
#'                     legend.position = "bottomright")
#'
#' @importFrom caret createFolds
#' @export

cv_decision_curve <- function(formula,
                           data,
                           family = binomial(link = "logit"),
                           thresholds = seq(0, 1, by = .01),
                           folds = 5,
                           study.design = c("cohort", "case-control"),
                           population.prevalence){
  call <- match.call()

  stopifnot(is.numeric(folds))
  stopifnot(folds >= 2)

  #check vars are in data
  if(any( names.check <- !is.element(all.vars(formula), names(data)))) stop(paste("variable(s)", paste( all.vars(formula)[names.check], collapse = ", ") , "not found in 'data'"))

  study.design <- match.arg(study.design)

  #throw out missing data
  data <- data[,all.vars(formula)]
  #complete case indicator
  cc.ind <- complete.cases(data)

  if(sum(cc.ind) < nrow(data)) warning(paste(sum(1-cc.ind), "observation(s) with missing data removed"))
  data <- data[cc.ind,]

  #retreive outcome
  outcome <- data[[all.vars(formula[[2]])]];
  if(length(unique(outcome)) != 2) stop("outcome variable is not binary (it does not take two unique values).")
  stopifnot(is.numeric(outcome))
  if(min(outcome) != 0 | max(outcome) != 1) stop("outcome variable must be binary taking on values 0 for control and 1 for case.")
  ####################
  ## done with checks
  ####################

  #create cross-validation folds using caret's 'createFolds'
  myfolds.ind <- createFolds(y = outcome, k = folds)
  #check to make sure there are cases and controls in each fold
  lapply(myfolds.ind, FUN = function(x){ if(length(table(outcome[x])) < 2) stop("Reduce number of folds requested: there are not enough cases to allocate across all folds")})
  #make sure there are at least 5 cases per fold.
  lapply(myfolds.ind, FUN = function(x){ if(min(table(outcome[x]))<5) stop("Reduce number of folds requested: there are not enough cases to allocate at least 5 cases into each fold.")})

  #now call `decision_curve` n = folds times and collect the results
  #call it once to allocate a spot for the results
  out <- list()
  out$derived.data <- decision_curve(formula = formula,
                        data = data,
                        fitted.risk = FALSE,
                        thresholds = thresholds,
                        confidence.intervals = "none",
                        study.design = study.design,
                        population.prevalence = population.prevalence)$derived.data
  out$derived.data[, 2:8] <- 0

  for(kk in 1:folds){

    #fit the model on -kk
    #cohort
    if(study.design == "cohort"){
      myglm <- do.call(glm, list("formula" = formula, "data" = data[-myfolds.ind[[kk]], ], "family" = family ))
      offset = 0
      }else{
      #case.control
      #offset by the relative observed outcome prevalence and the provided population rho
      obs.rho = mean(outcome)
      offset = - log((population.prevalence)/ (1-(population.prevalence))) + log((obs.rho)/(1-obs.rho))
      myglm <- do.call(glm, list("formula" = formula, "data" = data[-myfolds.ind[[kk]], ], "family" = family, "offset" = rep(offset, nrow(data[-myfolds.ind[[kk]], ])) ))

    }

    #predict on fold kk
    y <- predict(myglm, newdata = data[myfolds.ind[[kk]], ], type = "link") - offset
    y <- exp(y)/(1+exp(y))

    dat.cv <- data.frame("outcome" = outcome[myfolds.ind[[kk]]], "risk.hat" = y)

    #add measures for each fold to the first 8 columns, bc these are the
    #the only numeric estimates.
    #we divide by number of folds later
    out$derived.data[,2:8] <- out$derived.data[,2:8] +
                             decision_curve(formula = outcome~risk.hat,
                             data = dat.cv,
                             family = family,
                             fitted.risk = TRUE,
                             thresholds = thresholds,
                             confidence.intervals = "none",
                             study.design = study.design,
                             population.prevalence = population.prevalence)$derived.data[,2:8]
  }
  #take the mean across folds as estimates
  out$derived.data[,2:8] <- out$derived.data[,2:8]/folds


  #return list of elements
  out$call <- call
  out$folds <- folds
  out$confidence.intervals <- 'none'

  class(out) = "decision_curve"
  invisible(out)

}


