#' Calculate decision curves
#'
#'Decision curves are a useful tool to evaluate the population impact of adopting a risk prediction instrument into clinical practice. Given one or more instruments (risk models) that predict the probability of a binary outcome, this package calculates and plots decision curves, which display estimates of the standardized net benefit by the probabilty threshold used to categorize observations as 'high risk.'  Bootstrap confidence intervals are displayed as well. This package is a companion to the manuscript '(put ref here.)'.
#'
#' @param formula an object of class 'formula' of the form outcome ~ predictors, giving the prediction model to be fitted using glm. The outcome must be a binary variable that equals '1' for cases and '0' for controls.
#' @param data data.frame containing outcome and predictors. Missing data on any of the predictors will cause the entire observation to be removed.
#' @param family a description of the error distribution and link function to pass to 'glm' used for model fitting. Defaults to binomial(link = "logit") for logistic regression.
#' @param fitted.risk logical (default FALSE) indicating whether the predictor provided are estimated risks from an already established model. If set to TRUE, no model fitting will be done and all estimates will be conditional on the risks provided.  Risks must fall between 0 and 1.
#' @param thresholds Numeric vector of high risk thresholds to use when plotting and calculating net benefit values.
#' @param confidence.intervals Numeric (default 0.95 for 95\% confidence bands) level of bootstrap confidence intervals to plot. Set as NA or 'none' to remove confidence intervals. See details for more information.
#' @param bootstraps Number of bootstrap replicates to use to calculate confidence intervals (default 500).
#' @param study.design Either 'cohort' (default) or 'case-control' describing the study design used to obtain data. See details for more information.
#' @param population.prevalence  Outcome prevalence rate in the population used to calculate decision curves when study.design = 'case-control'.
#' @details  Confidence intervals for (standardized) net benefit are calculated pointwise at each risk threshold. For when data come from an observational cohort, bootstrap sampling is done without stratifying on outcome, so disease prevalence varies within bootstrap samples. For case-control data, bootstrap sampling is done stratified on outcome.
#' @return List with components
#' \itemize{
#'   \item derived.data: A data frame in long form showing the following for each predictor and each 'threshold', 'FPR':false positive rate, 'TPR': true positive rate, 'NB': net benefit, 'sNB': standardized net benefit, 'rho': outcome prevalence, 'prob.high.risk': percent of the population considered high risk. DP': detection probability = TPR*rho, 'model': name of prediction model or 'all' or 'none', cost.benefit.ratio, and 'xx_lower', 'xx_upper': the lower and upper confidence bands for all measures (if calculated).
#'   \item confidence.intervals: Level of confidence intervals returned.
#'   \item call: matched function call.
#' }
#'
#' @seealso \code{\link{summary.decision_curve}}, \code{\link{cv_decision_curve}}, \code{\link{Add_CostBenefit_Axis}}
#' @examples
#'#helper function
#' expit <- function(xx) exp(xx)/ (1+exp(xx))
#'
#'#load simulated cohort data
#'data(dcaData)
#'baseline.model <- decision_curve(Cancer~Age + Female + Smokes,
#'                                 data = dcaData,
#'                                 thresholds = seq(0, .4, by = .01),
#'                                 study.design = "cohort",
#'                                 bootstraps = 10) #number of bootstraps should be higher
#'
#'full.model <- decision_curve(Cancer~Age + Female + Smokes + Marker1 + Marker2,
#'                             data = dcaData,
#'                             thresholds = seq(0, .4, by = .01),
#'                             bootstraps = 10)
#'
#'#simulated case-control data with same variables as above
#'data(dcaData_cc)
#'
#'table(dcaData_cc$Cancer)
#'
#'#estimated from the population where the
#'#case-control sample comes from.
#'population.rho = 0.11
#'
#'full.model_cc <- decision_curve(Cancer~Age + Female + Smokes + Marker1 + Marker2,
#'                                data = dcaData,
#'                                thresholds = seq(0, .4, by = .01),
#'                                bootstraps = 10,
#'                                study.design = "case-control",
#'                                population.prevalence = population.rho)
#'
#'
#' @import MASS
#' @importFrom grDevices rainbow
#' @importFrom graphics axis grid legend lines mtext par  plot
#' @importFrom stats binomial complete.cases fitted glm predict quantile
#' @export

decision_curve <- function(formula,
                          data,
                          family = binomial(link = "logit"),
                          fitted.risk = FALSE,
                          thresholds = seq(0, 1, by = .01),
                          confidence.intervals = 0.95,
                          bootstraps = 500,
                          study.design = c("cohort", "case-control"),
                          population.prevalence){
  call <- match.call()

  #################
  ### Begin Checks
  stopifnot(class(formula) == "formula")  #check formula
  stopifnot(is.data.frame(data)) #check data
  stopifnot(is.logical(fitted.risk))
  stopifnot(is.numeric(thresholds))
  stopifnot(all(thresholds >= 0)); stopifnot(all(thresholds <= 1));
  if(is.numeric(confidence.intervals)) stopifnot(confidence.intervals > 0 & confidence.intervals < 1)
  stopifnot(is.numeric(bootstraps))
  study.design <- match.arg(study.design)

  if(!missing(population.prevalence)) {
    stopifnot(is.numeric(population.prevalence))
    stopifnot(population.prevalence > 0 & population.prevalence < 1)
  }

  #check vars are in data
  if(any( names.check <- !is.element(all.vars(formula), names(data)))) stop(paste("variable(s)", paste( all.vars(formula)[names.check], collapse = ", ") , "not found in 'data'"))

  #throw out missing data
  data <- data[,all.vars(formula)]
  #complete case indicator
  cc.ind <- complete.cases(data)
  if(sum(cc.ind) < nrow(data)) warning(paste(sum(1-cc.ind), "observation(s) with missing data removed"))

  data <- data[cc.ind,]

  #study design
  if(missing(population.prevalence)) population.prevalence <- NULL

  if(study.design == "cohort"){
    if(!is.null(population.prevalence)){
      warning("population.prevalence was provided, but study.design = 'cohort'. The value input for population.prevalence will be ignored. If you are using case-control data, please set study.design = 'case-control'.")
    }
  }else{
    if(missing(population.prevalence)){
      stop("Need to set population.prevalence to calculate decision curves using case-control data.")
      if(family$family != "binomial") stop("Calculations for case-control data are done assuming logistic regression (family = binomial(link = 'logit'))")
    }else{
      stopifnot(0< population.prevalence & population.prevalence <1)
      message("Calculating net benefit curves for case-control data. All calculations are done conditional on the outcome prevalence provided.")

    }
  }

  #retreive outcome and check
  outcome <- data[[all.vars(formula[[2]])]];
  if(length(unique(outcome)) != 2) stop("outcome variable is not binary (it does not take two unique values).")
  stopifnot(is.numeric(outcome))
  if(min(outcome) != 0 | max(outcome) != 1) stop("outcome variable must be binary taking on values 0 for control and 1 for case.")


   #if fitted risks are provided, then there can only be one term on the rhs of formula
   #and the provided risks must be
  if(fitted.risk){
    #message("Fitted risks are provided, no model fitting will be done by DecisionCurve. Bootstrap confidence intervals are conditional on the model used to fit risks.")
    if(length(all.vars(formula[[3]])) > 1) stop("When fitted.risk = TRUE, there can only be one term  (denoting the fitted risks) on the right hand side of the formula provided.")

    provided.risks <-  data[[Reduce(paste, deparse(formula[[3]]))]] #get the name of the fitted risk variable from formula.
    if(min(provided.risks) < 0 | max(provided.risks) > 1) stop("When fitted.risks = TRUE, all risks provided must be between 0 and 1.")

  }
  #########
  ## End Checks
  #########

  #calculate curves
  #first we fit the model

  #extract the model name from formula
  predictors <- c(Reduce(paste, deparse(formula[[3]])), "All", "None")
  predictor.names <- c(Reduce(paste, deparse(formula)), "All", "None")

  #indicate whether we are fitting a model with a formula or not
  #the last two are FALSE since they correspond to 'all' and 'none'
  formula.ind <- c(ifelse(fitted.risk, FALSE, TRUE), FALSE, FALSE)

  data[["All"]] <- 1
  data[["None"]] <- 0

  n.preds <- length(predictors) #should always be three

  n.out <- length(predictors)*length(thresholds)
  dc.data <- data.frame("thresholds" = numeric(n.out),
                    "FPR" = numeric(n.out),"TPR" = numeric(n.out),
                    "NB" = numeric(n.out), "sNB" = numeric(n.out),
                    "rho" = numeric(n.out),"prob.high.risk" = numeric(n.out),
                    "DP" = numeric(n.out),
                    "model"= numeric(n.out))


  #if ci's
  #set up vars for bootstrap ci's and calculate bootstrap indices
  if(is.numeric(confidence.intervals))  {
    if(bootstraps < 1 ) stop("bootstraps must be greater than 0. If no confidence intervals are needed, set `confidence.intervals = 'none'`")
    #bootstrap sampling indices
    B.ind <- matrix(nrow = nrow(data), ncol = bootstraps)

    ## cohort design: don't stratify by outcome status.
    if(study.design == "cohort"){
        for(b in 1:bootstraps) B.ind[,b] <- sample.int(nrow(data), replace = TRUE)
    }else{
      #case-control design: stratify on outcome status for bootstrapping
        all.ind  <- 1:nrow(data)
        uu <- unique(outcome) #unique outcome levels

        for(b in 1:bootstraps){
          ind.1 <- sample(all.ind[outcome == uu[1] ], replace = TRUE)
          ind.2 <- sample(all.ind[outcome == uu[2] ], replace = TRUE)

          B.ind[,b] <- c(ind.1, ind.2)
        }
    }
    dc.data <- add.ci.columns(dc.data)
  }

  index = 1
  n.pred = 1
  #cycle through the model provided, all and none
  for(i in 1:n.preds){

    tmpNBdata <- calculate.nb(d = outcome,
                              y = data[[predictors[[i]]]],
                              rH = thresholds,
                              formula = formula,
                              family = family,
                              data = data,
                              formula.ind = formula.ind[i],
                              casecontrol.rho = population.prevalence)

    tmpNBdata$model <- predictor.names[[i]]

    #################
    ## Bootstrapping for CI's
    ###############

    if(is.numeric(confidence.intervals)){
      #calculate measures in each bootstrap
      boot.data <- apply(B.ind, 2, function(x){

      calculate.nb(d = outcome[x],
                   y = data[[predictors[[i]] ]][x],
                   rH = thresholds,
                   formula = formula,
                   family = family,
                   data = data[x,],
                   formula.ind = formula.ind[i],
                   casecontrol.rho = population.prevalence)})

      alpha = 1- confidence.intervals

      xx <- NULL #appease check

      #go through each measure and get the quantiles from the bootstrap distribution at each threshold
      for(rtn in names(boot.data[[1]][-1])){
        #collate the data from the measure estimates across bootstrap replicates
        #so I can call apply next.
        tmpdat <- sapply(boot.data, function(xx) xx[,rtn])

        tmpNBdata[[paste(rtn, "_lower", sep = "")]] <- apply(tmpdat, 1, quantile, probs = alpha/2, type = 1, na.rm = TRUE)
        tmpNBdata[[paste(rtn, "_upper", sep = "")]] <- apply(tmpdat, 1, quantile, probs = 1-alpha/2, type = 1, na.rm = TRUE)
      }
    }


    dc.data[index:(length(thresholds)*n.pred),] <-  tmpNBdata
    index = index + length(thresholds)
    n.pred = n.pred + 1
  }


  dc.data$cost.benefit.ratio <- as.character(fractions(threshold_to_costbenefit(dc.data$thresholds)))
  #find indices without a fraction and make them "xx/1"
  add.dash1 <-  which(!is.element(1:nrow(dc.data), grep("/", dc.data$cost.benefit.ratio)))
  dc.data$cost.benefit.ratio[add.dash1] <- paste(dc.data$cost.benefit.ratio[add.dash1], "/1", sep = "")
  dc.data$cost.benefit.ratio <- gsub("/", ":", dc.data$cost.benefit.ratio)



  #return list of elements
  out <- list("derived.data"  = dc.data,
              "confidence.intervals" = confidence.intervals,
              "call" = call)
  class(out) = "decision_curve"
  invisible(out)

}



