##' This function calculates scalar weights for use in other models
##'
##' This function calculates the scalar weights
##'
##'
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param numerator_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.  i.e.
##' \code{treat ~ 1}.
##' @param data is a dataframe containing \code{treat}, and
##' \code{X}.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Poisson"} for Poisson model,
##' \code{"Sqrt"} for square-root transformation
##' to a normal treatment,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model.
##' @param link_function is either "log", "inverse", or "identity" for the
##' "Gamma" \code{treat_mod}.
##' @param ...	additional arguments to be passed to the treatment regression fitting function.
##'
##' @return \code{scalar_wts}  returns an object of class "causaldrf_wts",
##' a list that contains the following components:
##'
##' \item{param}{summary of estimated weights.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{num_mod}{the result of the numerator model fit.}
##' \item{weights}{estimated weights for each unit.}
##' \item{call}{the matched call.}
##'
##' @seealso \code{\link{iptw_est}}, \code{\link{ismw_est}},
##'  \code{\link{reg_est}}, \code{\link{aipwee_est}}, \code{\link{wtrg_est}},
##'   etc. for other estimates.
##'
##' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
##' for use in the different estimates.
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' scalar_wts_list <- scalar_wts(treat = T,
##'                      treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'                      numerator_formula = T ~ 1,
##'                      data = example_data,
##'                      treat_mod = "Normal")
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'      scalar_wts_list$weights[sample_index],
##'      xlab = "T",
##'      ylab = "weights",
##'      main = "scalar_wts")
##'
##'
##' rm(example_data, scalar_wts_list, sample_index)
##'
##'
##' @export
##'
##'
##' @usage
##'
##' scalar_wts(treat,
##'            treat_formula,
##'            numerator_formula,
##'            data,
##'            treat_mod,
##'            link_function,
##'            ...)
##'

scalar_wts <- function(treat,
                       treat_formula,
                       numerator_formula,
                       data,
                       treat_mod,
                       link_function,
                       ...){

  # treat is the name of the treatment variable
  # treat_formula is the formula for the covariates model of the form: ~ X.1 + ....
  # numerator_formula is the formula for the numerator of the weights of the form: ~ X.1 + .... or ~ 1
  # data will contain all the data: X, treat
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is the estimated parameters.


  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("numerator_formula" %in% names(tempcall))) stop("No numerator_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) & !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Normal")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "NegBinom") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = NegBinom (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}




  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    treat = data[,as.character(tempcall$treat)]

  )


  formula_t <- tempcall$treat_formula
  formula_numerator = tempcall$numerator_formula

  samp <- data

  if (treat_mod == "Poisson"){

    result <- glm(formula_t,
                  family = "poisson",
                  data = samp,
                  ...)

    cond_mean <- result$fitted.values

    samp$gps_vals <- dpois(x = tempdat$treat,
                           lambda = cond_mean)

    # calculate numerator of the weight
    num_mod <-  glm(formula_numerator,
                    family = "poisson",
                    data = samp,
                    ...)

    cond_mean_num <- num_mod$fitted.values

    samp$num_weight <- dpois(x = tempdat$treat,
                             lambda = cond_mean_num)

    # weights
    est_import_wt <- samp$num_weight/samp$gps_vals

  } else if (treat_mod == "NegBinom"){

    result <- MASS::glm.nb(formula_t,
                           link = link_function,
                           data = samp,
                           ...)

    cond_mean <- result$fitted.values
    cond_var <- cond_mean + cond_mean^2/result$theta
    prob_nb_est <- (cond_var - cond_mean) / cond_var

    samp$gps_vals <- dnbinom(x = tempdat$treat,
                             size = result$theta,
                             mu = result$fitted.values,
                             log = FALSE)


    # calculate numerator of the weight
    num_mod <-  MASS::glm.nb(formula_numerator,
                             link = link_function,
                             data = samp,
                             ...)

    cond_mean_num <- num_mod$fitted.values
    cond_var_num <- cond_mean_num + cond_mean_num^2/num_mod$theta

    prob_nb_est_num <- (cond_var_num - cond_mean_num) / cond_var_num

    samp$num_weight <- dnbinom(x = tempdat$treat,
                               size = num_mod$theta,
                               prob = prob_nb_est_num,
                               mu = num_mod$fitted.values,
                               log = FALSE)

    # weights
    est_import_wt <- samp$num_weight/samp$gps_vals

  } else if (treat_mod == "Gamma"){

    result <- glm(formula_t,
                  family = Gamma(link = link_function),
                  data = samp,
                  ...)


    shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    theta_given_X <- result$fitted.values/shape_gamma

    samp$gps_vals <- dgamma(tempdat$treat,
                            shape = shape_gamma,
                            scale = theta_given_X)


    # calculate numerator of the weight
    num_mod <- glm(formula_numerator,
                   family = Gamma(link = "log"),
                   data = samp,
                   ...)

    shape_gamma_2 <- as.numeric(MASS::gamma.shape(num_mod)[1])
    theta_given_X_2 <- num_mod$fitted.values/shape_gamma

    samp$num_weight <- dgamma(tempdat$treat,
                              shape = shape_gamma_2,
                              scale = theta_given_X_2)

    # weights
    est_import_wt <- samp$num_weight/samp$gps_vals


  } else  if (treat_mod == "LogNormal"){


    samp[, as.character(tempcall$treat)] <- log(samp[, as.character(tempcall$treat)])


    result <- lm(formula_t,
                 data = samp,
                 ...)


    samp$gps_vals <- dnorm(samp[, as.character(tempcall$treat)],
                           mean = result$fitted,
                           sd = summary(result)$sigma   )


    num_mod <- lm(formula_numerator,
                  data = samp,
                  ...)

    samp$num_weight <- dnorm(samp[, as.character(tempcall$treat)],
                             mean = coef(num_mod)[1],
                             sd = summary(num_mod)$sigma  )
    est_import_wt <- samp$num_weight/samp$gps_vals


  } else  if (treat_mod == "Sqrt"){


    samp[, as.character(tempcall$treat)] <- sqrt(samp[, as.character(tempcall$treat)])


    result <- lm(formula_t,
                 data = samp,
                 ...)


    samp$gps_vals <- dnorm(samp[, as.character(tempcall$treat)],
                           mean = result$fitted,
                           sd = summary(result)$sigma   )


    num_mod <- lm(formula_numerator,
                  data = samp,
                  ...)

    samp$num_weight <- dnorm(samp[, as.character(tempcall$treat)],
                             mean = coef(num_mod)[1],
                             sd = summary(num_mod)$sigma  )
    est_import_wt <- samp$num_weight/samp$gps_vals


  } else if (treat_mod == "Normal"){

    result <- lm(formula_t,
                 data = samp,
                 ...)

    samp$gps_vals <- dnorm(tempdat$treat,
                           mean = result$fitted,
                           sd = summary(result)$sigma   )
    num_mod <- lm(formula_numerator,
                  data = samp,
                  ...)

    samp$num_weight <- dnorm(tempdat$treat,
                             mean = coef(num_mod)[1],
                             sd = summary(num_mod)$sigma  )
    est_import_wt <- samp$num_weight/samp$gps_vals

  } else {
    stop("Error: treat_mod not recognized!!!")
  }






  z_object <- list(params = summary(est_import_wt),
                   t_mod = result,
                   num_mod = num_mod,
                   weights = est_import_wt,
                   call = tempcall)

  class(z_object) <- "causaldrf_wts"
  z_object

}




##' @export
summary.causaldrf_wts <- function(object,...){
  cat("\nSummary values:\n")
  print(object$param)
  cat("\nTreatment Summary:\n")
  print(summary(object$t_mod))
  cat("\nNumerator Summary:\n")
  print(summary(object$num_mod))

}


##' @export
print.summary.causaldrf_wts <- function(x, ...){
  cat("\nTreatment Summary:\n")
  print(summary(x$t_mod))
  cat("\nNumerator Summary:\n")
  print(summary(x$num_mod))
}


##' @export
print.causaldrf_wts <- function(x,...){
  #   cat("Call:\n")
  #   print(x$call)
  cat("\nEstimated values:\n")
  print(x$param)
}



