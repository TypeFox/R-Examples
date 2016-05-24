##' The additive spline estimator
##'
##' This function estimates the ADRF with an additive spline estimator described
##' in Bia et al. (2014).
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param grid_val contains the treatment values to be evaluated.
##' @param knot_num is the number of knots used in outcome model
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model.
##' @param link_function is either "log", "inverse", or "identity" for the
##' "Gamma" \code{treat_mod}.
##' @param ...	additional arguments to be passed to the outcome regression fitting function.
##'
##' @details This function estimates the ADRF using additive splines in the
##' outcome model described in Bia et al. (2014).
##'
##' @return \code{add_spl_est}  returns an object of class "causaldrf",
##' a list that contains the following components:
##'
##' \item{param}{parameter estimates for a add_spl fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{out_mod}{the result of the outcome model fit.}
##' \item{call}{the matched call.}
##'
##' @seealso \code{\link{nw_est}}, \code{\link{iw_est}}, \code{\link{hi_est}}, \code{\link{gam_est}},
##' \code{\link{bart_est}}, etc. for other estimates.
##'
##' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
##' for use in the different estimates.
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' Bia, Michela, et al. (2014). A Stata package for the application of
##' semiparametric estimators of dose response functions.  \emph{Stata Journal}
##' \bold{14.3}, 580-604.
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##' example_data <- sim_data
##' add_spl_list <- add_spl_est(Y = Y,
##'             treat = T,
##'             treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'             data = example_data,
##'             grid_val = seq(8, 16, by = 1),
##'             knot_num = 3,
##'             treat_mod = "Normal")
##'
##'
##' sample_index <- sample(1:1000, 100)
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "additive spline estimate")
##'
##' lines(seq(8, 16, by = 1),
##'       add_spl_list$param,
##'       lty = 2,
##'       lwd = 2,
##'       col = "blue")
##' legend('bottomright',
##'         "additive spline estimate",
##'         lty=2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y', cex=1)
##'
##' rm(example_data, add_spl_list, sample_index)
##'
##' ## See Vignette for more examples.
##'
##'
##' @usage
##'
##' add_spl_est(Y,
##'             treat,
##'             treat_formula,
##'             data,
##'             grid_val,
##'             knot_num,
##'             treat_mod,
##'             link_function,
##'             ...)
##'
##'
##'
##' @export


add_spl_est <- function(Y,
                        treat,
                        treat_formula,
                        data,
                        grid_val,
                        knot_num,
                        treat_mod,
                        link_function,
                        ...){

  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # treat_formula is the formula for the treatment model
  # data will contain all the data: X, treat, and Y
  # grid_val is the set of grid points on T
  # knot_num is the number of knots used in outcome model
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is a list of 3 objects:
  #    (1) estimated values of ADRF at gridvalues
  #    (2) result of treatment model
  #    (3) result from the outcome model: outcome regressed on


  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("grid_val" %in% names(tempcall)))  stop("No grid_val specified")
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) & !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Normal")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}


  sample_data <- data
  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(Y = sample_data[, as.character(tempcall$Y)])
  tempdat$treat <- sample_data[,as.character(tempcall$treat)]

  #-------------------------------------------------------------------------------

  # make a formula for the treatment model
  #   formula_t <- eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$treat_formula, width.cutoff = 500), sep = "")))
  formula_t <- tempcall$treat_formula

  if (treat_mod == "NegBinom") {
    samp_dat <- data
    result <- MASS::glm.nb(formula = formula_t, link = log, data = samp_dat)
    cond_mean <- result$fitted.values
    cond_var <- cond_mean + cond_mean^2/result$theta
    prob_nb_est <- (cond_var - cond_mean)/cond_var
    gps_fun_NB <- function(tt) {
      dnbinom(x = tt, size = result$theta, mu = result$fitted.values,
              log = FALSE)
    }
    gps_fun <- gps_fun_NB
  }
  else if (treat_mod == "Poisson") {
    samp_dat <- data
    result <- glm(formula = formula_t, family = "poisson", data = samp_dat)
    cond_mean <- result$fitted.values
    samp_dat$gps_vals <- dpois(x = tempdat$treat, lambda = cond_mean)
    gps_fun_Pois <- function(t) {
      dpois(t, lambda = cond_mean)
    }
    gps_fun <- gps_fun_Pois
  }
  else if (treat_mod == "Gamma") {
    samp_dat <- data
    result <- glm(formula = formula_t, family = Gamma(link = link_function),
                  data = samp_dat)
    est_treat <- result$fitted
    shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    theta_given_X <- result$fitted.values/shape_gamma
    theta_treat_X <- tempdat$treat/shape_gamma
    gps_fun_Gamma <- function(t) {
      dgamma(t, shape = shape_gamma, scale = theta_given_X)
    }
    gps_fun <- gps_fun_Gamma
  }
  else if (treat_mod == "LogNormal") {

    samp_dat <- data
    samp_dat[, as.character(tempcall$treat)] <- log(samp_dat[, as.character(tempcall$treat)])
    result <- lm(formula = formula_t, data = samp_dat)
    est_log_treat <- result$fitted
    sigma_est <- summary(result)$sigma
    gps_fun_Log <- function(tt) {
      dnorm(log(tt), mean = est_log_treat, sd = sigma_est)
    }
    gps_fun <- gps_fun_Log
  }
  else if (treat_mod == "Sqrt") {

    samp_dat <- data
    samp_dat[, as.character(tempcall$treat)] <- sqrt(samp_dat[, as.character(tempcall$treat)])
    result <- lm(formula = formula_t, data = samp_dat)
    est_sqrt_treat <- result$fitted
    sigma_est <- summary(result)$sigma
    gps_fun_sqrt <- function(tt) {
      dnorm(sqrt(tt), mean = est_sqrt_treat, sd = sigma_est)
    }
    gps_fun <- gps_fun_sqrt
  }
  else if (treat_mod == "Normal") {
    samp_dat <- data
    result <- lm(formula = formula_t, data = samp_dat)
    gps_fun_Normal <- function(tt) {
      dnorm(tt, mean = result$fitted, sd = summary(result)$sigma)
    }
    gps_fun <- gps_fun_Normal
  }
  else {
    stop("No valid treat_mod specified.  Please try again.")
  }



  #--------------------------------------------------------



  #------------------------------------------------------------------

  # create additive bases ----------------------------------------------------

  # # need to use indicator functions creating additive spline bases
  indic_terms <-function(t, knot_val) ifelse(t > knot_val, t - knot_val, 0)


  knot_values_treat <- attr(splines::ns(tempdat$treat, df = (knot_num + 1)), "knots")   # pick knot_num interior knots
  knot_values_gps <- attr(splines::ns(gps_fun(tempdat$treat), df = (knot_num + 1)), "knots")   # pick knot_num interior knots


  # generate treat bases
  treat_mat <- matrix(numeric(length(knot_values_treat) * nrow(sample_data) ), nrow = nrow(sample_data) )
  for (i in 1:length(knot_values_treat)){
    treat_mat[,i] <- indic_terms(tempdat$treat, knot_values_treat[i])

  }

  # generate gps bases
  gps_mat <- matrix(numeric(length(knot_values_gps) * nrow(sample_data) ), nrow = nrow(sample_data) )
  for (i in 1:length(knot_values_gps)){
    gps_mat[,i] <- indic_terms(gps_fun(tempdat$treat), knot_values_gps[i])

  }


  # fit the model ------------------------------------------------------------
  out_data <- data.frame( cbind(tempdat$Y, tempdat$treat, gps_fun(tempdat$treat), treat_mat, gps_mat) )


  # function to give generic names to columns of data frame
  name_gen <- function(mydf, var_name){paste(rep(var_name,ncol(mydf)),
                                             c(1:ncol(mydf)),
                                             sep="")
  }

  colnames(out_data) <- c("Y", "treat", "gps", name_gen(treat_mat, "treat_"), name_gen(gps_mat, "gps_"))


  add_spl_fit <- lm(Y ~ . , data = out_data, ...  )

  coef_add_spl <- add_spl_fit$coefficients


  # predict the outcome values at the grid points ----------------------------
  gps_at_grid <- matrix(numeric(length(grid_val) * nrow(out_data)), nrow = nrow(out_data) )
  for (i in 1:length(grid_val)) {
    gps_at_grid[, i] <-  gps_fun(grid_val[i])
  }


  add_spl_pred <- numeric( length(grid_val) )

  for ( i in 1:length(grid_val) ) {
    # for each grid_val, generate new bases, new gps_at_grid

    # generate treat bases based on a fixed grid_val
    treat_mat_pred <- matrix(numeric(length(knot_values_treat) * nrow(out_data) ), nrow = nrow(out_data) )
    for (j in 1:length(knot_values_treat)){
      treat_mat_pred[,j] <- indic_terms(grid_val[i], knot_values_treat[j])

    }

    # generate gps bases based on a fixed grid_val
    gps_mat_pred <- matrix(numeric(length(knot_values_gps) * nrow(out_data) ), nrow = nrow(out_data) )
    for (j in 1:length(knot_values_gps)){
      gps_mat_pred[,j] <- indic_terms(gps_fun(grid_val[i]), knot_values_gps[j])

    }



    add_spl_data <- data.frame(cbind(1,
                                     grid_val[i],
                                     gps_at_grid[,i],
                                     treat_mat_pred,
                                     gps_mat_pred)
    )


    colnames(add_spl_data) <- c("Y", "treat", "gps", name_gen(treat_mat, "treat_"), name_gen(gps_mat, "gps_"))

    pred_est <-predict(add_spl_fit,
                       type = "response",
                       newdata = add_spl_data
    )

    add_spl_pred[i] <- mean(pred_est)

  }

  add_spl_estimate <- add_spl_pred

  ###Saving results


  z_object <- list(param = add_spl_estimate,
                   t_mod = result,
                   out_mod = add_spl_fit,
                   call = tempcall)

  class(z_object) <- "causaldrf"
  z_object


}




##' @export
summary.causaldrf <- function(object,...){
  cat("\nEstimated values:\n")
  print(object$param)
  cat("\nTreatment Summary:\n")
  print(summary(object$t_mod))
  cat("\nOutcome Summary:\n")
  print(summary(object$out_mod))

}


##' @export
print.summary.causaldrf <- function(x, ...){
  cat("\nTreatment Summary:\n")
  print(summary(x$t_mod))
  cat("\nOutcome Summary:\n")
  print(summary(x$out_mod))
}


##' @export
print.causaldrf <- function(x,...){
#   cat("Call:\n")
#   print(x$call)
  cat("\nEstimated values:\n")
  print(x$param)
}


