#' A function to estimate conditional expected values and higher order moments
#'
#' This function fits a glm regression specified by the user to
#' estimate conditional moments.
#'
#'
#'
#' @param treat is the name of the treatment variable contained in \code{data}.
#' @param treat_formula an object of class "formula" (or one that can be coerced to that class)
#' that regresses \code{treat} on a linear combination of \code{X}:
#' a symbolic description of the model to be fitted.
#' @param data is a dataframe containing \code{Y}, \code{treat}, and \code{X}.
#' @param treat_mod a description of the error distribution
#' to be used in the model for treatment. Options include:
#' \code{"Normal"} for normal model,
#' \code{"LogNormal"} for lognormal model,
#' \code{"Sqrt"} for square-root transformation
#' to a normal treatment,
#' \code{"Poisson"} for Poisson model,
#' \code{"NegBinom"} for negative binomial model,
#' \code{"Gamma"} for gamma model.
#' @param link_function is either "log", "inverse", or "identity" for the "Gamma" \code{treat_mod}.
#'
#' @param ...	additional arguments to be passed to the low level treatment regression fitting functions.
#'
#' @return
#' \code{t_mod} returns a list containing the following elements:
#' \item{T_data}{a dataframe containing estimated treatment, estimated treatment squared,
#' estimated treatment cube, estimated treatment quartic, and estimated gps.}
#' \item{T_result}{the result of the treatment model fit.}
#'
#' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
#' continuous treatment and outcome: alternative estimators for parametric
#' dose-response models. \emph{Manuscript in preparation}.
#'
#' @seealso   \code{\link{ismw_est}}, \code{\link{reg_est}},
#' \code{\link{wtrg_est}}, \code{\link{aipwee_est}}, etc. for other estimates.
#'
#'
#' \code{\link{overlap_fun}} to prepare the \code{data} for use in the different estimates.
#'
#' @examples
#'
#' ## Example from Schafer (2015).
#'
#' example_data <- sim_data
#'
#' t_mod_list <- t_mod(treat = T,
#'                     treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
#'                     data = example_data,
#'                     treat_mod = "Normal")
#'
#' cond_exp_data <- t_mod_list$T_data
#'
#' full_data <- cbind(example_data, cond_exp_data)
#'
#' rm(example_data, t_mod_list, cond_exp_data, full_data)
#'
#'
#' @usage
#'
#' t_mod(treat,
#'       treat_formula,
#'       data,
#'       treat_mod,
#'       link_function,
#'       ...)
#'
#' @export
#'


t_mod <- function(treat,
                  treat_formula,
                  data,
                  treat_mod,
                  link_function,
                  ...){

  # treat is the name of the treatment variable
  # treat_formula is the formula for the treatment model
  # data will contain all the data: X, treat
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is a list of 2 objects:
  #    (1) a data.frame containing the estimated conditional higher order moments.
  #    (2) the result from the treament model fit.



  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) & !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Normal")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}



  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    treat = data[,as.character(tempcall$treat)]

  )

  samp_dat <- data

  # make a formula for the treatment model
  #   formula_t <- eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$treat_formula, width.cutoff = 500), sep = "")))
  formula_t <- tempcall$treat_formula


  if (treat_mod == "NegBinom"){


    result <- MASS::glm.nb(formula_t,
                           data = samp_dat,
                           ...)

    cond_mean <- result$fitted.values
    cond_var <- cond_mean + cond_mean^2/result$theta
    prob_nb_est <- (cond_var - cond_mean) / cond_var

    pp <- (cond_var^2 - cond_mean) / cond_var^2
    rr <- cond_mean^2 / (cond_var^2 - cond_mean)

    # Need to get estimates of higher moments of NegBinom distributed RV using the MGF
    # fill in the formulas....
    est_treat <- pp * rr /(1 - pp)
    est_treat_sq <- pp * rr * (pp * rr + 1) / (pp - 1)^2
    est_treat_cube <- - pp * rr * (pp^2 * rr^2 + 3 * pp * rr + pp + 1) / (pp - 1)^3
    est_treat_quartic <- (pp^4 * rr^4 + 6 * pp^3 * rr^3 + pp^2 * (4 * pp + 7) * rr^2 + pp * (pp^2 + 4 * pp + 1) * rr) / (pp - 1)^4

    gps <- dnbinom(x = tempdat$treat,
                   size = result$theta,
                   mu = result$fitted.values,
                   log = FALSE)


    t_mod_data <- data.frame(cbind(est_treat,
                                   est_treat_sq,
                                   est_treat_cube,
                                   est_treat_quartic,
                                   gps))

  } else if (treat_mod == "Poisson"){


    result <- glm(formula_t,
                  family = "poisson",
                  data = samp_dat,
                  ...)

    cond_mean <- result$fitted.values

    # Need to get estimates of higher moments of gamma distributed RV using the MGF
    est_treat <- cond_mean
    est_treat_sq <- cond_mean^2 + cond_mean
    est_treat_cube <- cond_mean^3 + 3*cond_mean^2 + cond_mean
    est_treat_quartic <- cond_mean^4 + 6*cond_mean^2 + 7*cond_mean^2 + cond_mean

    gps <- dpois(x = tempdat$treat,
                 lambda = cond_mean)


    t_mod_data <- data.frame(cbind(est_treat,
                                   est_treat_sq,
                                   est_treat_cube,
                                   est_treat_quartic,
                                   gps))


  } else if (treat_mod == "Gamma"){


    result <- glm(formula_t,
                  family = Gamma(link = link_function),
                  data = samp_dat,
                  ... )

    est_treat <- result$fitted

    shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    theta_given_X <- result$fitted.values/shape_gamma

    GPS_gamma <- dgamma(tempdat$treat,
                        shape = shape_gamma,
                        scale = theta_given_X)


    # Need to get estimates of higher moments of gamma distributed RV using the MGF
    est_treat <- shape_gamma * theta_given_X
    est_treat_sq <- theta_given_X^2 * shape_gamma * (shape_gamma + 1)
    est_treat_cube <- theta_given_X^3 * shape_gamma * (shape_gamma + 1) * (shape_gamma + 2)
    est_treat_quartic <- theta_given_X^4 * shape_gamma * (shape_gamma + 1) * (shape_gamma + 2) * (shape_gamma + 3)



    gps <- dgamma(tempdat$treat,
                  shape = shape_gamma,
                  scale = theta_given_X)

    t_mod_data <- data.frame(cbind(est_treat,
                                   est_treat_sq,
                                   est_treat_cube,
                                   est_treat_quartic,
                                   gps))

  } else if (treat_mod == "LogNormal"){


  samp_dat[, as.character(tempcall$treat)] <- log(samp_dat[, as.character(tempcall$treat)])



  result <- lm(formula_t,
               data = samp_dat,
               ... )

  est_log_treat <- result$fitted
  sigma_est <- summary(result)$sigma


  # moment function calculator for lognormal RV
  # mu is the mean from normal and sigma is from the normal density
  moment_fun <- function(mu,
                         sigma,
                         k){
    return(exp(k * mu + k^2 * sigma^2 / 2) )
  }

  # Need to get estimates of higher moments of normally distributed RV using the MGF
  est_treat <- moment_fun(est_log_treat,
                          sigma_est,
                          1)
  est_treat_sq <- moment_fun(est_log_treat,
                             sigma_est,
                             2)
  est_treat_cube <- moment_fun(est_log_treat,
                               sigma_est,
                               3)
  est_treat_quartic <- moment_fun(est_log_treat,
                                  sigma_est,
                                  4)


  gps <- dnorm((samp_dat[, as.character(tempcall$treat)] - result$fitted)/summary(result)$sigma,
               mean = 0,
               sd = 1 )

  t_mod_data <- data.frame(cbind(est_treat,
                                 est_treat_sq,
                                 est_treat_cube,
                                 est_treat_quartic,
                                 gps))


  }  else if (treat_mod == "Sqrt"){


    samp_dat[, as.character(tempcall$treat)] <- sqrt(samp_dat[, as.character(tempcall$treat)])



    result <- lm(formula_t,
                 data = samp_dat,
                 ... )

    est_sqrt_treat <- result$fitted
    sigma_est <- summary(result)$sigma




    # Need to get estimates of higher moments of normally distributed RV using the MGF
    est_treat <- result$fitted^2 +
      var(result$fitted)

    est_treat_sq <- result$fitted^4 +
      6 * result$fitted^2 * var(result$fitted) +
      3 * var(result$fitted)^2

    est_treat_cube <- result$fitted^6 +
      15 * result$fitted^4 * var(result$fitted) +
      45 * result$fitted^2 * var(result$fitted)^2 +
      15 * var(result$fitted)^3

    est_treat_quartic <- result$fitted^8 +
      28 * result$fitted^6 * var(result$fitted) +
      210 * result$fitted^4 * var(result$fitted)^2 +
      420 * result$fitted^2 * var(result$fitted)^3 +
      105 * var(result$fitted)^4


    gps <- dnorm((samp_dat[, as.character(tempcall$treat)] - result$fitted)/summary(result)$sigma,
                 mean = 0,
                 sd = 1 )

    t_mod_data <- data.frame(cbind(est_treat,
                                   est_treat_sq,
                                   est_treat_cube,
                                   est_treat_quartic,
                                   gps))


  } else if (treat_mod == "Normal"){



  result <- lm(formula_t,
               data = samp_dat,
               ... )


  est_treat <- result$fitted
  sigma_est <- summary(result)$sigma


  # Need to get estimates of higher moments of normally distributed RV using the MGF
  est_treat <- result$fitted
  est_treat_sq <- result$fitted^2 + var(result$fitted)
  est_treat_cube <- result$fitted^3 + 3 * result$fitted * var(result$fitted)
  est_treat_quartic <- result$fitted^4 + 6 * result$fitted^2 * var(result$fitted) + 3 * var(result$fitted)^2


  gps <- dnorm(tempdat$treat,
               mean = est_treat,
               sd = summary(result)$sigma )

  t_mod_data <- data.frame(cbind(est_treat,
                                 est_treat_sq,
                                 est_treat_cube,
                                 est_treat_quartic,
                                 gps))



  } else {
    stop("No valid treatment model chosen.  Please try again.")
  }


  return(list(T_data = t_mod_data, T_result = result))

}


