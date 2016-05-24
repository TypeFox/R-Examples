##' The inverse weighting estimator (nonparametric method)
##'
##'
##' This is a nonparametric method that estimates the ADRF by using a local linear
##' regression of \code{Y} on \code{treat} with weighted kernel function.  For
##' details, see Flores et. al. (2012).
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
##' @param bandw is the bandwidth.  Default is 1.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model.
##' @param link_function is either "log", "inverse", or "identity" for the
##' "Gamma" \code{treat_mod}.
##' @param ... additional arguments to be passed to the treatment regression function.
##'
##' @details
##'
##' The ADRF is estimated by
##' \deqn{(D_{0}(t) S_{2}(t) - D_{1}(t) S_{1}(t)) / (S_{0}(t) S_{2}(t) - S_{1}^{2}(t)) }
##' where \deqn{D_{j}(t) = \sum_{i = 1}^{N} \tilde{K}_{h, X} (T_i - t)  (T_i - t)^j Y_i} and
##' \eqn{S_{j}(t) = \sum_{i = 1}^{N} \tilde{K}_{h, X} (T_i - t)  (T_i - t)^j}
##' \eqn{\tilde{K}_{h, X}(t) = K_{h}(t) / \hat{R}_i(t)}
##' which is a local linear regression.
##' More details are given in Flores (2012).
##'
##' @return \code{iw_est} returns an object of class "causaldrf",
##' a list that contains the following components:
##' \item{param}{parameter estimates for a iw fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{call}{the matched call.}
##'
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' Flores, Carlos A., et al. "Estimating the effects of length of exposure to
##' instruction in a training program: the case of job corps."
##' \emph{Review of Economics and Statistics} \bold{94.1} (2012): 153-171.
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' iw_list <- iw_est(Y = Y,
##'                 treat = T,
##'                 treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'                 data = example_data,
##'                 grid_val = seq(8, 16, by = 1),
##'                 bandw = bw.SJ(example_data$T),
##'                 treat_mod = "Normal")
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "iw estimate")
##'
##' lines(seq(8, 16, by = 1),
##'         iw_list$param,
##'         lty = 2,
##'         lwd = 2,
##'         col = "blue")
##'
##' legend('bottomright',
##'         "iw estimate",
##'         lty=2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y',
##'         cex=1)
##'
##' rm(example_data, iw_list, sample_index)
##'
##' ## Example from Imai & van Dyk (2004).
##'
##' data("nmes_data")
##' head(nmes_data)
##' # look at only people with medical expenditures greater than 0
##' nmes_nonzero <- nmes_data[which(nmes_data$TOTALEXP > 0), ]
##'
##'
##' iw_list <- iw_est(Y = TOTALEXP,
##'                   treat = packyears,
##'                   treat_formula = packyears ~ LASTAGE + I(LASTAGE^2) +
##'                     AGESMOKE + I(AGESMOKE^2) + MALE + RACE3 + beltuse +
##'                     educate + marital + SREGION + POVSTALB,
##'                   data = nmes_nonzero,
##'                   grid_val = seq(5, 100, by = 5),
##'                   bandw = bw.SJ(nmes_nonzero$packyears),
##'                   treat_mod = "LogNormal")
##'
##' set.seed(307)
##' sample_index <- sample(1:nrow(nmes_nonzero), 1000)
##'
##' plot(nmes_nonzero$packyears[sample_index],
##'      nmes_nonzero$TOTALEXP[sample_index],
##'      xlab = "packyears",
##'      ylab = "TOTALEXP",
##'      main = "iw estimate",
##'      ylim = c(0, 10000),
##'      xlim = c(0, 100))
##'
##' lines(seq(5, 100, by = 5),
##'       iw_list$param,
##'       lty = 2,
##'       lwd = 2,
##'       col = "blue")
##'
##' legend('topright',
##'        "iw estimate",
##'        lty=2,
##'        lwd = 2,
##'        col = "blue",
##'        bty='Y',
##'        cex = 1)
##' abline(0, 0)
##'
##'
##' @usage
##'
##' iw_est(Y,
##'        treat,
##'        treat_formula,
##'        data,
##'        grid_val,
##'        bandw,
##'        treat_mod,
##'        link_function,
##'        ...)
##'
##' @seealso \code{\link{nw_est}}, \code{\link{iw_est}}, \code{\link{hi_est}}, \code{\link{gam_est}},
##' \code{\link{add_spl_est}},
##' \code{\link{bart_est}}, etc. for other estimates.
##'
##' @export
##'



iw_est <- function (Y,
                    treat,
                    treat_formula,
                    data,
                    grid_val,
                    bandw,
                    treat_mod,
                    link_function,
                    ...)
{
  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # treat_formula is the formula for the treatment model
  # data will contain all the data: X, treat, and Y
  # grid_val is the set of grid points on T
  # bandw is the bandwidth
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is a list of 2 objects:
  #    (1) estimated values of ADRF at gridvalues
  #    (2) result of treatment model


  #save input
  tempcall <- match.call()


  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("grid_val" %in% names(tempcall)))  stop("No grid_val specified")
  if (!("bandw" %in% names(tempcall)))  stop("No bandw specified")
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) & !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Normal")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}



  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]

  )


  # make a formula for the treatment model
  #   formula_t <- eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$treat_formula, width.cutoff = 500), sep = "")))
  formula_t <- tempcall$treat_formula

  if (treat_mod == "NegBinom") {
    samp_dat <- data
    result <- MASS::glm.nb(formula = formula_t, link = log, data = samp_dat, ...)
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
    result <- glm(formula = formula_t, family = "poisson", data = samp_dat, ...)
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
                  data = samp_dat, ...)
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
    result <- lm(formula = formula_t, data = samp_dat, ...)
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
    result <- lm(formula = formula_t, data = samp_dat, ...)
    est_sqrt_treat <- result$fitted
    sigma_est <- summary(result)$sigma
    gps_fun_sqrt <- function(tt) {
      dnorm(sqrt(tt), mean = est_sqrt_treat, sd = sigma_est)
    }
    gps_fun <- gps_fun_sqrt
  }
  else if (treat_mod == "Normal") {
    samp_dat <- data
    result <- lm(formula = formula_t, data = samp_dat, ...)
    gps_fun_Normal <- function(tt) {
      dnorm(tt, mean = result$fitted, sd = summary(result)$sigma)
    }
    gps_fun <- gps_fun_Normal
  }
  else {
    stop("No valid treat_mod specified.  Please try again)")
  }


#--------------------------------
if (treat_mod == "LogNormal"){

# K_h <- function(t, h){dnorm( (log(treat) - log(t) )/h, 0, 1) / h }
#
# K_h_tild <- function(t, h){K_h(t, h) / gps_fun(t)}

  K_h <- function(t, h) {
    dnorm((log(tempdat$treat) - log(t))/h, 0, 1)/h
  }
  K_h_tild <- function(t, h) {
    K_h(t, h)/gps_fun(t)
  }

S_func <- function(grid_val, h, j){

  S_function <- sum( K_h_tild(grid_val, h) * (tempdat$treat - grid_val) ^ j )
  return(S_function)
}


D_func <- function(grid_val, h, j){

  D_function <- sum( K_h_tild(grid_val, h) * (tempdat$treat - grid_val) ^ j * tempdat$Y)
  return(D_function)
}


# iw_est uses S_func and D_func to estimate the ADRF.
iw_estimator <- function(grid_val, h){

  iw_function <- numeric(length(grid_val))

  for (i in 1:length(grid_val)){
    iw_function[i] <- (D_func(grid_val[i], h, j = 0) *
                         S_func(grid_val[i], h, j = 2) -
                         D_func(grid_val[i], h, j = 1) *
                         S_func(grid_val[i], h, j = 1)
    )    /

      (S_func(grid_val[i], h, j = 0)
       * S_func(grid_val[i], h, j = 2) -
         S_func(grid_val[i], h, j = 1)^2

      )

  }

  return(iw_function)

}

} else if (treat_mod == "Sqrt"){

  # K_h <- function(t, h){dnorm( (sqrt(treat) - sqrt(t) )/h, 0, 1) / h }
  #
  # K_h_tild <- function(t, h){K_h(t, h) / gps_fun(t)}

  K_h <- function(t, h) {
    dnorm((sqrt(tempdat$treat) - sqrt(t))/h, 0, 1)/h
  }
  K_h_tild <- function(t, h) {
    K_h(t, h)/gps_fun(t)
  }

  S_func <- function(grid_val, h, j){

    S_function <- sum( K_h_tild(grid_val, h) * (tempdat$treat - grid_val) ^ j )
    return(S_function)
  }


  D_func <- function(grid_val, h, j){

    D_function <- sum( K_h_tild(grid_val, h) * (tempdat$treat - grid_val) ^ j * tempdat$Y)
    return(D_function)
  }


  # iw_est uses S_func and D_func to estimate the ADRF.
  iw_estimator <- function(grid_val, h){

    iw_function <- numeric(length(grid_val))

    for (i in 1:length(grid_val)){
      iw_function[i] <- (D_func(grid_val[i], h, j = 0) *
                           S_func(grid_val[i], h, j = 2) -
                           D_func(grid_val[i], h, j = 1) *
                           S_func(grid_val[i], h, j = 1)
      )    /

        (S_func(grid_val[i], h, j = 0)
         * S_func(grid_val[i], h, j = 2) -
           S_func(grid_val[i], h, j = 1)^2

        )

    }

    return(iw_function)

  }

} else {


  K_h <- function(t, h) {
    dnorm((tempdat$treat - t)/h, 0, 1)/h
  }
  K_h_tild <- function(t, h) {
    K_h(t, h)/gps_fun(t)
  }

  S_func <- function(grid_val, h, j){

    S_function <- sum( K_h_tild(grid_val, h) * (tempdat$treat - grid_val) ^ j )
    return(S_function)
  }


  D_func <- function(grid_val, h, j){

    D_function <- sum( K_h_tild(grid_val, h) * (tempdat$treat - grid_val) ^ j * tempdat$Y)
    return(D_function)
  }


  # iw_est uses S_func and D_func to estimate the ADRF.
  iw_estimator <- function(grid_val, h){

    iw_function <- numeric(length(grid_val))

    for (i in 1:length(grid_val)){
      iw_function[i] <- (D_func(grid_val[i], h, j = 0) *
                           S_func(grid_val[i], h, j = 2) -
                           D_func(grid_val[i], h, j = 1) *
                           S_func(grid_val[i], h, j = 1)
      )    /

        (S_func(grid_val[i], h, j = 0)
         * S_func(grid_val[i], h, j = 2) -
           S_func(grid_val[i], h, j = 1)^2

        )

    }

    return(iw_function)

  }

}


  z_object <- list(param = iw_estimator(grid_val, bandw),
                   t_mod = result,
                   call = tempcall)

  class(z_object) <- "causaldrf"
  z_object

}
