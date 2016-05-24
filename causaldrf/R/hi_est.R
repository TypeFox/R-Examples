##' The Hirano and Imbens estimator
##'
##' This function estimates the GPS function and estimates the
##' ADRF. The GPS score is based on different treatment models.  The treatment is
##' linearly related to Xs.
##'
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param outcome_formula is the formula used for fitting the outcome surface.
##' gps is one of the independent variables to use in the outcome_formula. ie.
##'
##' \code{Y ~ treat+ I(treat^2) + gps + I(gps^2) + treat * gps}
##'
##'  or a variation
##' of this.  Use \code{gps} as the name of the variable representing the gps
##' in \code{outcome_formula}.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param grid_val contains the treatment values to be evaluated.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model, \code{"Binomial"} for binomial model.
##' @param link_function For \code{treat_mod = "Gamma"} (fitted using glm) alternatives are "log" or "inverse".
##' For \code{treat_mod = "Binomial"} (fitted using glm) alternatives are "logit", "probit", "cauchit", "log" and "cloglog".
##' @param ... additional arguments to be passed to the outcome lm() function.
##'
##' @details
##' Hirano (2004) (HI) introduced this imputation-type method that
##' includes a GPS component.  The idea is to fit a parametric observable
##' (outcome) model, which includes the estimated GPS as a covariate, to impute
##' missing potential outcomes.
##'
##' The method requires several steps.  First, a model is used to relate
##' treatment to the recorded covariates.  For example,
##' \eqn{T_i|\textbf{X}_i \sim \mathcal{N}(\textbf{X}_i^T \boldsymbol{\beta}, \sigma^2)}
##' and then estimate the \eqn{\boldsymbol{\beta}} parameters.
##' Next, the GPS for each unit is estimated
##'
##' \deqn{
##' \hat{R}_i(t) = \frac{1}{\sqrt{2 \pi \hat{\sigma}^2}  } e^{-\frac{(t - \textbf{X}_i^T \boldsymbol{\hat{\beta}})^2}{2 \hat{\sigma}^2}}
##' }
##' These GPS estimates are used in the outcome or observable model.
##' The outcome is modeled as a function of \eqn{T_i} and \eqn{\hat{R}_i}
##' parametrically.  For example,
##' \deqn{
##' E[Y_i | T_i, R_i] = \alpha_0 + \alpha_1 T_i + \alpha_2 T_i^2 + \alpha_3 \hat{R}_i + \alpha_4 \hat{R}_i^2 + \alpha_5 \hat{R}_i\cdot T_i
##' }
##' After collecting the estimated parameters in the outcome and treatment
##' models, plug-in the treatment values into the model to estimate the
##' missing potential outcomes of each individual at that treatment level.
##' For example, if we plug in \eqn{T_i = t} into the estimated models, then
##' each unit will have a potential outcome estimated at treatment level
##' \eqn{T_i= t}.
##'
##' \deqn{
##' \hat{Y}_i(t) = \hat{\alpha}_0 + \hat{\alpha}_1 t + \hat{\alpha}_2 t^2 + \hat{\alpha}_3 \hat{R}_i(t) + \hat{\alpha}_4 \hat{R}_i^2(t) + \hat{\alpha}_5 \hat{R}_i(t) \cdot t
##' }
##' The next step is to aggregate these estimated potential outcomes
##' to get an average treatment effect at dose level \eqn{T_i = t}.
##' The mean outcome at dose-level \eqn{T_i = t} is given by:
##' \deqn{
##' \hat{\mu}(t) = \frac{1}{N}\sum_i^N \hat{\alpha}_0 + \hat{\alpha}_1 t + \hat{\alpha}_2 t^2 + \hat{\alpha}_3 \hat{R}_i(t) + \hat{\alpha}_4 \hat{R^2}_i(t) + \hat{\alpha}_5 \hat{R}_i(t) \cdot t
##' }
##'
##' Different treatment levels are plugged into the previous equation
##' to estimate the missing potential outcomes.  If many \eqn{t} values
##' are evaluated, then it is possible to trace out an ADRF.
##'
##'
##'
##' @return \code{hi_est} returns an object of class "causaldrf",
##' a list that contains the following components:
##' \item{param}{parameter estimates for a hi fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{out_mod}{the result of the outcome model fit.}
##' \item{call}{the matched call.}
##'
##' @seealso \code{\link{nw_est}}, \code{\link{iw_est}}, \code{\link{hi_est}}, \code{\link{gam_est}},
##' \code{\link{add_spl_est}},
##' \code{\link{bart_est}}, etc. for other estimates.
##'
##' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
##' for use in the different estimates.
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##'
##' Hirano, Keisuke, Imbens, Guido W (2004).  The propensity score with
##' continuous treatments.  \emph{Applied Bayesian modeling and
##' causal inference from incomplete-data perspectives.}
##'
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' hi_list <- hi_est(Y = Y,
##'               treat = T,
##'               treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'               outcome_formula = Y ~ T + I(T^2) + gps + I(gps^2) + T * gps,
##'               data = example_data,
##'               grid_val = seq(8, 16, by = 1),
##'               treat_mod = "Normal")
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "hi estimate")
##'
##' lines(seq(8, 16, by = 1),
##'       hi_list$param,
##'       lty = 2,
##'       lwd = 2,
##'       col = "blue")
##'
##' legend('bottomright',
##'         "hi estimate",
##'         lty=2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y',
##'         cex=1)
##'
##' rm(example_data, hi_list, sample_index)
##'
##'
##' ## Example from van der Wal, Willem M., and Ronald B. Geskus. (2011)
##' #Simulate data with continuous confounder and outcome, binomial exposure.
##' #Marginal causal effect of exposure on outcome: 10.
##' n <- 1000
##' simdat <- data.frame(l = rnorm(n, 10, 5))
##' a.lin <- simdat$l - 10
##' pa <- exp(a.lin)/(1 + exp(a.lin))
##' simdat$a <- rbinom(n, 1, prob = pa)
##' simdat$y <- 10*simdat$a + 0.5*simdat$l + rnorm(n, -10, 5)
##' simdat[1:5,]
##' temp_hi <- hi_est(Y = y,
##'                  treat = a,
##'                  treat_formula = a ~ l,
##'                  outcome_formula = y ~ gps,
##'                  data = simdat,
##'                  grid_val = c(0, 1),
##'                  treat_mod = "Binomial",
##'                  link_function = "logit")
##'
##' temp_hi[[1]]  # estimated coefficients
##'
##'
##'
##'
##' @usage
##' hi_est(Y,
##'        treat,
##'        treat_formula,
##'        outcome_formula,
##'        data,
##'        grid_val,
##'        treat_mod,
##'        link_function,
##'        ...)
##'
##' @export
##'



hi_est <-  function (Y,
                     treat,
                     treat_formula,
                     outcome_formula,
                     data,
                     grid_val,
                     treat_mod,
                     link_function,
                     ...)
{
  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # treat_formula is the formula for the treatment model
  # data will contain all the data: X, treat, and Y
  # grid_val is the set of grid points on T
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is a list of 3 objects:
  #    (1) estimated values of ADRF at gridvalues
  #    (2) result of treatment model
  #    (3) result of the outcome model


  #save input
  tempcall <- match.call()


  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("outcome_formula" %in% names(tempcall))) stop("No outcome_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("grid_val" %in% names(tempcall)))  stop("No grid_val specified")
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) & !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Binomial", "Normal")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Binomial\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "Binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = Binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "Ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = Ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}




  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]
  )




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

  else if (treat_mod == "Binomial") {
    samp_dat <- data

    if(tempcall$link_function == "logit") lf <- binomial(link = logit)
    if(tempcall$link_function == "probit") lf  <- binomial(link = probit)
    if(tempcall$link_function == "cauchit") lf  <- binomial(link = cauchit)
    if(tempcall$link_function == "log") lf  <- binomial(link = log)
    if(tempcall$link_function == "cloglog") lf  <- binomial(link = cloglog)


      result <- glm(formula_t, family = lf, data = samp_dat)

      samp_dat$prob_1[tempdat$treat == 0] <- 1 - predict.glm(result, type = "response")[tempdat$treat == 0]
      samp_dat$prob_1[tempdat$treat == 1] <- predict.glm(result, type = "response")[tempdat$treat == 1]

      gps_fun_Binomial <- function(tt) ifelse(tt == 1, samp_dat$prob_1, 1 - samp_dat$prob_1)
      gps_fun <- gps_fun_Binomial

  }

#   else if (treat_mod == "Multinomial") {
#     samp_dat <- data
#       result1 <- nnet::multinom(
#                                   formula = formula_t,
#                                   data = samp_dat)
#       pred1 <- as.data.frame(predict(result1, type = "probs"))
#       samp_dat$prob_1 <- vector("numeric", nrow(tempdat))
#       for (i in 1:length(unique(tempdat$treat)))samp_dat$prob_1[with(tempdat, treat == sort(unique(tempdat$treat))[i])] <- pred1[tempdat$treat == sort(unique(tempdat$treat))[i],i]
#
#
#   }
#
#   else if (tempcall$treat_mod == "Ordinal") {
#     samp_dat <- data
#
#     if(tempcall$link_function == "logit") m <- "logistic"
#     if(tempcall$link_function == "probit") m  <- "probit"
#     if(tempcall$link_function == "cloglog") m  <- "cloglog"
#     if(tempcall$link_function == "cauchit") m  <- "cauchit"
#
#       result1 <- MASS::polr(
#         formula = eval(numerator_formula),
#         data = samp_dat,
#         method = m)
#       pred1 <- as.data.frame(predict(result1, type = "probs"))
#       samp_dat$prob_1  <- vector("numeric", nrow(tempdat))
#       for (i in 1:length(unique(tempdat$treat)))samp_dat$prob_1 [with(tempdat, treat == sort(unique(tempdat$treat))[i])] <- pred1[tempdat$treat == sort(unique(tempdat$treat))[i],i]
#
#
#   }
#




  else {
    stop("No valid treat_mod specified.  Please try again.")
  }


#---------------------------------

  # using the estimated gps, get the estimate for the ADRF using the hi method.

  tempdat$gps <- gps_fun(tempdat$treat)

  # this allows for the user to use flexible outcome models
colnames(tempdat) <- c(as.character(tempcall$Y), as.character(tempcall$treat), "gps" )


outcome_model <- lm(outcome_formula, data = tempdat, ...)

outcome_coef <- outcome_model$coef


# need to create model matrix with gps evaluated at grid values and then average over all the units.


mean_outcome_grid <- numeric(length(grid_val))
for (i in 1:length(grid_val)) {

  temp_matrix <- cbind(numeric(nrow(tempdat)), grid_val[i], gps_fun(grid_val[i]))
  colnames(temp_matrix) <- c(as.character(tempcall$Y), as.character(tempcall$treat), "gps" )
  temp_matrix <- data.frame(temp_matrix)

  # utils::str(m_frame <- model.frame(outcome_formula, tempdat))
  m_frame <- model.frame(outcome_formula, temp_matrix)
  covar_temp <- model.matrix(outcome_formula, m_frame)
  covar_grid_mat <- covar_temp

  potential_outcomes_at_t <- t(outcome_coef %*% t(covar_grid_mat))

  mean_outcome_grid[i] <- mean(potential_outcomes_at_t)

}




z_object <- list(param = mean_outcome_grid,
                 t_mod = result,
                 out_mod = outcome_model,
                 call = tempcall)

class(z_object) <- "causaldrf"
z_object


}
