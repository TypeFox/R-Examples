##' The inverse probability of treatment weighting (iptw) estimator
##'
##' The iptw method or importance weighting method estimates the ADRF by
##' weighting the data with stabilized or non-stabilized weights.
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param numerator_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.  i.e.
##' \code{treat ~ 1}.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param degree is 1 for linear and 2 for quadratic outcome model.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model, \code{"Binomial"} for binomial model, \code{"Ordinal"} for ordinal model,
##' \code{"Multinomial"} for multinomial model.
##' @param link_function specifies the link function between the variables in
##' numerator or denominator and exposure, respectively.
##' For \code{treat_mod = "Gamma"} (fitted using glm) alternatives are "log" or "inverse".
##' For \code{treat_mod = "Binomial"} (fitted using glm) alternatives are "logit", "probit", "cauchit", "log" and "cloglog".
##' For \code{treat_mod = "Multinomial"} this argument is ignored, and
##' multinomial logistic regression models are always used (fitted using multinom).
##' For \code{treat_mod = "Ordinal"} (fitted using polr) alternatives are "logit", "probit", "cauchit", and "cloglog".
##'
##' @param ...	additional arguments to be passed to the low level treatment regression fitting functions.
##'
##' @details
##' This method uses inverse probability of treatment weighting to adjust
##' for possible biases.  For more details see Schafer and Galagate (2015) and
##' Robins, Hernan, and Brumback (2000).
##'
##'
##'
##'
##' @return \code{iptw_est} returns an object of class "causaldrf",
##' a list that contains the following components:
##' \item{param}{parameter estimates for a iptw fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{num_mod}{the result of the numerator model fit.}
##' \item{weights}{the estimated weights.}
##' \item{weight_data}{the weights.}
##' \item{out_mod}{the outcome model.}
##' \item{call}{the matched call.}
##'
##'
##'
##'
##' @seealso \code{\link{iptw_est}}, \code{\link{ismw_est}},
##'  \code{\link{reg_est}}, \code{\link{aipwee_est}}, \code{\link{wtrg_est}},
##'    etc. for other estimates.
##'
##' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
##' for use in the different estimates.
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' van der Wal, Willem M., and Ronald B. Geskus.
##' "IPW: an R package for inverse probability weighting."
##' \emph{Journal of Statistical Software} \bold{43.13} (2011): 1-23.
##'
##' Robins, James M and Hernan, Miguel Angel and Brumback, Babette.
##' Marginal structural models and causal inference in epidemiology.
##' \emph{Epidemiology} \bold{11.5} (2000): 550--560.
##'
##' Zhu, Yeying and Coffman, Donna L and Ghosh, Debashis.
##' A Boosting Algorithm for Estimating Generalized Propensity Scores
##' with Continuous Treatments.
##' \emph{Journal of Causal Inference} \bold{3.1} (2015): 25--40.
##'
##'
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' iptw_list <- iptw_est(Y = Y,
##'               treat = T,
##'               treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'               numerator_formula = T ~ 1,
##'               data = example_data,
##'               degree = 1,
##'               treat_mod = "Normal")
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "iptw estimate")
##'
##' abline(iptw_list$param[1],
##'         iptw_list$param[2],
##'         lty=2,
##'         lwd = 2,
##'         col = "blue")
##'
##' legend('bottomright',
##'         "iptw estimate",
##'         lty=2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y',
##'         cex=1)
##'
##' rm(example_data, iptw_list, sample_index)
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
##' temp_iptw <- iptw_est(Y = y,
##'                       treat = a,
##'                       treat_formula = a ~ l,
##'                       numerator_formula = a ~ 1,
##'                       data = simdat,
##'                       degree = 1,
##'                       treat_mod = "Binomial",
##'                       link_function = "logit")
##'
##'
##' temp_iptw[[1]]  # estimated coefficients
##'
##'
##' @usage
##'
##' iptw_est(Y,
##'          treat,
##'          treat_formula,
##'          numerator_formula,
##'          data,
##'          degree,
##'          treat_mod,
##'          link_function,
##'          ...)
##'
##' @importFrom survey svydesign svyglm
##' @export
##'


iptw_est <- function (Y,
                      treat,
                      treat_formula,
                      numerator_formula,
                      data,
                      degree,
                      treat_mod,
                      link_function,
                      ...){

  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # treat_formula is the formula for the treatment model
  # numerator_formula is the formula for the numerator
  # data will contain all the data: X, treat, and Y
  # degree is either 1 or 2
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is a list of 4 objects:
  #    (1) estimated parameters of ADRF
  #    (2) names of input values
  #    (3) result of treatment model
  #    (4) result of numerator model



  #save input
  tempcall <- match.call()


  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("numerator_formula" %in% names(tempcall))) stop("No numerator_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("degree" %in% names(tempcall)))  {if(!(tempcall$degree %in% c(1, 1))) stop("degree must be 1 or 2")}
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) &
                                             !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Normal", "Binomial", "Ordinal", "Multinomial")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Binomial\", \"Ordinal\", \"Multinomial\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "NegBinom") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = NegBinom (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "Binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "Ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}



  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]

  )


  # make a formula for the treatment model
  #   formula_t <- eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$treat_formula, width.cutoff = 500), sep = "")))
  formula_t <- tempcall$treat_formula


  if (treat_mod == "Poisson") {
    samp <- data
    result2 <- stats::glm(formula_t, family = "poisson", data = samp, ...)
    cond_mean <- result2$fitted.values
    samp$gps_vals <- dpois(x = tempdat$treat, lambda = cond_mean)
    result1 <- glm(numerator_formula, family = "poisson", data = samp, ...)
    cond_mean_num <- result1$fitted.values
    samp$num_weight <- dpois(x = tempdat$treat, lambda = cond_mean_num)
    est_import_wt <- samp$num_weight/samp$gps_vals
  }
  else if (treat_mod == "NegBinom") {
    samp <- data
    result2 <- MASS::glm.nb(formula_t, link = log, data = samp, ...)
    cond_mean <- result2$fitted.values
    cond_var <- cond_mean + cond_mean^2/result2$theta
    prob_nb_est <- (cond_var - cond_mean)/cond_var
    samp$gps_vals <- dnbinom(x = tempdat$treat, size = result2$theta,
                             mu = result2$fitted.values, log = FALSE)
    result1 <- MASS::glm.nb(numerator_formula, link = link_function,
                            data = samp, ...)
    cond_mean_num <- result1$fitted.values
    cond_var_num <- cond_mean_num + cond_mean_num^2/result1$theta
    prob_nb_est_num <- (cond_var_num - cond_mean_num)/cond_var_num
    samp$num_weight <- dnbinom(x = tempdat$treat, size = result1$theta,
                               prob = prob_nb_est_num, mu = result1$fitted.values,
                               log = FALSE)
    est_import_wt <- samp$num_weight/samp$gps_vals
  }
  else if (treat_mod == "Gamma") {
    samp <- data
    result2 <- glm(formula_t, family = Gamma(link = link_function),
                  data = samp, ...)
    shape_gamma <- as.numeric(MASS::gamma.shape(result2)[1])
    theta_given_X <- result2$fitted.values/shape_gamma
    samp$gps_vals <- dgamma(tempdat$treat, shape = shape_gamma, scale = theta_given_X)
    result1 <- glm(numerator_formula, family = Gamma(link = "log"),
                   data = samp, ...)
    shape_gamma_2 <- as.numeric(MASS::gamma.shape(result1)[1])
    theta_given_X_2 <- result1$fitted.values/shape_gamma
    samp$num_weight <- dgamma(tempdat$treat, shape = shape_gamma_2,
                              scale = theta_given_X_2)
    est_import_wt <- samp$num_weight/samp$gps_vals
  }
  else if (treat_mod == "LogNormal") {

    samp <- data
    samp[, as.character(tempcall$treat)] <- log(samp[, as.character(tempcall$treat)])
    result2 <- lm(formula_t, data = samp, ...)
    samp$gps_vals <- dnorm(samp[, as.character(tempcall$treat)], mean = result2$fitted,
                           sd = summary(result2)$sigma)

    result1 <- lm(numerator_formula, data = samp, ...)
    samp$num_weight <- dnorm(samp[, as.character(tempcall$treat)], mean = result1$fitted,
                             sd = summary(result1)$sigma)
    est_import_wt <- samp$num_weight/samp$gps_vals
  }
  else if (treat_mod == "Sqrt") {

    samp <- data
    samp[, as.character(tempcall$treat)] <- sqrt(samp[, as.character(tempcall$treat)])
    result2 <- lm(formula_t, data = samp, ...)
    samp$gps_vals <- dnorm(samp[, as.character(tempcall$treat)], mean = result2$fitted,
                           sd = summary(result2)$sigma)

    result1 <- lm(numerator_formula, data = samp, ...)
    samp$num_weight <- dnorm(samp[, as.character(tempcall$treat)], mean = result1$fitted,
                             sd = summary(result1)$sigma)
    est_import_wt <- samp$num_weight/samp$gps_vals
  }
  else if (treat_mod == "Normal") {
    samp <- data
    result2 <- lm(formula_t, data = samp, ...)
    samp$gps_vals <- dnorm(tempdat$treat, mean = result2$fitted, sd = summary(result2)$sigma)
    result1 <- lm(numerator_formula, data = samp, ...)
    samp$num_weight <- dnorm(tempdat$treat, mean = result1$fitted,
                             sd = summary(result1)$sigma)
    est_import_wt <- samp$num_weight/samp$gps_vals
  }

  else if (treat_mod == "Binomial") {
    samp <- data

    if(tempcall$link_function == "logit") lf <- binomial(link = logit)
    if(tempcall$link_function == "probit") lf  <- binomial(link = probit)
    if(tempcall$link_function == "cauchit") lf  <- binomial(link = cauchit)
    if(tempcall$link_function == "log") lf  <- binomial(link = log)
    if(tempcall$link_function == "cloglog") lf  <- binomial(link = cloglog)
    if (is.null(tempcall$numerator_formula)) samp$num_weight <- 1
    else {

      result1 <- glm(numerator_formula, family = lf, data = samp, ...)

      samp$num_weight[tempdat$treat == 0] <- 1 - predict.glm(result1, type = "response")[tempdat$treat == 0]
      samp$num_weight[tempdat$treat == 1] <- predict.glm(result1, type = "response")[tempdat$treat == 1]
      result1$call$formula <- eval(numerator_formula)
      result1$call$family <- tempcall$link_function
      result1$call$data <- eval(tempcall$data)
    }
    result2 <- glm(formula_t, family = lf, data = samp, ...)

    samp$denom_weight[tempdat$treat == 0] <- 1 - predict.glm(result2, type = "response")[tempdat$treat == 0]
    samp$denom_weight[tempdat$treat == 1] <- predict.glm(result2, type = "response")[tempdat$treat == 1]
    result2$call$formula <- eval(formula_t)
    result2$call$family <- tempcall$link_function
    result2$call$data <- eval(tempcall$data)
    est_import_wt <- samp$num_weight/samp$denom_weight
  }

  else if (treat_mod == "Multinomial") {
    samp <- data
    if (is.null(tempcall$numerator_formula)) samp$num_weight <- 1
    else {
      result1 <- nnet::multinom(
                          formula = numerator_formula,
                          data = samp, ...)
      pred1 <- as.data.frame(predict(result1, type = "probs"))
      samp$num_weight <- vector("numeric", nrow(tempdat))
      for (i in 1:length(unique(tempdat$treat)))samp$num_weight[with(tempdat, treat == sort(unique(tempdat$treat))[i])] <- pred1[tempdat$treat == sort(unique(tempdat$treat))[i],i]
      result1$call$formula <- eval(numerator_formula)
      result1$call$data <- eval(tempcall$data)
    }
    result2 <- nnet::multinom(
                      formula = formula_t,
                      data = samp, ...)

    pred2 <- as.data.frame(predict(result2, type = "probs"))
    samp$denom_weight <- vector("numeric", nrow(tempdat))
    for (i in 1:length(unique(tempdat$treat)))samp$denom_weight[with(tempdat, treat == sort(unique(tempdat$treat))[i])] <- pred2[tempdat$treat == sort(unique(tempdat$treat))[i],i]
    result2$call$formula <- eval(formula_t)
    result2$call$data <- eval(tempcall$data)
    est_import_wt <- samp$num_weight/samp$denom_weight
  }

  else if (tempcall$treat_mod == "Ordinal") {
    if(tempcall$link_function == "logit") m <- "logistic"
    if(tempcall$link_function == "probit") m  <- "probit"
    if(tempcall$link_function == "cloglog") m  <- "cloglog"
    if(tempcall$link_function == "cauchit") m  <- "cauchit"
    if (is.null(tempcall$numerator_formula)) samp$num_weight <- 1
    else {
      result1 <- MASS::polr(
        formula = eval(numerator_formula),
        data = samp,
        method = m, ...)
      pred1 <- as.data.frame(predict(result1, type = "probs"))
      samp$num_weight  <- vector("numeric", nrow(tempdat))
      for (i in 1:length(unique(tempdat$treat)))samp$num_weight [with(tempdat, treat == sort(unique(tempdat$treat))[i])] <- pred1[tempdat$treat == sort(unique(tempdat$treat))[i],i]
      result1$call$formula <- eval(numerator_formula)
      result1$call$data <- eval(tempcall$data)
      result1$call$method <- m
    }
    result2 <- MASS::polr(
      formula = eval(formula_t),
      data = samp,
      method = m,
      na.action = na.fail, ...)
    pred2 <- as.data.frame(predict(result2, type = "probs"))
    samp$denom_weight <- vector("numeric", nrow(tempdat))
    for (i in 1:length(unique(tempdat$treat)))samp$denom_weight[with(tempdat, treat == sort(unique(tempdat$treat))[i])] <- pred2[tempdat$treat == sort(unique(tempdat$treat))[i],i]
    result2$call$formula <- eval(formula_t)
    result2$call$data <- eval(tempcall$data)
    result2$call$method <- m
    est_import_wt <- samp$num_weight/samp$denom_weight
  }


  else {
    stop("Treatment model specified is not valid.  Please try again.")
  }
  if (degree == 2) {
    weight_dat <- data.frame(tempdat$Y, tempdat$treat, tempdat$treat^2, est_import_wt)
    colnames(weight_dat) <- c("Y", "treat", "treat.2", "est_import_wt")
    design.weight <- survey::svydesign(ids = ~1, weights = ~est_import_wt,
                                       data = weight_dat)
    import_weight_mod <- survey::svyglm(Y ~ treat + treat.2,
                                        design = design.weight)
    se_est_iptw_coef <- sqrt(diag(vcov(import_weight_mod)))
    import_coefs <- coef(import_weight_mod)
  }
  else if (degree == 1) {
    weight_dat <- data.frame(tempdat$Y, tempdat$treat, est_import_wt)
    colnames(weight_dat) <- c("Y", "treat", "est_import_wt")

    design.weight <- survey::svydesign(ids = ~1, weights = ~est_import_wt,
                                       data = weight_dat)
    import_weight_mod <- survey::svyglm(Y ~ treat, design = design.weight)
    se_est_iptw_coef <- sqrt(diag(vcov(import_weight_mod)))
    import_coefs <- coef(import_weight_mod)
  }
  else {
    stop("Error: degree needs to be 1 or 2")
  }




  z_object <- list(param = import_coefs,
                   t_mod = result2,
                   num_mod = result1,
                   weights = est_import_wt,
                   weight_data = weight_dat,
                   out_mod = import_weight_mod,
                   call = tempcall)

  class(z_object) <- "causaldrf"
  z_object


}
