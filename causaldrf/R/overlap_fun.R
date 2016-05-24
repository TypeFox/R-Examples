##' This function creates an overlapping dataset
##'
##' This function ensures that the units overlap according to the estimated gps
##' values. The overlapping dataset depends on the number of classes
##' \code{n_class} to subclassify on.
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param data_set is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param n_class is the number of classes to split \code{gps} into.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model,  \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model.
##' @param link_function is either "log", "inverse", or "identity" for the
##' "Gamma" \code{treat_mod}.
##' @param ... additional arguments to be passed to the treatment regression function
##'
##' @return \code{overlap_fun} returns a list containing the following
##' elements: \item{overlap_dataset}{dataframe containing overlapping data.}
##' \item{median_vec}{a vector containing median values.}
##' \item{overlap_treat_result}{the resulting treatment fit.}
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
##' Bia, Michela, et al.
##' "A Stata package for the application of semiparametric estimators of dose response functions."
##'  \emph{Stata Journal} \bold{14.3} (2014): 580-604.
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' overlap_list <- overlap_fun(Y = Y,
##'                   treat = T,
##'                   treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'                   data_set = example_data,
##'                   n_class = 3,
##'                   treat_mod = "Normal")
##'
##' overlapped_data <- overlap_list$overlap_dataset
##' summary(overlapped_data)
##'
##' rm(example_data, overlap_list, overlapped_data)
##'
##'
##'
##' @usage
##'
##' overlap_fun(Y,
##'             treat,
##'             treat_formula,
##'             data_set,
##'             n_class,
##'             treat_mod,
##'             link_function,
##'             ...)
##'
##' @export
##'
##'





overlap_fun <- function(Y,
                        treat,
                        treat_formula,
                        data_set,
                        n_class,
                        treat_mod,
                        link_function,
                        ...){

  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # treat_formula is the formula for the treatment model
  # data_set will contain all the data_set: X, treat, and Y
  # n_class is the number of classes for which to use in the overlap
  # treat_mod is th treatment model to fit
  # link_function is the link function used, if needed

  # The outcome is a list of 2 objects:
  #    (1) The overlapping dataset
  #    (2) a vector containing median values
  #    (3) the resulting treatment fit
  #    (4) the cutoff gps values of each subclass



  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("treat_formula" %in% names(tempcall))) stop("No treat_formula model specified")
  if (!("data_set" %in% names(tempcall))) stop("No data_set specified")
  if (!("n_class" %in% names(tempcall))) stop("No n_class specified")
  if (!("treat_mod" %in% names(tempcall)) | ("treat_mod" %in% names(tempcall) & !(tempcall$treat_mod %in% c("NegBinom", "Poisson", "Gamma", "LogNormal", "Sqrt", "Normal")))) stop("No valid family specified (\"NegBinom\", \"Poisson\", \"Gamma\", \"Log\", \"Sqrt\", \"Normal\")")
  if (tempcall$treat_mod == "Gamma") {if(!(tempcall$link_function %in% c("log", "inverse"))) stop("No valid link function specified for family = Gamma (\"log\", \"inverse\")")}
  if (tempcall$treat_mod == "binomial") {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$treat_mod == "ordinal" ) {if(!(tempcall$link_function %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = ordinal (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}

  sample_dataset <- eval(tempcall$data_set)

  full_dataset_ord <- sample_dataset[order(sample_dataset[, as.character(tempcall$treat)]), ]
  treat <- full_dataset_ord[, as.character(tempcall$treat)]

  ### Use n_div_groups for dividing groups.
  n_div_groups <- n_class
  n_sequence <- 1:nrow(full_dataset_ord)

  support_indices <- list(n_div_groups)
  for (i in 1:n_div_groups){
    support_indices[[i]] <- which(cut(n_sequence, breaks = floor(quantile(n_sequence, probs=seq(0,1, by=1/n_div_groups))),
                                      labels = FALSE,
                                      include.lowest = TRUE) == i)
  }

  full_dataset_ord$support_indices <- cut(n_sequence, breaks = floor(quantile(n_sequence, probs=seq(0,1, by=1/n_div_groups))),
                                          labels = FALSE,
                                          include.lowest = TRUE)

  CS_list <- list(numeric(n_div_groups))
  median_list <- numeric(n_div_groups)

  treat_val <- full_dataset_ord[, as.character(tempcall$treat)]


  #-------------------------------------------------------
  #-------------------------------------------------------
  #-------------------------------------------------------


  if (treat_mod == "NegBinom"){


    result <- MASS::glm.nb(treat_formula,
                           data = full_dataset_ord,
                           ...)

    cond_mean <- result$fitted.values
    cond_var <- cond_mean + cond_mean^2/result$theta
    prob_nb_est <- (cond_var - cond_mean) / cond_var

    gps_fun_NB <- function(tt) {dnbinom(x = tt,
                                        size = result$theta,
                                        mu = result$fitted.values,
                                        log = FALSE)}

    gps_fun <- gps_fun_NB


  } else if (treat_mod == "Poisson"){


    result <- glm(treat_formula,
                  family = "poisson",
                  data = full_dataset_ord,
                  ...)

    cond_mean <- result$fitted.values

    samp_dat$gps_vals <- dpois(x = treat,
                               lambda = cond_mean)


    gps_fun_Pois <- function(t){  dpois(t,
                                        lambda = cond_mean) }

    gps_fun <- gps_fun_Pois

  } else if (treat_mod == "Gamma"){



    # this seems to be the best simple model with main effects
    result <- glm(treat_formula,
                  family = Gamma(link = link_function),
                  data = full_dataset_ord,
                  ...)

    est_treat <- result$fitted


    #     shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    #     theta_given_X <- result$fitted.values/shape_gamma
    #
    #     sigma_est <-  sqrt(sum((result$residuals)^2)/ df.residual(result))

    shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    theta_given_X <- result$fitted.values/shape_gamma
    theta_treat_X <-treat/shape_gamma


    gps_fun_Gamma <- function(t){ dgamma(t,
                                         shape = shape_gamma,
                                         scale = theta_given_X)  }

    gps_fun <- gps_fun_Gamma

  } else if (treat_mod == "LogNormal"){
    full_dataset_ord[, as.character(tempcall$treat)] <- log(treat)



    result <- lm(treat_formula,
                 data = full_dataset_ord,
                 ...)

    est_log_treat <- result$fitted
    sigma_est <- summary(result)$sigma

    gps_fun_Log <- function(tt) {dnorm(log(tt), mean = est_log_treat, sd = sigma_est )}

    gps_fun <- gps_fun_Log


  } else if (treat_mod == "Sqrt"){
    full_dataset_ord[, as.character(tempcall$treat)] <- sqrt(treat)



    result <- lm(treat_formula,
                 data = full_dataset_ord,
                 ...)

    est_sqrt_treat <- result$fitted
    sigma_est <- summary(result)$sigma

    gps_fun_sqrt <- function(tt) {dnorm(sqrt(tt), mean = est_sqrt_treat, sd = sigma_est )}

    gps_fun <- gps_fun_sqrt


  } else if (treat_mod == "Normal"){

    result <- lm(treat_formula,
                 data = full_dataset_ord,
                 ...)

    gps_fun_Normal <- function(tt) {dnorm(tt, mean = result$fitted, sd = summary(result)$sigma )}

    gps_fun <- gps_fun_Normal


  }


  # else {print("Treatment model specified is not valid.  Please try again.")}



  #-------------------------------------------------------
  #-------------------------------------------------------
  #-------------------------------------------------------

  for (i in 1:n_div_groups) {

    med_1 <- median(full_dataset_ord[support_indices[[i]], as.character(tempcall$treat)])

    R_hat <- gps_fun(med_1)
    left_cut <- max(min(R_hat[support_indices[[i]]]), min(R_hat[-support_indices[[i]]]))
    right_cut <- min(max(R_hat[support_indices[[i]]]), max(R_hat[-support_indices[[i]]]))

    CS_list[[i]] <- which((R_hat >= left_cut) & (R_hat <= right_cut))
    median_list[i] <- med_1
  }

  full_dataset_ord[, as.character(tempcall$treat)] <- treat


  CS_final <- Reduce(intersect, CS_list)
  length(intersect(intersect(CS_list[[1]], CS_list[[2]]), CS_list[[3]]))
  full_dataset_ord_1 <- full_dataset_ord[CS_final, ]


  return(list(overlap_dataset = full_dataset_ord_1, median_vec = median_list, overlap_treat_result = result))

}

