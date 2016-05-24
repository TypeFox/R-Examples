##' Prediction with a residual bias correction estimator
##'
##' This method combines the regression estimator with a residual bias correction
##' for estimating a parametric ADRF.
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param covar_formula is the formula to describe the covariates needed
##' to estimate the constant term:
##' \code{~ X.1 + ....}. Can include higher order terms or interactions.  i.e.
##' \code{~ X.1 + I(X.1^2) + X.1 * X.2 + ....}.  Don't forget the tilde before
##' listing the covariates.
##'
##' @param covar_lin_formula is the formula to describe the covariates needed
##' to estimate the linear term, t:
##' \code{~ X.1 + ....}. Can include higher order terms or interactions.  i.e.
##' \code{~ X.1 + I(X.1^2) + X.1 * X.2 + ....}.  Don't forget the tilde before
##' listing the covariates.
##'
##' @param covar_sq_formula is the formula to describe the covariates needed
##' to estimate the quadratic term, t^2:
##' \code{~ X.1 + ....}. Can include higher order terms or interactions.  i.e.
##' \code{~ X.1 + I(X.1^2) + X.1 * X.2 + ....}.  Don't forget the tilde before
##' listing the covariates.
##'
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##'
##' @param e_treat_1 a vector, representing the conditional expectation of
##' \code{treat} from \code{T_mod}.
##' @param e_treat_2 a vector, representing the conditional expectation of
##' \code{treat^2} from \code{T_mod}.
##' @param e_treat_3 a vector, representing the conditional expectation of
##' \code{treat^3} from \code{T_mod}.
##' @param e_treat_4 a vector, representing the conditional expectation of
##' \code{treat^4} from \code{T_mod}.
##'
##' @param degree is 1 for linear and 2 for quadratic outcome model.
##' @param wt is weight used in lsfit for outcome regression.
##' Default is wt = NULL.
##' @param method is "same" if the same set of covariates are used to estimate
##' the constant, linear, and/or quadratic term.  If method = "different", then
##' different sets of covariates can be used to estimate the constant, linear,
##' and/or quadratic term.  covar_lin_formula and covar_sq_formula must be specified
##' if method = "different".
##' @param spline_df degrees of freedom. The default, spline_df = NULL, corresponds to no knots.
##' @param spline_const is the number of spline terms needed to estimate the constant term.
##' @param spline_linear is the number of spline terms needed to estimate the linear term.
##' @param spline_quad is the number of spline terms needed to estimate the quadratic term.
##'
##' @details
##'
##'   This estimator bears a strong
##'   resemblance to general regression estimators in the survey
##'   literature, part of a more general class of calibration
##'   estimators (Deville and Sarndal, 1992). It is
##'   doubly robust, which means that it is consistent if
##'   either of the models is true (Scharfstein, Rotnitzky and Robins
##'   1999).  If the Y-model is correct, then the first term in
##'   the previous equation is unbiased for \eqn{\xi} and the second term has mean
##'   zero even if the T-model is wrong. If the Y-model is incorrect, the
##'   first term is biased, but the second term gives a consistent estimate
##'   of (minus one times) the bias from the Y-model if the T-model is
##'   correct.
##'
##' This function is a doubly-robust estimator that fits an outcome regression
##' model with a bias correction term.  For details see Schafer and Galagate (2015).
##'
##' @return \code{aipwee_est}  returns an object of class "causaldrf_lsfit",
##' a list that contains the following components:
##' \item{param}{parameter estimates for a add_spl fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{out_mod}{the result of the outcome model fit.}
##' \item{call}{the matched call.}
##'
##' @seealso \code{\link{iptw_est}}, \code{\link{ismw_est}},
##'  \code{\link{reg_est}}, \code{\link{wtrg_est}},
##'   ##'    etc. for other estimates.
##'
##' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
##' for use in the different estimates.
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' Schafer, Joseph L, Kang, Joseph (2008).  Average causal effects from
##' nonrandomized studies: a practical guide and simulated example.
##' \emph{Psychological methods}, \bold{13.4}, 279.
##'
##' Robins, James M and Rotnitzky, Andrea (1995).
##' Semiparametric efficiency in multivariate regression models with missing data
##' \emph{Journal of the American Statistical Association}, \bold{90.429}, 122--129.
##'
##' Scharfstein, Daniel O and Rotnitzky, Andrea and Robins, James M (1999).
##' Adjusting for nonignorable drop-out using semiparametric nonresponse models
##' \emph{Journal of the American Statistical Association}, \bold{94.448}, 1096--1120.
##'
##' Deville, Jean-Claude and Sarndal, Carl-Erik (1992).
##' Calibration estimators in survey sampling
##' \emph{Journal of the American Statistical Association}, \bold{87.418}, 376--380.
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##'
##' t_mod_list <- t_mod(treat = T,
##'               treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'               data = example_data,
##'               treat_mod = "Normal")
##'
##' cond_exp_data <- t_mod_list$T_data
##' full_data <- cbind(example_data, cond_exp_data)
##'
##' aipwee_list <- aipwee_est(Y = Y,
##'                          treat = T,
##'                          covar_formula = ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'                          covar_lin_formula = ~ 1,
##'                          covar_sq_formula = ~ 1,
##'                          data = example_data,
##'                          e_treat_1 = full_data$est_treat,
##'                          e_treat_2 = full_data$est_treat_sq,
##'                          e_treat_3 = full_data$est_treat_cube,
##'                          e_treat_4 = full_data$est_treat_quartic,
##'                          degree = 1,
##'                          wt = NULL,
##'                          method = "same",
##'                          spline_df = NULL,
##'                          spline_const = 1,
##'                          spline_linear = 1,
##'                          spline_quad = 1)
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "aipwee estimate")
##'
##' abline(aipwee_list$param[1],
##'         aipwee_list$param[2],
##'         lty = 2,
##'         lwd = 2,
##'         col = "blue")
##'
##' legend('bottomright',
##'         "aipwee estimate",
##'         lty = 2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y',
##'         cex=1)
##'
##' rm(example_data, t_mod_list, cond_exp_data, full_data, aipwee_list, sample_index)
##'
##'
##'
##' @usage
##'
##' aipwee_est(Y,
##'            treat,
##'            covar_formula = ~ 1,
##'            covar_lin_formula = ~ 1,
##'            covar_sq_formula = ~ 1,
##'            data,
##'            e_treat_1 = NULL,
##'            e_treat_2 = NULL,
##'            e_treat_3 = NULL,
##'            e_treat_4 = NULL,
##'            degree = 1,
##'            wt = NULL,
##'            method = "same",
##'            spline_df = NULL,
##'            spline_const = 1,
##'            spline_linear = 1,
##'            spline_quad = 1)
##'
##'
##' @export
##'
##'
##'





aipwee_est <- function(Y,
                       treat,
                       covar_formula = ~ 1,
                       covar_lin_formula = ~ 1,
                       covar_sq_formula = ~ 1,
                       data,
                       e_treat_1 = NULL,
                       e_treat_2 = NULL,
                       e_treat_3 = NULL,
                       e_treat_4 = NULL,
                       degree = 1,
                       wt = NULL,
                       method = "same",
                       spline_df = NULL,
                       spline_const = 1,
                       spline_linear = 1,
                       spline_quad = 1){

  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # covar_formula is the formula for the covariates model of the form: ~ X.1 + ....
  # data will contain all the data: X, treat, and Y
  # degree is 1 for linear and 2 for quadratic outcome model

  # The outcome is the estimated parameters.


  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("degree" %in% names(tempcall)))  {if(!(tempcall$degree %in% c(1, 1))) stop("degree must be 1 or 2")}
  if (!("e_treat_1" %in% names(tempcall))) stop("No e_treat_1 specified")
  if (!("e_treat_2" %in% names(tempcall))) stop("No e_treat_2 specified")
  if (!("e_treat_3" %in% names(tempcall))) stop("No e_treat_3 specified")
  if (!("e_treat_4" %in% names(tempcall))) stop("No e_treat_4 specified")



  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]

  )


  #====================================================================


  # make a formula for the treatment model
  if (method == "same"){  # this restricts the covariates for the linear and quadratic interactions to be the same.

    formula_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(covar_formula, width.cutoff = 500), sep = "")))

    m_frame <- model.frame(formula_covar, data)
    covar_sq_mat <- covar_lin_mat <- covar_mat <- model.matrix(formula_covar, m_frame)

  } else if (method == "different" & is.null(spline_df) ){  # this allows the covariates for the constant, linear, and quadratic interactions to be different.

    formula_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_formula, width.cutoff = 500), sep = "")))
    formula_lin_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_lin_formula, width.cutoff = 500), sep = "")))
    formula_sq_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_sq_formula, width.cutoff = 500), sep = "")))

    m_frame <- model.frame(formula_covar, data)
    covar_mat <- model.matrix(formula_covar, m_frame)  # matrix containing covariates needed in constant term parameter.

    m_frame_lin <- model.frame(formula_lin_covar, data)
    covar_lin_mat <- model.matrix(formula_lin_covar, m_frame_lin)  # matrix containing covariates needed in linear term parameter.

    m_frame_sq <- model.frame(formula_sq_covar, data)
    covar_sq_mat <- model.matrix(formula_sq_covar, m_frame_sq)  # matrix containing covariates needed in square term parameter.

  } else if (method == "different" & is.numeric(spline_df)){  # this allows the covariates for the constant, linear, and quadratic interactions to be different and includes spline terms.
    spline_df <- round(spline_df)
    if (spline_df < 1) stop("spline_df must be an integer greater than zero!!!  Please try again.")
    if (is.null(e_treat_1)) stop("No e_treat_1 specified!!!  Please try again.")


    #====================================================================
    # create spline basis matrix for a natural cubic spline
    spline_basis <- splines::ns(e_treat_1,
                                df = spline_df,
                                intercept = FALSE,
                                Boundary.knots = range(e_treat_1))


    colnames(spline_basis) <- paste("b",
                                    1:ncol(spline_basis),
                                    sep="")



    # add spline terms to constant, linear, and quadratic interactions.
    const_spline <-  spline_basis[, 1:spline_const]
    linear_spline <-  spline_basis[, 1:spline_linear]
    quad_spline  <-  spline_basis[, 1:spline_quad]


    formula_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_formula, width.cutoff = 500), sep = "")))
    formula_lin_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_lin_formula, width.cutoff = 500), sep = "")))
    formula_sq_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_sq_formula, width.cutoff = 500), sep = "")))

    m_frame <- model.frame(formula_covar, data)
    covar_mat <- model.matrix(formula_covar, m_frame)  # matrix containing covariates needed in constant term parameter.
    covar_mat <- cbind(covar_mat, const_spline)

    m_frame_lin <- model.frame(formula_lin_covar, data)
    covar_lin_mat <- model.matrix(formula_lin_covar, m_frame_lin)  # matrix containing covariates needed in linear term parameter.
    covar_lin_mat <- cbind(covar_lin_mat, linear_spline)

    m_frame_sq <- model.frame(formula_sq_covar, data)
    covar_sq_mat <- model.matrix(formula_sq_covar, m_frame_sq)  # matrix containing covariates needed in square term parameter.
    covar_sq_mat <- cbind(covar_sq_mat, quad_spline)


  } else {
    stop("Error: method is not recognized.")
  }


  if (degree == 2){
    if (is.null(e_treat_1)) stop("No e_treat_1 specified!!!  Please try again.")
    if (is.null(e_treat_2)) stop("No e_treat_2 specified!!!  Please try again.")
    if (is.null(e_treat_3)) stop("No e_treat_3 specified!!!  Please try again.")
    if (is.null(e_treat_4)) stop("No e_treat_4 specified!!!  Please try again.")

    B <- cbind(1,
               tempdat$treat,
               tempdat$treat^2)

    n_B <- ncol(B)

    E_mat <- matrix(numeric(n_B * n_B),
                    nrow = n_B)

    wBB <- matrix(numeric(n_B * n_B),
                  nrow = n_B)

    # create the covariate matrix that will have:
    # intercept, main effects, treat, treat * (main effects), treat^2, treat^2 * (main effects)

    basis_mat <- cbind(covar_mat,
                       tempdat$treat * covar_lin_mat,
                       tempdat$treat^2 * covar_sq_mat)

    colnames(basis_mat) <- paste("x",
                                 1:ncol(basis_mat),
                                 sep="")

    ls_mod <- lsfit(basis_mat,
                    tempdat$Y,
                    wt,
                    intercept = FALSE)

    coef_reg_pred <- ls_mod$coefficients

    n_covars <- ncol(covar_mat)
    n_covars_lin <- ncol(covar_lin_mat)
    n_covars_sq <- ncol(covar_sq_mat)

    coef_0_reg <- coef_reg_pred[ c(1:n_covars) ]
    coef_1_reg <- coef_reg_pred[ c( (n_covars + 1):(n_covars + n_covars_lin) ) ]
    coef_2_reg <- coef_reg_pred[ c( (n_covars + n_covars_lin + 1):(n_covars + n_covars_lin + n_covars_sq) ) ]

    theta_0_covars <- basis_mat[,  c(1:n_covars) ]
    theta_1_covars <- basis_mat[,  c( (n_covars + 1):(n_covars + n_covars_lin) )   ]
    theta_2_covars <- basis_mat[,  c( (n_covars + n_covars_lin + 1):(n_covars + n_covars_lin + n_covars_sq) )  ]


    theta_0_est_reg <- mean(as.matrix(theta_0_covars) %*% coef_0_reg )
    theta_1_est_reg <- mean( (as.matrix(theta_1_covars) %*% coef_1_reg / tempdat$treat)[which(tempdat$treat != 0)] )
    theta_2_est_reg <- mean( (as.matrix(theta_2_covars) %*% coef_2_reg / tempdat$treat^2 )[which(tempdat$treat != 0)] )

    reg_coefs <- c(theta_0_est_reg,
                   theta_1_est_reg,
                   theta_2_est_reg)

    Y_hat <- theta_0_est_reg + theta_1_est_reg * tempdat$treat +
             theta_2_est_reg * tempdat$treat^2


    E_mat <- matrix(numeric(n_B * n_B),
                    nrow = n_B)

    wBB_3 <- matrix(numeric(n_B),
                    nrow = n_B)

    for (i in 1:nrow(covar_mat) ) {
      E_mat[1, 1] <- 1
      E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
      E_mat[2, 2] <- E_mat[1, 3] <- E_mat[3, 1] <-  e_treat_2[i]
      E_mat[2, 3] <- E_mat[3, 2] <- e_treat_3[i]
      E_mat[3, 3] <- e_treat_4[i]

      wBB_3 <- wBB_3 + solve(E_mat)%*% as.matrix(B[i,]) %*%(tempdat$Y[i] - Y_hat[i])
    }


    aipwee_coefs <- reg_coefs + (1 / nrow(covar_mat)) * wBB_3

  } else if (degree == 1 ){
    if (is.null(e_treat_1)) stop("No e_treat_1 specified!!!  Please try again.")
    if (is.null(e_treat_2)) stop("No e_treat_2 specified!!!  Please try again.")


    B <- cbind(1,
               tempdat$treat)

    n_B <- ncol(B)

    E_mat <- matrix(numeric(n_B * n_B),
                    nrow = n_B)

    wBB <- matrix(numeric(n_B * n_B),
                  nrow = n_B)

    # create the covariate matrix that will have:
    # intercept, main effects, treat, treat * (main effects)

    basis_mat <- cbind(covar_mat,
                       tempdat$treat * covar_lin_mat)

    colnames(basis_mat) <- paste("x",
                                 1:ncol(basis_mat),
                                 sep="")

    ls_mod <- lsfit(basis_mat,
                    tempdat$Y,
                    wt,
                    intercept = FALSE)

    coef_reg_pred <- ls_mod$coefficients

    n_covars <- ncol(covar_mat)
    n_covars_lin <- ncol(covar_lin_mat)

    coef_0_reg <- coef_reg_pred[ c(1:n_covars) ]
    coef_1_reg <- coef_reg_pred[ c( (n_covars + 1):(n_covars + n_covars_lin) ) ]

    theta_0_covars <- basis_mat[,  c(1:n_covars) ]
    theta_1_covars <- basis_mat[,  c( (n_covars + 1):(n_covars + n_covars_lin) )   ]

    theta_0_est_reg <- mean(as.matrix(theta_0_covars) %*% coef_0_reg )
    theta_1_est_reg <- mean( (as.matrix(theta_1_covars) %*% coef_1_reg / tempdat$treat )[which(tempdat$treat != 0)] )

    reg_coefs <- c(theta_0_est_reg,
                   theta_1_est_reg)


    Y_hat <- theta_0_est_reg + theta_1_est_reg * tempdat$treat

    E_mat <- matrix(numeric(n_B * n_B),
                    nrow = n_B)

    wBB_3 <- matrix(numeric(n_B),
                    nrow = n_B)

    for (i in 1:nrow(covar_mat) ) {
      E_mat[1, 1] <- 1
      E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
      E_mat[2, 2] <- e_treat_2[i]

      wBB_3 <- wBB_3 + solve(E_mat)%*% as.matrix(B[i,]) %*%(tempdat$Y[i] - Y_hat[i])
    }


    aipwee_coefs <- reg_coefs + (1 / nrow(covar_mat)) * wBB_3



  }else{
    stop("Error: degree needs to be 1 or 2")
  }




  z_object <- list(param = aipwee_coefs,
                   out_mod = ls_mod,
                   reg_param = reg_coefs,
                   call = tempcall)

  class(z_object) <- "causaldrf_lsfit"
  z_object


}


