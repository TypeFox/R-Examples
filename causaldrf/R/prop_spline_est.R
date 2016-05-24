##' The propensity-spline prediction estimator
##'
##' This method estimates the linear or quadratic parameters of the ADRF by
##' estimating a least-squares fit on the basis functions which are composed of
##' combinations of the covariates, propensity spline basis, and treatment values.
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
##' \code{treat} from \code{T_mod}.  Or, plug in gps estimates here to create
##' splines from the gps values.
##'
##' @param degree is 1 for linear and 2 for quadratic outcome model.
##' @param wt is weight used in lsfit for outcome regression.
##' Default is wt = NULL.
##' @param method is "same" if the same set of covariates are used to estimate
##' the constant, linear, and/or quadratic term with no spline terms.  If method = "different", then
##' different sets of covariates can be used to estimate the constant, linear,
##' and/or quadratic term.  To use spline terms, it is necessary to set
##' method = "different".  covar_lin_formula and covar_sq_formula must be specified
##' if method = "different".
##' @param spline_df degrees of freedom. The default, spline_df = NULL, corresponds to no knots.
##' @param spline_const is the number of spline terms to include when estimating the constant term.
##' @param spline_linear is the number of spline terms to include when estimating the linear term.
##' @param spline_quad is the number of spline terms to include when estimating the quadratic term.
##'
##' @details
##' This function estimates the ADRF by the method described in Schafer and Galagate (2015),
##' that fits an outcome model using a function of the covariates and spline basis
##' functions derived from the propensity function component.
##'
##' @return \code{prop_spline_est} returns an object of class "causaldrf_lsfit",
##' a list that contains the following components:
##' \item{param}{the estimated parameters.}
##' \item{out_mod}{the result of the outcome model fit using lsfit.}
##' \item{call}{the matched call.}
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
##' Little, Roderick and An, Hyonggin (2004).  ROBUST LIKELIHOOD-BASED ANALYSIS
##'  OF MULTIVARIATE DATA WITH MISSING VALUES. \emph{Statistica Sinica}.
##'  \bold{14}: 949--968.
##'
##' Schafer, Joseph L, Kang, Joseph (2008).  Average causal effects from
##' nonrandomized studies: a practical guide and simulated example.
##' \emph{Psychological methods}, \bold{13.4}, 279.
##'
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' t_mod_list <- t_mod(treat = T,
##'               treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'               data = example_data,
##'               treat_mod = "Normal")
##'
##' cond_exp_data <- t_mod_list$T_data
##' full_data <- cbind(example_data, cond_exp_data)
##'
##' prop_spline_list <- prop_spline_est(Y = Y,
##'                             treat = T,
##'                             covar_formula = ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'                             covar_lin_formula = ~ 1,
##'                             covar_sq_formula = ~ 1,
##'                             data = example_data,
##'                             e_treat_1 = full_data$est_treat,
##'                             degree = 1,
##'                             wt = NULL,
##'                             method = "different",
##'                             spline_df = 5,
##'                             spline_const = 4,
##'                             spline_linear = 4,
##'                             spline_quad = 4)
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "propensity spline estimate")
##'
##' abline(prop_spline_list$param[1],
##'         prop_spline_list$param[2],
##'         lty = 2,
##'         col = "blue",
##'         lwd = 2)
##'
##' legend('bottomright',
##'         "propensity spline estimate",
##'         lty = 2,
##'         bty = 'Y',
##'         cex = 1,
##'         col = "blue",
##'         lwd = 2)
##'
##' rm(example_data, prop_spline_list, sample_index)
##'
##'
##'
##' @usage
##'
##' prop_spline_est(Y,
##'                 treat,
##'                 covar_formula = ~ 1,
##'                 covar_lin_formula = ~ 1,
##'                 covar_sq_formula = ~ 1,
##'                 data,
##'                 e_treat_1 = NULL,
##'                 degree = 1,
##'                 wt = NULL,
##'                 method = "same",
##'                 spline_df = NULL,
##'                 spline_const = 1,
##'                 spline_linear = 1,
##'                 spline_quad = 1)
##'
##'
##' @export
##'
##'
##'





prop_spline_est <- function(Y,
                    treat,
                    covar_formula = ~ 1,
                    covar_lin_formula = ~ 1,
                    covar_sq_formula = ~ 1,
                    data,
                    e_treat_1 = NULL,
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



  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]

  )





  #====================================================================


  # make a formula for the treatment model
  if (method == "same"){  # this restricts the covariates for the linear and quadratic interactions to be the same.

    formula_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(covar_formula, width.cutoff = 500), sep = "")))

  # utils::str(m_frame <- model.frame(formula_covar, data))
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
    if (spline_df < 1) stop("spline_df must be a whole number greater than zero!!!  Please try again.")
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

  } else if (degree == 1 ){
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
    theta_1_est_reg <- mean((as.matrix(theta_1_covars) %*% coef_1_reg / tempdat$treat)[which(tempdat$treat != 0)] )

    reg_coefs <- c(theta_0_est_reg,
                   theta_1_est_reg)


  }else{
    stop("Error: degree needs to be 1 or 2")
  }




  z_object <- list(param = reg_coefs,
                   out_mod = ls_mod,
                   call = tempcall)

  class(z_object) <- "causaldrf_lsfit"
  z_object

}


