#' The weighted regression estimator
#'
#' This method uses weight matrices to
#' estimate parameters for an ADRF with quadratic or linear fits.
#'
#'
#' @param Y is the output
#' @param treat is the treatment variable
#' @param covar_formula is the formula for the covariates model of the form: ~ X.1 + ....
#' @param data will contain all the data: X, treat, and Y
#' @param e_treat_1 is estimated treatment
#' @param e_treat_2 is estimated treatment squared
#' @param e_treat_3 is estimated treatment cubed
#' @param e_treat_4 is estimated treatment to the fourth
#' @param degree is 1 for linear fit and 2 for quadratic fit
#'
#' @details
#' This function estimates the ADRF by the method described in Schafer and Galagate (2015)
#' which uses weight matrices to adjust for possible bias.
#'
#' @return \code{wtrg_est} returns an object of class "causaldrf",
#' a list that contains the following components:
#' \item{param}{the estimated parameters.}
#' \item{call}{the matched call.}
#'
#' @seealso \code{\link{iptw_est}}, \code{\link{ismw_est}},
#'  \code{\link{reg_est}}, \code{\link{aipwee_est}}, \code{\link{wtrg_est}},
#'    etc. for other estimates.
#'
#' \code{\link{t_mod}}, \code{\link{overlap_fun}} to prepare the \code{data}
#' for use in the different estimates.
#'
#'
#' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
#' continuous treatment and outcome: alternative estimators for parametric
#' dose-response models. \emph{Manuscript in preparation}.
#'
#'
#' @examples
#'
#' ## Example from Schafer (2015).
#'
#' example_data <- sim_data
#'
#'
#' t_mod_list <- t_mod(treat = T,
#'               treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
#'               data = example_data,
#'               treat_mod = "Normal")
#'
#' cond_exp_data <- t_mod_list$T_data
#' full_data <- cbind(example_data, cond_exp_data)
#'
#' wtrg_list <- wtrg_est(Y = Y,
#'                       treat = T,
#'                       covar_formula = ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
#'                       data = example_data,
#'                       e_treat_1 = full_data$est_treat,
#'                       e_treat_2 = full_data$est_treat_sq,
#'                       e_treat_3 = full_data$est_treat_cube,
#'                       e_treat_4 = full_data$est_treat_quartic,
#'                       degree = 1)
#'
#' sample_index <- sample(1:1000, 100)
#'
#' plot(example_data$T[sample_index],
#'       example_data$Y[sample_index],
#'       xlab = "T",
#'       ylab = "Y",
#'       main = "weighted regression estimate")
#'
#' abline(wtrg_list$param[1],
#'         wtrg_list$param[2],
#'         lty = 2,
#'         lwd = 2,
#'         col = "blue")
#'
#' legend('bottomright',
#'         "weighted regression estimate",
#'         lty = 2,
#'         lwd = 2,
#'         col = "blue",
#'         bty='Y',
#'         cex=1)
#'
#' rm(example_data, t_mod_list, cond_exp_data, full_data, wtrg_list, sample_index)
#'
#'
#' @usage
#' wtrg_est(Y,
#'          treat,
#'          covar_formula,
#'          data,
#'          e_treat_1,
#'          e_treat_2,
#'          e_treat_3,
#'          e_treat_4,
#'          degree)
#'
#' @export
#'
#'



wtrg_est <-     function(Y,
                         treat,
                         covar_formula,
                         data,
                         e_treat_1,
                         e_treat_2,
                         e_treat_3,
                         e_treat_4,
                         degree){


  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # covar_formula is the formula for the covariates model of the form: ~ X.1 + ....
  # data will contain all the data: X, treat, and Y
  # e_treat_1 is the vector containing the estimated means
  # e_treat_2 is the vector containing the estimated second moments
  # e_treat_3 is the vector containing the estimated third moments
  # e_treat_4 is the vector containing the estimated fourth moments
  # degree is either 1 or 2

  # The outcome is the estimated parameters.


  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("covar_formula" %in% names(tempcall))) stop("No covar_formula model specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("e_treat_1" %in% names(tempcall))) stop("No e_treat_1 specified")
  if (!("e_treat_2" %in% names(tempcall))) stop("No e_treat_2 specified")
  if (!("e_treat_3" %in% names(tempcall))) stop("No e_treat_3 specified")
  if (!("e_treat_4" %in% names(tempcall))) stop("No e_treat_4 specified")
  if (!("degree" %in% names(tempcall)))  {if(!(tempcall$degree %in% c(1, 1))) stop("degree must be 1 or 2")}




  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]

  )


  # create the formula used to create covar_mat.  covar_mat will contain the variables used in the basis functions.
  formula_covar = eval(parse(text = paste(deparse(tempcall$treat, width.cutoff = 500), deparse(tempcall$covar_formula, width.cutoff = 500), sep = "")))



  # utils::str(m_frame <- model.frame(formula_covar, data))
  m_frame <- model.frame(formula_covar, data)
  covar_matrix_temp <- model.matrix(formula_covar, m_frame)
  covar_mat <- covar_matrix_temp[, -1]   # remove the first column




  if (degree == 2){

  # three basis functions: constant, T, and T^2
  B <- cbind(1,
             tempdat$treat,
             tempdat$treat^2)

  n_B <- ncol(B)

  Xstar <- cbind(1,
                 covar_mat)

  n_Xstar <- ncol(Xstar)


  # E_mat is a d x d matrix containing conditional expectations of T given Xs
  E_mat <- matrix(numeric(n_B * n_B),
                  nrow = n_B)

  wBB <- matrix(numeric(n_B * n_Xstar * n_B * n_Xstar),
                nrow = (n_B * n_Xstar) )


  for (  i in 1:nrow(covar_mat)  ) {
    E_mat[1, 1] <- 1
    E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
    E_mat[1, 3] <- E_mat[2, 2] <- E_mat[3, 1] <- e_treat_2[i]
    E_mat[2, 3] <- E_mat[3, 2] <- e_treat_3[i]
    E_mat[3, 3] <- e_treat_4[i]

#     # try these three lines...
#     W_Z_vec <- ( solve(E_mat) %*% B[i,] ) %x% t( as.matrix(Xstar[i,]) )
#     Z_vec <- B[i,]  %x% t( as.matrix(Xstar[i,]) )
#     wBB <- wBB + W_Z_vec %x% t(Z_vec)

    # alternative way to get wBB
    wBB <- wBB + (solve(E_mat) %*% B[i,] %*% t(B[i,])) %x%
      (as.matrix(Xstar[i,]) %*% t(as.matrix(Xstar[i,])))
  }
  # add some noise to prevent singularity?
  # wBB <- wBB + matrix(runif(prod(dim(wBB)), min = -0.00001, max = 0.0001), nrow = dim(wBB)[1])
  # should be solve, but there is a singularity or wBB is not invertible
  # xi_1st_factor <- ginv(wBB)

  xi_1st_factor <- solve(wBB[, ])

  #     temp_xi_1st <- solve(wBB)
  #     summary(temp_xi_1st)

  E_mat <- matrix(numeric(n_B * n_B),
                  nrow = n_B)

  wBB_2 <- matrix(numeric(n_B * n_Xstar ),
                  nrow = (n_B * n_Xstar ) )

  for (i in 1:nrow(covar_mat)) {
    E_mat[1, 1] <- 1
    E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
    E_mat[1, 3] <- E_mat[2, 2] <- E_mat[3, 1] <- e_treat_2[i]
    E_mat[2, 3] <- E_mat[3, 2] <- e_treat_3[i]
    E_mat[3, 3] <- e_treat_4[i]

    wBB_2 <- wBB_2 +
      ( ( solve(E_mat) %*% as.matrix(B[i,]) )  %x%  as.matrix(Xstar[i,]) ) %*% tempdat$Y[i]
  }

  xi_2nd_factor <- wBB_2



  vec_gam_star <- xi_1st_factor %*% xi_2nd_factor

  xi_alpha_hat <- vec_gam_star[1:n_Xstar] %*%  colMeans(Xstar[,])
  xi_beta_hat <- vec_gam_star[(n_Xstar + 1):(2 * n_Xstar) ] %*%  colMeans(Xstar[,])
  xi_gamma_hat <- vec_gam_star[(2 * n_Xstar + 1):(3 * n_Xstar)] %*%  colMeans(Xstar[,])

  wtrg_coefs <- c(xi_alpha_hat, xi_beta_hat, xi_gamma_hat)

  } else if (degree == 1 ){

    B <- cbind(1, tempdat$treat)
    n_B <- ncol(B)
    Xstar <- cbind(1, covar_mat)
    n_Xstar <- ncol(Xstar)

    # E_mat is a d x d matrix containing conditional expectations of T given Xs
    E_mat <- matrix(numeric(n_B * n_B), nrow = n_B)

    wBB <- matrix(numeric(n_B * n_Xstar * n_B * n_Xstar), nrow = (n_B * n_Xstar) )

    for (  i in 1:nrow(covar_mat)  ) {
      E_mat[1, 1] <- 1
      E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
      E_mat[2, 2] <- e_treat_2[i]

#       # try these three lines...
#       W_Z_vec <- ( solve(E_mat) %*% B[i,] ) %x% t( as.matrix(Xstar[i,]) )
#       Z_vec <- B[i,]  %x% t( as.matrix(Xstar[i,]) )
#       wBB <- wBB + W_Z_vec %x% t(Z_vec)

          wBB <- wBB + (solve(E_mat) %*% B[i,] %*% t(B[i,])) %x%
         (as.matrix(Xstar[i,]) %*% t(as.matrix(Xstar[i,])))
    }

    xi_1st_factor <- solve(wBB)
    E_mat <- matrix(numeric(n_B * n_B), nrow = n_B)
    wBB_2 <- matrix(numeric(n_B * n_Xstar ), nrow = (n_B * n_Xstar ) )


    for (i in 1:nrow(covar_mat)) {
      E_mat[1, 1] <- 1
      E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
      E_mat[2, 2] <- e_treat_2[i]

#       W_Z_vec <- c(  solve(E_mat) %*% B[i,]  %o%  as.matrix(Xstar[i,]) )   # tensor product and then turned into a vector
#       wBB_2 <- wBB_2 + W_Z_vec * tempdat$Y[i]

      wBB_2 <- wBB_2 +  ( ( solve(E_mat) %*% as.matrix(B[i,]) )  %x%  as.matrix(Xstar[i,]) ) %*% tempdat$Y[i]
    }


    xi_2nd_factor <- wBB_2

    vec_gam_star <- xi_1st_factor %*% xi_2nd_factor

    xi_alpha_hat <- vec_gam_star[1:n_Xstar] %*%  colMeans(Xstar[,])
    xi_beta_hat <- vec_gam_star[(n_Xstar + 1):(2 * n_Xstar) ] %*%  colMeans(Xstar[,])

    wtrg_coefs <- c(xi_alpha_hat, xi_beta_hat)

  }else{
    stop("Error: degree needs to be 1 or 2")
  }



  z_object <- list(param = wtrg_coefs,
                   call = tempcall)

  class(z_object) <- "causaldrf_simple"
  class(z_object) <- append(class(z_object), "causaldrf")
  z_object

}




##' @export
summary.causaldrf_simple <- function(object, ...){
  cat("\nEstimated outcomes at grid values:\n")
  print(object$param)
#   cat("\nOutcome Summary:\n")
#   print(ls.print(object$out_mod))

}


##' @export
print.summary.causaldrf_simple <- function(x, ...){
  #   cat("\nTreatment Summary:\n")
  #   print(summary(x$t_mod))
  cat("\nEstimated values:\n")
  print(x$param)
}


##' @export
print.causaldrf_simple <- function(x,...){
  #   cat("Call:\n")
  #   print(x$call)
  cat("\nEstimated values:\n")
  print(x$param)
}


