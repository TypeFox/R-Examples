##' The inverse second moment weighting (ismw) estimator
##'
##' This method estimates the ADRF by using weighting matrices instead of
##' scalars.  The weight matrices require conditional expectations of the
##' treatment and higher order conditional expectations. It uses outputs from
##' the \code{t_mod} function.
##'
##'
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param e_treat_1 a vector, representing the conditional expectation of
##' \code{treat} from \code{t_mod}.
##' @param e_treat_2 a vector, representing the conditional expectation of
##' \code{treat^2} from \code{t_mod}.
##' @param e_treat_3 a vector, representing the conditional expectation of
##' \code{treat^3} from \code{t_mod}.
##' @param e_treat_4 a vector, representing the conditional expectation of
##' \code{treat^4} from \code{t_mod}.
##' @param degree is 1 for linear and 2 for quadratic outcome model.
##'
##' @details
##' This function estimates the ADRF requires estimated moments and
##' uses the outputs of the t_mod function as inputs.  For more details,
##' see Schafer and Galagate (2015).
##'
##' @return \code{ismw_est} returns an object of class "causaldrf_simple",
##' a list that contains the following components:
##' \item{param}{the estimated parameters.}
##' \item{call}{the matched call.}
##'
##' @seealso  \code{\link{iptw_est}}, \code{\link{ismw_est}},
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
##' @examples
##'
##' ## Example from Schafer (2015).
##'
##' example_data <- sim_data
##'
##' t_mod_list <- t_mod(treat = T,
##'                 treat_formula = T ~ B.1 + B.2 + B.3 + B.4 + B.5 + B.6 + B.7 + B.8,
##'                 data = example_data,
##'                 treat_mod = "Normal")
##'
##' cond_exp_data <- t_mod_list$T_data
##'
##' full_data <- cbind(example_data, cond_exp_data)
##'
##' ismw_list <- ismw_est(Y = Y,
##'                       treat = T,
##'                       data = full_data,
##'                       e_treat_1 = full_data$est_treat,
##'                       e_treat_2 = full_data$est_treat_sq,
##'                       e_treat_3 = full_data$est_treat_cube,
##'                       e_treat_4 = full_data$est_treat_quartic,
##'                       degree = 1)
##'
##' sample_index <- sample(1:1000, 100)
##'
##' plot(example_data$T[sample_index],
##'       example_data$Y[sample_index],
##'       xlab = "T",
##'       ylab = "Y",
##'       main = "ismw estimate")
##'
##' abline(ismw_list$param[1],
##' ismw_list$param[2],
##'         lty=2,
##'         lwd = 2,
##'         col = "blue")
##'
##' legend('bottomright',
##'         "ismw estimate",
##'         lty=2,
##'         lwd = 2,
##'         col = "blue",
##'         bty='Y',
##'         cex=1)
##'
##' rm(example_data, t_mod_list, cond_exp_data, full_data, ismw_list, sample_index)
##'
##'
##' @usage
##' ismw_est(Y,
##'          treat,
##'          data,
##'          e_treat_1,
##'          e_treat_2,
##'          e_treat_3,
##'          e_treat_4,
##'          degree )
##'
##' @export
##'
##'
##'


ismw_est <- function(Y,
                     treat,
                     data,
                     e_treat_1,
                     e_treat_2,
                     e_treat_3,
                     e_treat_4,
                     degree ){
  # Y is the name of the Y variable
  # treat is the name of the treatment variable
  # data will contain all the data: X, treat, and Y
  # e_treat_1 is the vector containing the estimated means
  # e_treat_2 is the vector containing the estimated second moments
  # e_treat_3 is the vector containing the estimated third moments
  # e_treat_4 is the vector containing the estimated fourth moments
  # degree is 1 for linear and 2 for quadratic outcome model.

  # The outcome is a vector containing the estimated parameters.



  #save input
  tempcall <- match.call()

  #some basic input checks
  if (!("Y" %in% names(tempcall))) stop("No Y variable specified")
  if (!("treat" %in% names(tempcall)))  stop("No treat variable specified")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!("e_treat_1" %in% names(tempcall))) stop("No e_treat_1 specified")
  if (!("e_treat_2" %in% names(tempcall))) stop("No e_treat_2 specified")
  if (!("e_treat_3" %in% names(tempcall))) stop("No e_treat_3 specified")
  if (!("e_treat_4" %in% names(tempcall))) stop("No e_treat_4 specified")
  if (!("degree" %in% names(tempcall)))  {if(!(tempcall$degree %in% c(1, 2))) stop("degree must be 1 or 2")}



  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    Y = data[,as.character(tempcall$Y)],
    treat = data[,as.character(tempcall$treat)]

  )

  if (degree == 2){

  B <- cbind(1,
             tempdat$treat,
             tempdat$treat^2 )

  n_B <- ncol(B)

  E_mat <- matrix(numeric(n_B * n_B),
                  nrow = n_B)

  # wBB will contain the first sum or first factor of the ismw estimator
  wBB <- matrix(numeric(n_B * n_B), nrow = n_B)

  for (i in 1:nrow(data) ) {
    E_mat[1, 1] <- 1
    E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
    E_mat[2, 2] <- E_mat[1, 3] <- E_mat[3, 1] <-  e_treat_2[i]
    E_mat[2, 3] <- E_mat[3, 2] <- e_treat_3[i]
    E_mat[3, 3] <- e_treat_4[i]
    # just use the data variance of the simulated dataset.
    wBB <- wBB +  solve(E_mat) %*% B[i,] %*% t(B[i,])
  }

  xi_1st_factor <- solve(wBB)


  E_mat <- matrix(numeric(n_B * n_B), nrow = n_B)
  wBB_2 <- matrix(numeric(n_B), nrow = n_B)

  for (i in 1:nrow(data) ) {
    E_mat[1, 1] <- 1
    E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
    E_mat[2, 2] <- E_mat[1, 3] <- E_mat[3, 1] <-  e_treat_2[i]
    E_mat[2, 3] <- E_mat[3, 2] <- e_treat_3[i]
    E_mat[3, 3] <- e_treat_4[i]
    ## just use the data variance of the simulated dataset.
    wBB_2 <- wBB_2 + solve(E_mat) %*% as.matrix(B[i,]) %*% tempdat$Y[i]
  }

  xi_2nd_factor <- wBB_2

  ismw_coefs <- xi_1st_factor %*% xi_2nd_factor

  } else if (degree == 1 ){

    B <- cbind(1, tempdat$treat )
    n_B <- ncol(B)

    # E_mat is a d x d matrix containing conditional expectations of T given Xs
    E_mat <- matrix(numeric(n_B * n_B), nrow = n_B)
    # wBB will contain the first sum or first factor of the ismw estimator
    wBB <- matrix(numeric(n_B * n_B), nrow = n_B)

    for (i in 1:nrow(data) ) {
      E_mat[1, 1] <- 1
      E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
      E_mat[2, 2] <- e_treat_2[i]
      wBB <- wBB +  solve(E_mat) %*% B[i,] %*% t(B[i,])
    }

    xi_1st_factor <- solve(wBB)

    E_mat <- matrix(numeric(n_B * n_B), nrow = n_B)
    wBB_2 <- matrix(numeric(n_B), nrow = n_B)

    for (i in 1:nrow(data) ) {
      E_mat[1, 1] <- 1
      E_mat[1, 2] <- E_mat[2, 1] <- e_treat_1[i]
      E_mat[2, 2] <- e_treat_2[i]

      wBB_2 <- wBB_2 + solve(E_mat) %*% as.matrix(B[i,]) %*% tempdat$Y[i]
    }

    xi_2nd_factor <- wBB_2

    ismw_coefs <- xi_1st_factor %*% xi_2nd_factor

  }else{
    stop("Error: degree needs to be 1 or 2")
  }


  z_object <- list(param = ismw_coefs,
                   call = tempcall)

  class(z_object) <- "causaldrf_simple"
  z_object
}




