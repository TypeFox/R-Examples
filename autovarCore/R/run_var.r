#' Calculate the VAR model and apply restrictions
#'
#' This function calls the \code{vars::var} function to calculate the VAR model and applies restrictions if needed. We set the intercept to 1 for restricted equations because calculations go wrong otherwise (this is a bug in the vars library).
#' @param endo_matrix A numeric matrix of endogenous data.
#' @param exo_matrix Either \code{NULL} or a numeric matrix of exogenous data.
#' @param lag A nonnegative integer specifying the lag length of the model. Specifying 0 for the lag results in calculating a lag 1 model with all lag-1 terms restricted.
#' @return A \code{varest} object with the VAR estimation result.
#' @examples
#' endo_matrix <- matrix(rnorm(120), ncol = 2, nrow = 60,
#'                       dimnames = list(NULL, c("rumination", "activity")))
#' autovarCore:::run_var(endo_matrix, NULL, 1)
#' @export
run_var <- function(endo_matrix, exo_matrix, lag) {
  varest <- VAR(y = endo_matrix, p = max(lag, 1), exogen = exo_matrix)
  resmat_ncol <- length(colnames(varest$datamat)[!(colnames(varest$datamat) %in% colnames(varest$y))])
  resmat_nrow <- ncol(varest$y)
  restrictions <- restrictions_for_lag(lag, resmat_nrow, resmat_ncol)
  if (!is.null(restrictions)) {
    varest <- restrict(varest, method = "manual", resmat = restrictions)
    for (i in 1:length(varest$varresult))
      attr(varest$varresult[[i]]$terms, "intercept") <- 1
  }
  varest
}

restrictions_for_lag <- function(lag, resmat_nrow, resmat_ncol) {
  restrictions <- NULL
  if (lag == 0) {
    # Assume a lag 1 model. Restrict all the lag 1 columns.
    resmat <- rep.int(1, resmat_ncol)
    resmat[1:resmat_nrow] <- 0
    resmat <- rep(resmat, resmat_nrow)
    restrictions <- matrix(resmat, nrow = resmat_nrow, byrow = TRUE)
  } else if (lag == 2) {
    # Restrict the second lag in all models except the one using it.
    resmat <- rep.int(1, resmat_ncol)
    resmat[(resmat_nrow + 1):(2 * resmat_nrow)] <- 0
    resmat <- rep(resmat, resmat_nrow)
    resmat[seq(resmat_nrow + 1, length(resmat), resmat_ncol + 1)] <- 1
    restrictions <- matrix(resmat, nrow = resmat_nrow, byrow = TRUE)
  }
  restrictions
}
