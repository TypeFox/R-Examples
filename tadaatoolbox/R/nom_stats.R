#' Simple Chi^2
#'
#' This is a very simple wrapper for \link{chisq.test}.
#' @param x Dependent variable. Alternatively a \code{table}.
#' @param y Independent variable
#' @param correct Apply correction, passed to \code{chisq.test}.
#' @return A \code{numeric} value
#' @export
#' @import stats
#' @examples
#' nom_chisqu(ngo$abschalt, ngo$geschl)
nom_chisqu <- function(x, y = NULL, correct = FALSE){
  if (is.table(x)) {
    as.numeric(chisq.test(x = x, correct = correct)$statistic)
  } else {
    as.numeric(chisq.test(x = x, y = y, correct = correct)$statistic)
  }
}

#' Phi coefficient
#'
#' Very simple wrapper for \link[vcd]{assocstats}.
#' @param x Dependent variable. Alternatively a \code{table}.
#' @param y Independent variable
#'
#' @return \code{numeric} value
#' @export
#' @importFrom vcd assocstats
#' @examples
#' nom_phi(ngo$abschalt, ngo$geschl)
nom_phi <- function(x, y = NULL){
  if (!is.table(x)) {
    x <- table(x, y)
  }
  vcd::assocstats(x)$phi
}

#' Cramer's V
#'
#' Very simple wrapper for \link[vcd]{assocstats}.
#' @param x Dependent variable. Alternatively a \code{table}.
#' @param y Independent variable
#'
#' @return \code{numeric} value
#' @export
#' @importFrom vcd assocstats
#' @examples
#' nom_v(ngo$abschalt, ngo$geschl)
nom_v <- function(x, y = NULL){
  if (!is.table(x)) {
    x <- table(x, y)
  }
  vcd::assocstats(x)$cramer
}

#' Contingency Coefficient C
#'
#' Very simple wrapper for \link[vcd]{assocstats}.
#' @param x Dependent variable. Alternatively a \code{table}.
#' @param y Independent variable
#'
#' @return \code{numeric} value
#' @export
#' @importFrom vcd assocstats
#' @examples
#' nom_c(ngo$abschalt, ngo$geschl)
nom_c <- function(x, y = NULL){
  if (!is.table(x)) {
    x <- table(x, y)
  }
  vcd::assocstats(x)$contingency
}

#' Lambda
#'
#' Very simple wrapper for \link[ryouready]{nom.lambda}.
#' @param x Dependent variable. Alternatively a \code{table}.
#' @param y Independent variable
#' @param symmetric If \code{TRUE}, symmetric lambda is returned. Default is \code{FALSE}.
#' @param reverse If \code{TRUE}, row and column variable are switched.
#'
#' @return \code{numeric} value
#' @export
#' @importFrom ryouready nom.lambda
#' @examples
#' nom_lambda(ngo$abschalt, ngo$geschl)
nom_lambda <- function(x, y = NULL, symmetric = FALSE, reverse = FALSE){
  if (!is.table(x)) {
    x <- table(x, y)
  }
  if (symmetric) {
    ryouready::nom.lambda(x)$lambda.symmetric
  } else if (!reverse) {
    ryouready::nom.lambda(x)$lambda.rc
  } else {
    ryouready::nom.lambda(x)$lambda.cr
  }
}

#' Get all the nominal stats
#'
#' @param x Dependent variable. Alternatively a \code{table}.
#' @param y Independent variable
#' @param round Ho many digits should be rounded. Default is 2.
#' @param print Print method. Passed to \link[pixiedust]{sprinkle_print_method} as of now.
#' @return A \code{dust} object, depending on \code{print}.
#' @export
#' @import pixiedust
#' @family Tadaa-functions
#' @examples
#' tadaa_nom(ngo$abschalt, ngo$geschl)
tadaa_nom <- function(x, y = NULL, round = 2, print = "console"){
  if (!is.table(x)) {
    x <- table(x, y)
  }
  chisq  <- round(nom_chisqu(x), round)
  v      <- round(nom_v(x), round)
  cc     <- round(nom_c(x), round)
  lmbd_x <- round(nom_lambda(x), round)
  lmbd_y <- round(nom_lambda(x, reverse = T), round)
  lmbd_s <- round(nom_lambda(x, symmetric = T), round)

  ret <- data.frame("chisq" = chisq, "cv" = v, "c" = cc,
                    "lmbd_x" = lmbd_x, "lmbd_y" = lmbd_y,
                    "lmbd_s" = lmbd_s)

  retprint <- pixiedust::sprinkle_colnames(pixiedust::dust(ret), chisq = "Chi^2",
                                      cv = "Cramer's V",
                                      lmbd_x = "Lambda (x dep.)",
                                      lmbd_y = "Lambda (y dep.)",
                                      lmbd_s = "Lambda (sym.)")
  return(pixiedust::sprinkle_print_method(retprint, print))
}
