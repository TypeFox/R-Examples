#' \code{FWDselect}: Selecting Variables in Regression Models.
#'
#' This package introduces a simple method to select the best model
#' using different types of data (binary, Gaussian or Poisson) and applying it
#' in different contexts (parametric or non-parametric). The proposed method is a
#' new forward stepwise-based selection procedure that  selects a model
#' containing a subset of variables according to an optimal criterion (obtained
#' by cross-validation) and also takes into account the computational cost.
#' Additionally,  bootstrap resampling techniques are used to implement tests
#' capable of detecting whether significant effects of the unselected variables
#' are present in the model.
#'
#' @name FWDselect
#' @docType package
#' @details \tabular{ll}{ Package: \tab FWDselect\cr Type: \tab Package\cr
#' Version: \tab 2.1.0\cr Date: \tab 2015-12-18\cr License: \tab MIT + file LICENSE\cr}
#'
#' \code{FWDselect} is just a shortcut for ``Forward selection'' and is a very
#' good summary of one of the package's major functionalities, i.e., that of
#' providing a forward stepwise-based selection procedure. This software helps
#' the user select relevant variables and evaluate how many of these need to be
#' included in a regression model. In addition, it enables both numerical and
#' graphical outputs to be displayed. The package includes several functions
#' that enable users to select the variables to be included in linear,
#' generalized linear or generalized additive regression models. Users can
#' obtain the best combinations of \code{q} variables by means of the main
#' function \code{\link{selection}}. Additionally, if one wants to obtain the
#' results for more than one size of subset, it is possible to apply the
#' \code{\link{qselection}} function, which returns a summary table showing the
#' different subsets, selected variables and information criterion values. The
#' object obtained when using this last function is the argument required for
#' \code{\link{plot.qselection}}, which provides a graphical output. Finally, to
#' determine the number of variables that should be introduced into the model,
#' only the \code{\link{test}} function needs to be applied.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @references Burnham, K., Anderson, D. (2002). Model selection and multimodel
#'   inference: a practical information-theoretic approach. 2nd Edition
#'   Springer.
#'
#'   Efron, B. (1979). Bootstrap methods: another look at the jackknife. Annals
#'   of Statistics, 7:1-26.
#'
#'   Efron, B. and Tibshirani, R. J. (1993). An introduction to the Bootstrap.
#'   Chapman and Hall, London.
#'
#'   Miller, A. (2002). Subset selection in regression. Champman and Hall.
#'
#'   Sestelo, M., Villanueva, N. M. and Roca-Pardinas, J. (2013). FWDselect:
#'   Variable selection algorithm in regression models. Discussion Papers in
#'   Statistics and Operation Research, University of Vigo, 13/02.
#'
NULL
