#' Plate Reader Data of Bacterial Growth
#'
#' Example data set from growth experiments with different
#' concentrations of antibiotics.
#'
#' This rather 'difficult' data set was intentionally selected to make model
#' fitting by the package more challenging.
#'
#' @format Data frame with the following columns:
#' \describe{
#'   \item{strain}{identifier of the bacterial strain, D=donor, R=recipient, T=transconjugant.}
#'   \item{replicate}{replicate of the trial.}
#'   \item{conc}{concentration of the antibiotics (Tetracycline).}
#'   \item{time}{time in hours.}
#'   \item{value}{bacteria concentration measured as optical density.}
#' }
#'
#' @source Claudia Seiler, TU Dresden, Institute of Hydrobiology.
#'
#' @name bactgrowth
#' @docType data
#' @keywords data
#'
#' @examples
#' ## plot data and determine growth rates
#' data(bactgrowth)
#'
#'
#' library(lattice)
#' xyplot(value ~ time | strain + as.factor(conc),
#'       data = bactgrowth, groups = replicate)

NULL
