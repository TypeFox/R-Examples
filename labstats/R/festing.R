#' @title Mouse liver enzyme randomised block design data set
#'
#' @description Sixteen mice from four strains were used; half were
#' assigned to the treated condition and half to the control. The
#' experiment was split into two batches, which were run two months
#' apart.
#'
#' @details The purpose of the experiment was to test if diallyl
#' sulphide (DS) affects the activity of the liver enzyme Gst. Four
#' strains of mice were used and the experiment was conducted in two
#' batches, where the housing conditions differed between
#' batches. Both batch and strain are considered blocking variables.
#'
#' @format A data frame with 16 rows and 4 variables:
#' \describe{
#' 
#'   \item{strain:}{Four strains of mice were used: NIH, BALB/c, A/J,
#' and 129/Ola.}
#' 
#'   \item{treatment:}{Presence or absence of diallyl sulphide.}
#' 
#'   \item{batch:}{The experiment was conducted in two batches.}
#' 
#'   \item{value:}{Levels of glutathione-S-transferase (Gst).}  }
#'
#' @references Festing MFW (2014). Randomized block experimental
#' designs can increase the power and reproducibility of laboratory
#' animal experiments. \emph{ILAR Journal} 55(3):472-476.
#' 
#' @docType data
#' @name festing
NULL
