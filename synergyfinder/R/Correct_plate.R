#' Baseline correction for the dose-response matrix of drug combinations
#'
#' A function to do baseline correction on the dose-response matrix for drug combinations with a weighted correction fator
#'
#' @param response.mat a dose-response matrix with concentrations as row names and column names
#' @return A list of the original dose-response matrix without correction and the corrected dose-response matrix
#' @author Liye He \email{liye.he@helsinki.fi}, Jing Tang \email{jing.tang@helsinki.fi}
BaselineCorrectionSD <- function (response.mat) {
  pm <- response.mat
  # check if response.mat has row names and column names
  if (is.null(rownames(response.mat)) | is.null(colnames(response.mat))) {
    stop("Please provide drug contrations as row names and column names!")
  }

  single.fitted <- FittingSingleDrug(response.mat)

  baseline <- (min(as.numeric(single.fitted$drug.row.fitted)) +
                 min(as.numeric(single.fitted$drug.col.fitted)))/2

  # correct matrix by a weighted correction fator
  pm.corrected <- pm - ((100 - pm) / 100 * baseline)
  output <- list(original.mat = pm, corrected.mat = pm.corrected)
  return(output)
}
