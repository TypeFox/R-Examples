#' HSA synergy score based on highest single agent (HSA) model
#'
#' A function to calculate HSA synergy score based on HSA model
#'
#' @param response.mat a dose-response matrix with concentrations as row names and column names
#' @param correction a parameter to specify if baseline correction is used or not. Defaults to TRUE.
#' @return A matrix of HSA synergy scores for all the dose pairs for a drug combination. For a
#' does pair with at least one zero concentration, 0 is used as the synergy score.
#' @author Liye He \email{liye.he@helsinki.fi}
HSA <- function(response.mat, correction = TRUE) {
  if(correction) {
    # correct the response data
    response.mat <- BaselineCorrectionSD(response.mat)$corrected.mat
  }
  drug1.response <- response.mat[, 1]
  drug2.response <- response.mat[1, ]
  # reference matrix
  ref.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      ref.mat[i, j] <- ifelse(drug1.response[i] > drug2.response[j],
                              drug1.response[i], drug2.response[j])
    }
  }
  # synergy matrix
  syn.mat <- response.mat - ref.mat
  syn.mat
}
