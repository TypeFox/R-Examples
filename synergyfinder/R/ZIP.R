#' Delta synergy score based on zero interaction potency (ZIP) model
#'
#' A function to calculate delta synergy score based on zero interaction potency (ZIP) model
#'
#' @param response.mat a dose-response matrix with concentrations as row names and column names
#' @param correction a parameter to specify if baseline correction is used or not. Defaults to TRUE.
#' @return A matrix of delta scores for all the dose pairs for a drug combination. For a does pair with at least one zero concentration, 0 is used as the synergy score.
#' @author Liye He \email{liye.he@helsinki.fi}, Jing Tang \email{jing.tang@helsinki.fi}
#' @references Yadav B, Wennerberg K, Aittokallio T, Tang J. Searching for Drug Synergy in Complex Dose-Response Landscape Using an Interaction Potency Model.
#' Computational and Structural Biotechnology Journal 2015; 13: 504-513.
ZIP <- function(response.mat, correction = TRUE) {
  if(correction) {
    # correct the response data
    response.mat <- BaselineCorrectionSD(response.mat)$corrected.mat
  }
  # Fitting single drugs using logistic functions
  # NA values treated
  single.fitted <- FittingSingleDrug(response.mat, fixed = c(NA, 0, 100, NA))
  drug.col.response <- single.fitted$drug.col.fitted
  drug.row.response <- single.fitted$drug.row.fitted
  # Update the first row and first column
  updated.single.mat <- mat.or.vec(nrow(response.mat),ncol(response.mat))
  colnames(updated.single.mat) <- colnames(response.mat)
  rownames(updated.single.mat) <- rownames(response.mat)
  updated.single.mat[1, c(2:ncol(response.mat))] <- drug.col.response
  updated.single.mat[c(2:nrow(response.mat)), 1] <- drug.row.response

  # Update the column2-column8
  updated.col.mat <- updated.single.mat
  for (i in 2:ncol(response.mat)){
    tmp <- as.data.frame(mat.or.vec(nrow(response.mat) - 1, 0))
    tmp$dose <- as.numeric(rownames(response.mat)[-1])
    tmp$inhibition <- response.mat[c(2:nrow(response.mat)), i]
    tmp.min <- updated.single.mat[1, i]
    if (var(tmp$inhibition, na.rm = T) == 0) { ## no variance in the drug responses
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }
    tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, tmp.min, 100,NA)),
                     ,na.action = na.omit)
    tmp$fitted.inhibition <- suppressWarnings(fitted(tmp.model))
    if(tmp$fitted.inhibition[nrow(response.mat) - 1] < 0) tmp$fitted.inhibition[nrow(response.mat) - 1] <- tmp.min
    updated.col.mat[c(2:nrow(response.mat)), i] <- tmp$fitted.inhibition
  }

  # Update the row2-row8
  updated.row.mat <- updated.single.mat
  for (i in 2:nrow(response.mat)){
    tmp <- as.data.frame(mat.or.vec(ncol(response.mat) - 1, 0))
    tmp$dose <- as.numeric(colnames(response.mat)[-1])
    tmp$inhibition <- response.mat[i, c(2:ncol(response.mat))]
    tmp.min <- updated.single.mat[i, 1]
    if (var(tmp$inhibition, na.rm = T) == 0) { ## no variance in the drug responses
      tmp$inhibition[1] <- tmp$inhibition[1] - 10^-10
    }
    tmp.model <- drm(inhibition ~ dose, data = tmp, fct = L.4(fixed = c(NA, tmp.min, 100,NA)),
                                      na.action = na.omit)
    tmp$fitted.inhibition <- suppressWarnings(fitted(tmp.model))
    if(tmp$fitted.inhibition[ncol(response.mat) - 1] < 0) tmp$fitted.inhibition[ncol(response.mat) - 1] <- tmp.min
    updated.row.mat[i, c(2:ncol(response.mat))] <- tmp$fitted.inhibition
  }

  # take average of updated.col.mat and updated.row.mat as fitted.mat
  fitted.mat <- (updated.col.mat + updated.row.mat) / 2

  # make zip.mat based on updated.single.mat
  zip.mat <- updated.single.mat
  for (i in 2:nrow(updated.single.mat)){
    for (j in 2:ncol((updated.single.mat))){
      zip.mat[i,j] <- updated.single.mat[i, 1] + updated.single.mat[1, j] - updated.single.mat[i, 1] * updated.single.mat[1, j] / 100
    }
  }
  # negative and positive controls are removed
  fitted.mat[1, 1] <- 0
  zip.mat[1, 1] <- 0
  # cannot be over 100 for the estimation
  fitted.mat <- apply(fitted.mat, c(1, 2), function(x) ifelse(x > 100, 100, x))
  delta.mat <- (fitted.mat - zip.mat)
  delta.mat
}
