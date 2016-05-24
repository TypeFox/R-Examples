#' Delta synergy score based on Loewe model
#'
#' A function to calculate Loewe synergy score based on Loewe model
#'
#' @param response.mat a dose-response matrix with concentrations as row names and column names
#' @param correction a parameter to specify if baseline correction is used or not. Defaults to TRUE.
#' @return A matrix of Loewe synergy scores for all the dose pairs for a drug combination. For a does pair with at least one zero concentration, 0 is used as the synergy score.
#' @author Liye He \email{liye.he@helsinki.fi}
#' @references Yadav B, Wennerberg K, Aittokallio T, Tang J. Searching for Drug Synergy in Complex Dose-Response Landscape Using an Interaction Potency Model.
#' Computational and Structural Biotechnology Journal 2015; 13: 504-513.
Loewe <- function(response.mat, correction = TRUE) {
  if(correction) {
    # correct the response data
    response.mat <- BaselineCorrectionSD(response.mat)$corrected.mat
  }

  single.fit <- FittingSingleDrug(response.mat)
  # column drug
  drug.col.model <- single.fit$drug.col.model
  drug.col.par <- coef(drug.col.model)
  d1.fun <- function(conc, drug.col.par) {
    (drug.col.par[3] + drug.col.par[2] *
       (conc / drug.col.par[4]) ^ drug.col.par[1]) /
      (1 + (conc / drug.col.par[4]) ^ drug.col.par[1])
  }

  # row drug
  drug.row.model <- single.fit$drug.row.model
  drug.row.par <- coef(drug.row.model)
  d2.fun <- function(conc, drug.row.par) {
    (drug.row.par[3] + drug.row.par[2] * (conc / drug.row.par[4]) ^ drug.row.par[1]) / (1 + (conc / drug.row.par[4]) ^ drug.row.par[1])
  }
  row.conc <- as.numeric(rownames(response.mat))[-1]
  col.conc <- as.numeric(colnames(response.mat))[-1]
  # construct the equation
  loewe.mat <- response.mat
  for (i in 1:length(col.conc)) {
    for (j in 1:length(row.conc)) {
      x1 <- col.conc[i]
      x2 <- row.conc[j]
      eq <- function (x) {
        x1 / (drug.col.par[4] * (((x - drug.col.par[3]) / (drug.col.par[2] - x)) ^ (1/drug.col.par[1]))) +
          x2 / (drug.row.par[4] * (((x - drug.row.par[3]) / (drug.row.par[2] - x)) ^ (1/drug.row.par[1]))) - 1
      }
      slv <- nleqslv(max(drug.col.par[2] + 1, drug.row.par[2] + 1), eq, method = "Newton")
      if (slv$termcd == 1) {
        loewe.mat[i + 1, j + 1] <- slv$x
      } else {
        y.loewe1 <- d1.fun(x1 + x2, drug.col.par)
        y.loewe2 <- d2.fun(x1 + x2, drug.row.par)
        loewe.mat[i + 1, j + 1] <- ifelse(y.loewe1 > y.loewe2, y.loewe1, y.loewe2)
      }

    }
  }
  return(response.mat - loewe.mat)
}
