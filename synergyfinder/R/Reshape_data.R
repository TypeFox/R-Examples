#' Transform the response data from the data frame format to dose-response matrixes
#'
#' A function to transform the response data from the data frame format to dose-response matrixes
#'
#' @param data drug combination response data in a data frame format
#' @param data.type a parameter to specify the response data type which can be either "viability" or "inhibition".
#' @return a list of the following components:
#' \item{dose.response.mats}{a list of the dose-response matrixes with \%inhibition as the response data. Row names and column names are drug concentrations.}
#' \item{drug.pairs}{a data frame contains the name of the row drug, the name of the column drug, concentration unit and block IDs.}
#' @details The input data must contain the following columns: BlockID, DrugRow, DrugCol, Row, Col, Response,
#' ConcRow, ConcCol, ConcUnit
#' @author Liye He \email{liye.he@helsinki.fi}
#' @examples 
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
ReshapeData <- function(data, data.type = "viability") {
  # check column names
  if (!all(c("BlockID", "DrugRow", "DrugCol", "Row", "Col", "Response", "ConcRow", "ConcCol",
             "ConcUnit") %in% colnames(data)))
    stop("The input data must contain the following columns: BlockID, DrugRow, DrugCol, Row, Col, Response,
         ConcRow, ConcCol, ConcUnit")
  # obtain BlockIDs
  id.drug.comb <- unique(data$BlockID)
  dose.response.mats <- list() ## store all the dose-response matrices
  drug.pairs <- data.frame(drug.row = character(length(id.drug.comb)),
                           drug.col = character(length(id.drug.comb)),
                           conc.unit = character(length(id.drug.comb)),
                           blockIDs = numeric(length(id.drug.comb)),
                           stringsAsFactors = F)
  for (i in 1:length(id.drug.comb)) {
    tmp.mat <- data[which(data$BlockID == id.drug.comb[i]), ]

    if (data.type == "viability") {
      tmp.mat$Inhibition <- 100 - tmp.mat$Response
    } else {
      tmp.mat$Inhibition <- tmp.mat$Response
    }



    # get single drug concentrations
    # first row of dose-response matrix: column concentrations
    conc.col <- tmp.mat$ConcCol[which(tmp.mat$Row == 1)]
    conc.col <- conc.col[order(tmp.mat$Col[which(tmp.mat$Row == 1)])]
    # first column of dose-response matrix: row concentrations
    conc.row <- tmp.mat$ConcRow[which(tmp.mat$Col == 1)]
    conc.row <- conc.row[order(tmp.mat$Row[which(tmp.mat$Col == 1)])]

    # response matrix for one drug combination
    response.mat <- acast(tmp.mat, Row ~ Col,value.var = "Inhibition")
    colnames(response.mat) <- conc.col
    rownames(response.mat) <- conc.row

    # adjust the dose-response matrix based the first concentration
    if (which.max(conc.row) == 1 & which.max(conc.col) == 1) {
      response.mat <- t(apply(apply(response.mat, 2, rev), 1, rev))
    } else if (which.max(conc.row) == length(conc.row) & which.max(conc.col) == 1) {
      response.mat <- t(apply(response.mat, 1, rev))
    } else if (which.max(conc.row) == 1 & which.max(conc.col) == length(conc.col)) {
      response.mat <- apply(response.mat, 2, rev)
    }
    conc.unit <- unique(tmp.mat$ConcUnit) ## concentration unit
    drug.row <- unique(tmp.mat$DrugRow)
    drug.col <- unique(tmp.mat$DrugCol)
    drug.pairs$drug.row[i] <- drug.row
    drug.pairs$drug.col[i] <- drug.col
    drug.pairs$concUnit[i] <- conc.unit
    # save dose-response matrix
    dose.response.mats[[i]] <- response.mat
  }
  drug.pairs$blockIDs <- id.drug.comb

  return (list(dose.response.mats = dose.response.mats, drug.pairs = drug.pairs))
}
