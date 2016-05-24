##################################################################################
#' Spot Meta Flatten Fbs Row
#' 
#' Helper function that replaces vector-valued columns in a final best solution
#' file row with N "flat" values, where N is the dimension of the vector replaced.
#' Also, each column will be converted to a character string.
#'
#' @param fbsRow The FBS row to flatten.
#' 
#' @return list \cr
#' - The FBS line, with vector-valued columns replaced with scalar columns.
#' @export
#' @keywords internal
###################################################################################
spotMetaFlattenFbsRow<-function(fbsRow) {
  replaceVectorsInFbsColumn <- function(fbsColumn, fbsColumnName) {
    if (is.list(fbsColumn)) {
      vectorValuedFbsColumn <- fbsColumn[[1]]
      vectorValuedFbsColumnAsList <- as.list(vectorValuedFbsColumn)
      names(vectorValuedFbsColumnAsList) <- paste(fbsColumnName,
                                                  1:length(vectorValuedFbsColumnAsList),
                                                  sep = "")
      return(vectorValuedFbsColumnAsList)
    } else {
      fbsColumnAsList <- list(fbsColumn)
      names(fbsColumnAsList) <- list(fbsColumnName)
      return(fbsColumnAsList)
    }
  }
  fbsColumnWithoutVectors <- Map(replaceVectorsInFbsColumn, fbsRow, names(fbsRow))
  flatFbsColumn <- Reduce(c, init = list(), fbsColumnWithoutVectors)
  Map(as.character, flatFbsColumn)
}

