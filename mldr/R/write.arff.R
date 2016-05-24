#' @title Writes an \code{mldr} object to a file
#' @description Save the \code{mldr} content to an ARFF file and the label data to an XML file
#' @param obj The \code{mldr} object whose content is going to be written
#' @param filename Base name for the files (without extension)
#' @param write.xml \code{TRUE} or \code{FALSE}, stating if the XML file has to be written
#' @seealso \code{\link{mldr_from_dataframe}}, \code{\link{mldr}}
#' @examples
#'
#' library(mldr)
#'
#' write_arff(emotions, "myemotions")
#' @export
write_arff <- function(obj, filename, write.xml = FALSE) {
  # Open file
  connection <- file(paste(filename, ".arff", sep = ""))

  # Create header an attribute lines
  header <- paste("@relation ", obj$name, sep = "")
  attributes <- paste("@attribute ", names(obj$attributes), " ", obj$attributes, sep = "")
  data <- obj$dataset[, 1:obj$measures$num.attributes]
  data[is.na(data)] <- '0' # NAs aren't missing values ('?') but 0
  data <- apply(data, 1, function(c) paste(c, collapse = ','))

  # Write header, attributes and data
  writeLines(c(header, "", attributes, "", "@data", data), connection)
  close(connection)

  if(write.xml) {
    # Open XML file
    connection <- file(paste(filename, ".xml", sep = ""))

    xmlheader <- '<?xml version="1.0" encoding="utf-8"?>'
    labelstag <- '<labels xmlns="http://mulan.sourceforge.net/labels">'
    labeltags <- paste(c('<label name="'), rownames(obj$labels), c('"></label>'), sep = "")
    labelsend <- '</labels>'

    writeLines(c(xmlheader, labelstag, labeltags, labelsend), connection)
    close(connection)
  }
}
