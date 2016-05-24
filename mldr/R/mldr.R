#' @title Creates an object representing a multilabel dataset
#' @description Reads a multilabel dataset from a file and returns an \code{mldr} object
#' containing the data and additional measures. The file has to be in ARFF format.
#' The label information could be in a separate XML file (MULAN style) or in the
#' the arff header (MEKA style)
#' @param filename Name of the dataset
#' @param use_xml Specifies whether to use an
#'  associated XML file to identify the labels. Defaults to TRUE
#' @param auto_extension Specifies whether to add
#'  the '.arff' and '.xml' extensions to the filename
#'  where appropriate. Defaults to TRUE
#' @param xml_file Path to the XML file. If not
#'  provided, the filename ending in ".xml" will be
#'  assumed
#' @param label_indices Optional vector containing the indices of the attributes
#'  that should be read as labels
#' @param label_names Optional vector containing the names of the attributes
#'  that should be read as labels
#' @param label_amount Optional parameter indicating the number of labels in the
#'  dataset, which will be taken from the last attributes of the dataset
#' @param force_read_from_file Set this parameter to TRUE to always read from a local file,
#'  or set it to FALSE to look for the dataset within the `mldr.datasets` package
#' @return An mldr object containing the multilabel dataset
#' @seealso \code{\link{mldr_from_dataframe}}, \code{\link{read.arff}}, \code{\link{summary.mldr}}
#' @examples
#'
#' library(mldr)
#'\dontrun{
#' # Read "yeast.arff" and labels from "yeast.xml"
#' mymld <- mldr("yeast")
#'
#' # Read "yeast-tra.arff" and labels from "yeast.xml"
#' mymld <- mldr("yeast-tra", xml_file = "yeast.xml")
#'
#' # Read "yeast.arff" specifying the amount of attributes to be used as labels
#' mymld <- mldr("yeast", label_amount = 14)
#'
#' # Read MEKA style dataset, without XML file and giving extension
#' mymld <- mldr("IMDB.arff", use_xml = FALSE, auto_extension = FALSE)
#'}
#' @import methods
#' @export
mldr <- function(filename,
                 use_xml = TRUE,
                 auto_extension = TRUE,
                 xml_file,
                 label_indices,
                 label_names,
                 label_amount,
                 force_read_from_file = !all(c(missing(xml_file), missing(label_indices), missing(label_names), missing(label_amount), use_xml, auto_extension))) {

  no_filename <- missing(filename)
  no_xml_file <- missing(xml_file)
  no_label_indices <- missing(label_indices)
  no_label_names <- missing(label_names)
  no_label_amount <- missing(label_amount)

  success <- FALSE

  if (!force_read_from_file) {
    if (requireNamespace("mldr.datasets", quietly = TRUE)) {
      if (exists(filename, mode = "list")) {
        ret_value <- get(filename)
        success <- TRUE
      } else if (existsFunction(filename)) {
        get(filename)()
        ret_value <- get(filename)
        success <- TRUE
      }
    }
  }

  if (success) {
    ret_value
  } else {
    if (!no_filename) {
      do.call(mldr_from_dataframe, read.arff(filename, use_xml, auto_extension, xml_file, label_indices, label_names, label_amount))
    } else {
      NULL
    }
  }
}

#' @title Generates an mldr object from a data.frame and a vector with label indices
#' @description This function creates a new \code{mldr} object from the data
#' stored in a \code{data.frame}, taking as labels the columns pointed by the
#' indexes given in a vector.
#' @param dataframe The \code{data.frame} containing the dataset attributes and labels.
#' @param labelIndices Vector containing the indices of attributes acting as labels. Usually the
#' labels will be at the end (right-most columns) or the beginning (left-most columns) of the \code{data.frame}
#' @param attributes Vector with the attributes type, as returned by the  \code{attributes} member of an \code{mldr}
#' object. By default the type of the data.frame columns will be used.
#' @param name Name of the dataset. The name of the dataset given as first parameter will be used by default
#' @return An mldr object containing the multilabel dataset
#' @seealso \code{\link{mldr}}, \code{\link{summary.mldr}}
#' @examples
#'
#' library(mldr)
#'
#' df <- data.frame(matrix(rnorm(1000), ncol = 10))
#' df$Label1 <- c(sample(c(0,1), 100, replace = TRUE))
#' df$Label2 <- c(sample(c(0,1), 100, replace = TRUE))
#' mymldr <- mldr_from_dataframe(df, labelIndices = c(11, 12), name = "testMLDR")
#'
#' summary(mymldr)
#'
#' @import stats
#' @export
mldr_from_dataframe <- function(dataframe, labelIndices, attributes, name) {
  if(!is.data.frame(dataframe))
    stop(paste(substitute(dataframe), "is not a valid data.frame"))

  if(missing(labelIndices))
    stop("labelIndices parameter is compulsory")

  if(!is.numeric(labelIndices))
    stop(paste(substitute(labelIndices), "is not a numeric vector"))

  new_mldr <- list()
  new_mldr$name <- if(missing(name)) substitute(dataframe) else name
  new_mldr$dataset <- dataframe

  #new_mldr$attributesIndexes <- 1:length(dataframe)
  #new_mldr$attributesIndexes <- new_mldr$attributesIndexes[! new_mldr$attributesIndexes %in% labelIndices]
  new_mldr$attributesIndexes <- which(!1:length(dataframe) %in% labelIndices)

  if (missing(attributes)) {
    new_mldr$attributes <- sapply(new_mldr$dataset, class)
    new_mldr$attributes[labelIndices] <- "{0,1}"
    factorIndexes <- which(new_mldr$attributes == "character")

    if (length(factorIndexes > 0))
      new_mldr$attributes[factorIndexes] <- sapply(factorIndexes, function(idx)
          paste("{", paste(
              levels(as.factor(new_mldr$dataset[, idx])),
              collapse = ","
            ), "}", sep = ""
          )
        )
  } else {
    new_mldr$attributes <- attributes
  }

  new_mldr$labels <- label_measures(dataframe, labelIndices)
  new_mldr$labelsets <- if(nrow(dataframe) > 0)
                          sort(table(as.factor(do.call(paste, c(dataframe[, new_mldr$labels$index], sep = "")))))
                        else
                          array()
  new_mldr <- dataset_measures(new_mldr)
  new_mldr$measures <- measures(new_mldr)

  class(new_mldr) <- "mldr"

  new_mldr
}
