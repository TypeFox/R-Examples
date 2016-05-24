##
#' Generic method to generate an APA style table with descriptives for MS Word.
#'
#' @param  data Raw dataset with variables.
#' @param  variables The variable names for in the table.
#' @param  report (optional) Specify which descriptive statistics to report. Use a subset from \code{c("M", "SD", "r")}.
#' @param  title (optional) Name of the table.
#' @param  filename (optional) Specify the filename (including valid '\code{.docx}' extension).
#' @param  note (optional) Add a footnote to the bottom of the table.
#' @param  position (optional) Specify whether the correlations should be displayed in the \code{upper}, or \code{lower} diagonal of the table.
#' @param  merge (optional) Set (\code{TRUE}) if the mean and standard deviation columns should be merged into one column.
#' @param  landscape (optional) Set (\code{TRUE}) if the table should be generated in landscape mode.
#' @param  save (optional) Set (\code{FALSE}) if the table should not be saved in a document.
#' @return \code{apa.descriptives} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{save}{flag which indicates whther the document is saved}
#' \item{data}{dataset with descriptive statistics}
#' \item{table}{\code{FlexTable {ReporteRs}} object}
#' @importFrom "stats" "sd"
#' @export
#'
#' @examples
#'
#' # Use apa.descriptives function
#' apa.descriptives(
#'   data = data.frame(
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1)
#'   ),
#'   variables = c("Column 1", "Column 2", "Column 3", "Column 4")
#' )
##
apa.descriptives = function(data=data.frame(), variables=NULL, report="", title="APA Table", filename="APA Table.docx", note=NULL, position="lower", merge=FALSE, landscape=FALSE, save=TRUE) UseMethod("apa.descriptives")

##
#' Default method to generate an APA style table with descriptives for MS Word.
#'
#' @param  data Raw dataset with variables.
#' @param  variables The variable names for in the table.
#' @param  report (optional) Specify which descriptive statistics to report. Use a subset from \code{c("M", "SD", "r")}.
#' @param  title (optional) Name of the table.
#' @param  filename (optional) Specify the filename (including valid \code{.docx} extension).
#' @param  note (optional) Add a footnote to the bottom of the table.
#' @param  position (optional) Specify whether the correlations should be displayed in the \code{upper}, or \code{lower} diagonal of the table.
#' @param  merge (optional) Set (\code{TRUE}) if the mean and standard deviation columns should be merged into one column.
#' @param  landscape (optional) Set (\code{TRUE}) if the table should be generated in landscape mode.
#' @param  save (optional) Set (\code{FALSE}) if the table should not be saved in a document.
#' @return \code{apa.descriptives} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{save}{flag which indicates whther the document is saved}
#' \item{data}{dataset with descriptive statistics}
#' \item{table}{\code{FlexTable {ReporteRs}} object}
#' @importFrom "stats" "sd"
#' @export
#'
#' @examples
#'
#' # Use apa.descriptives function
#' apa.descriptives(
#'   data = data.frame(
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1),
#'     rnorm(100, mean = 0, sd = 1)
#'   ),
#'   variables = c("Column 1", "Column 2", "Column 3", "Column 4")
#' )
##
apa.descriptives.default = function(data=data.frame(), variables=NULL, report="", title="APA Table", filename="APA Table.docx", note=NULL, position="lower", merge=FALSE, landscape=FALSE, save=TRUE) {

  est = apaStyleDescriptives(data, variables, report, title, filename, note, position, merge, landscape, save)
  est$call = match.call()
  class(est) = "apa.descriptives"
  est

}

##
#' Define a print method
#'
#' @param  x A \code{apa.descriptives} object
#' @export
##
print.apa.descriptives = function(x, ...) {
  if(x$succes == TRUE) {
    cat("\n")
    if (x$save == TRUE) {
      cat("Word document succesfully generated in: ")
      cat(getwd())
    } else {
      cat("Succesfully generated the APA table")
    }
    cat("\n\n")
  }
}

# The main function
apaStyleDescriptives = function(data, variables, report, title, filename, note, position, merge, landscape, save) {

  # Initialize function
  options(warn = 0)

  # Check if a valid data frame is supplied
  if ((!is.data.frame(data)) || (is.data.frame(data) && nrow(data) == 0)) {
    error = "Invalid data is supplied."
    warning(error)
    return(list(succes = error))
  }

  # Define variables
  apa.report = c("M", "SD", "r")
  apa.data = data.frame(rep(NA, ncol(data)))

  # Check if valid variable names are supplied
  if(!is.character(variables)) {
    error = "No valid variable names are specified."
    warning(error)
    return(list(succes = error))
  } else {
    # Check if the length of the specified variables correspond with the length of the data
    if (length(variables) != length(data)) {
      error = "The supplied data doesn't match the specified number of variables."
      warning(error)
      return(list(succes = error))
    }
  }

  # Check if the items to report are valid
  if((!is.character(report))  || ("" %in% report && length(report) > 1)) {
    error = "The specified descriptives to report are not valid."
    warning(error)
    return(list(succes = error))
  } else {
    if (!"" %in% report) {
      if ((length(report) > 3) || (all(report %in% apa.report) == FALSE)) {
        error = "The specified descriptives to report are not valid. Only 'M', 'SD', and 'r' are allowed."
        warning(error)
        return(list(succes = error))
      } else {
        apa.report = report
      }
    }
  }

  # Check if a valid filename is supplied
  if((!is.character(filename)) || (!grepl(".docx", filename))) {
    error = "The supplied filename is not valid. Please specify a valid 'docx' file."
    warning(error)
    return(list(succes = error))
  } else {
    apa.filename = filename
  }

  # Check if a valid correlation matrix position is specified
  if((!is.character(position)) || (length(position) > 1) || (!"upper" %in% position && !"lower" %in% position)) {
    error = "The supplied display position for the correlation matrix is not valid. Only 'upper' or 'lower' position is allowed."
    warning(error)
    return(list(succes = error))
  } else {
    apa.position = position
  }

  # Check if the merge argument is a valid type
  if(!is.logical(merge)) {
    error = "The merge argument is not of logical type."
    warning(error)
    return(list(succes = error))
  } else {
    if ((merge == TRUE) && (!"M" %in% apa.report || !"SD" %in% apa.report)) {
      error = "Can not merge the mean and standard deviation into one column if they're not specified."
      warning(error)
      return(list(succes = error))
    }
  }

  # Check if the landscape argument is a valid type
  if(!is.logical(landscape)) {
    error = "The landscape argument is not of logical type."
    warning(error)
    return(list(succes = error))
  }

  # Check if the save argument is a valid type
  if(!is.logical(save)) {
    error = "The save argument is not of logical type."
    warning(error)
    return(list(succes = error))
  }

  # Check the size of the dataset
  if (ncol(data) > 20) {
    error = "The supplied data has too many variables to generate an APA formatted table."
    warning(error)
    return(list(succes = error))
  } else {

    # Calculate values for in the table
    if ("M" %in% apa.report) {
      apa.mean = apply(data, 2, function(x) mean(as.vector(x), na.rm = TRUE))
      if (merge == FALSE) {
        apa.data$mean = sprintf("%3.2f", round(apa.mean, digits = 2))
        header.m = c("M")
      } else {
        header.m = NULL
      }

    }

    if ("SD" %in% apa.report) {
      apa.sd = apply(data, 2, function(x) stats::sd(as.vector(x), na.rm = TRUE))
      if (merge == FALSE) {
        apa.data$sd = sprintf("%3.2f", round(apa.sd, digits = 2))
        header.sd = c("SD")
      } else {
        header.sd = NULL
      }
    }

    if (merge == TRUE) {
      apa.data$merge = apa.merge(apa.mean, apa.sd, c("M", "SD"))$data
      header.merge = apa.merge(apa.mean, apa.sd, c("M", "SD"))$header
    } else {
      header.merge = NULL
    }

    if ("r" %in% apa.report) {
      apa.r = as.data.frame(apa.cor.matrix(data, position = apa.position)$data, stringsAsFactors = FALSE)
      if (nrow(apa.r) != nrow(apa.data)) {
        error = "Arguments imply differing number of rows."
        warning(error)
        return(list(succes = error))
      }
      apa.data = data.frame(apa.data, apa.r)
      header.r = as.vector(rbind(c(1:ncol(data)), rep("*", ncol(data))))
    } else {
      header.r = NULL
    }

    colnames(apa.data)[1] = "Variable"

    # Put numbers in front of the variable names
    apa.data[[1]] = strsplit(paste(seq(1:length(variables)), ". ", variables, collapse = ";", sep = ""), ";")[[1]]
    apa.header = c("Variable", header.m, header.sd, header.merge, header.r)

    apa.table = apaStyle::apa.table(data = data.frame(apa.data), level1.header = apa.header, title = title, filename = filename, note = note, landscape = landscape, save = save)

    # Check if the document was succesfully generated
    if(apa.table$succes != TRUE) {
      return(list(succes = apa.table$succes))
    }

    return(list(succes = TRUE, save = save, data = data.frame(apa.data), table = apa.table$table))

  }

}
