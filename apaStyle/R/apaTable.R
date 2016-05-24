##
#' Generic method to generate an APA style table for MS Word
#'
#' @param  data Dataset with statistics.
#' @param  level1.header The column names for the first header in the table.
#' @param  level1.colspan (optional) The colspan for the first header column.
#' @param  level2.header (optional) The column names for the second header in the table.
#' @param  number (optional) The table number in the document.
#' @param  title (optional) Name of the table.
#' @param  filename (optional) Specify the filename (including valid '\code{.docx}' extension).
#' @param  note (optional) Add a footnote to the bottom of the table.
#' @param  landscape (optional) Set (\code{TRUE}) if the table should be generated in landscape mode.
#' @param  save (optional) Set (\code{FALSE}) if the table should not be saved in a document.
#' @details
#'
#' This method can generate tables with two headers. If two headers are required, it is necesary to
#' specifify the colspan for the upper level (\code{level1.colspan}). If only one header is required
#' only the header items need to be specified for \code{level1.header}, and \code{level1.colspan} and
#' \code{level2.header} do not need be specified.
#'
#' This method allows users to specify a column in which either the level of significance (header:
#' \code{"*"}), or a subscript (header: \code{"_"}) is given. For example, when there is a column
#' with a F-value and there shouldn't be an additional column with the corresponding p-values, it
#' is possible to specify an additional column with significant values (i.e., +p < .10; *p < .05;
#' **p < .01; ***p < .001) which will be merged as one column in the final table.
#'
#' Often it is necesary to provide a table with the means from different groups or conditions. Using
#' the subscript header (\code{"_"}) it is possible to supply a column with subscripts which indicates
#' which means on a row significantly differ from each other.
#'
#' @return \code{apa.table} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{save}{flag which indicates whether the document is saved}
#' \item{table}{\code{FlexTable {ReporteRs}} object}
#' @importFrom "ReporteRs" "docx" "writeDoc" "FlexTable" "addHeaderRow" "addSection" "addParagraph" "cellProperties" "parProperties" "textProperties" "pot" "chprop" "textNormal" "parLeft" "parRight" "parCenter"
#' @importFrom "utils" "head"
#' @export
#'
#' @examples
#'
#' # Use apa.table function with a minimum of parameters
#' # Specify statistics
#' example <- data.frame(
#'   c("Column 1", "Column 2", "Column 3"),
#'   c(3.45, 5.21, 2.64),
#'   c(1.23, 1.06, 1.12)
#' )
#'
#' # Create table
#' apa.table(data = example, level1.header = c("Variable", "M", "SD"))
#'
#' # Create a table with two headers
#' # Specify statistics
#' example <- data.frame(
#'   c("Column 1", "Column 2", "Column 3"),
#'   c(3.45, 5.21, 2.64),
#'   c(1.23, 1.06, 1.12),
#'   c(8.22, 25.12, 30.27),
#'   c("+", "**", "***")
#' )
#'
#' # Run method and preview table
#' apa.table(
#'   data = example,
#'   level1.header = c("", "Descriptives", "Inferential"),
#'   level1.colspan = c(1, 2, 2),
#'   level2.header = c("Variable", "M", "SD", "t-value", "*")
#' )$table
##
apa.table = function(data=data.frame(), level1.header=NULL, level1.colspan=NULL, level2.header=NULL, number="XX", title="APA Table", filename="APA Table.docx", note=NULL, landscape=FALSE, save=TRUE) UseMethod("apa.table")

##
#' Default method to generate an APA style table for MS Word
#'
#' @param  data Dataset with statistics.
#' @param  level1.header The column names for the first header in the table.
#' @param  level1.colspan (optional) The colspan for the first header column.
#' @param  level2.header (optional) The column names for the second header in the table.
#' @param  number (optional) The table number in the document.
#' @param  title (optional) Name of the table.
#' @param  filename (optional) Specify the filename (including valid '\code{.docx}' extension).
#' @param  note (optional) Add a footnote to the bottom of the table.
#' @param  landscape (optional) Set (\code{TRUE}) if the table should be generated in landscape mode.
#' @param  save (optional) Set (\code{FALSE}) if the table should not be saved in a document.
#' @details
#'
#' This method can generate tables with two headers. If two headers are required, it is necesary to
#' specifify the colspan for the upper level (\code{level1.colspan}). If only one header is required
#' only the header items need to be specified for \code{level1.header}, and \code{level1.colspan} and
#' \code{level2.header} do not need be specified.
#'
#' This method allows users to specify a column in which either the level of significance (header:
#' \code{"*"}), or a subscript (header: \code{"_"}) is given. For example, when there is a column
#' with a F-value and there shouldn't be an additional column with the corresponding p-values, it
#' is possible to specify an additional column with significant values (i.e., +p < .10; *p < .05;
#' **p < .01; ***p < .001) which will be merged as one column in the final table.
#'
#' Often it is necesary to provide a table with the means from different groups or conditions. Using
#' the subscript header (\code{"_"}) it is possible to supply a column with subscripts which indicates
#' which means on a row significantly differ from each other.
#'
#' @return \code{apa.table} object; a list consisting of
#' \item{succes}{message in case of an error}
#' \item{save}{flag which indicates whether the document is saved}
#' \item{table}{\code{FlexTable {ReporteRs}} object}
#' @importFrom "ReporteRs" "docx" "writeDoc" "FlexTable" "addHeaderRow" "addSection" "addParagraph" "cellProperties" "parProperties" "textProperties" "pot" "chprop" "textNormal" "parLeft" "parRight" "parCenter"
#' @importFrom "utils" "head"
#' @export
#'
#' @examples
#'
#' # Use apa.table function with a minimum of parameters
#' # Specify statistics
#' example <- data.frame(
#'   c("Column 1", "Column 2", "Column 3"),
#'   c(3.45, 5.21, 2.64),
#'   c(1.23, 1.06, 1.12)
#' )
#'
#' # Create table
#' apa.table(data = example, level1.header = c("Variable", "M", "SD"))
#'
#' # Create a table with two headers
#' # Specify statistics
#' example <- data.frame(
#'   c("Column 1", "Column 2", "Column 3"),
#'   c(3.45, 5.21, 2.64),
#'   c(1.23, 1.06, 1.12),
#'   c(8.22, 25.12, 30.27),
#'   c("+", "**", "***")
#' )
#'
#' # Run method and preview table
#' apa.table(
#'   data = example,
#'   level1.header = c("", "Descriptives", "Inferential"),
#'   level1.colspan = c(1, 2, 2),
#'   level2.header = c("Variable", "M", "SD", "t-value", "*")
#' )$table
##
apa.table.default = function(data=data.frame(), level1.header=NULL, level1.colspan=NULL, level2.header=NULL, number="XX", title="APA Table", filename="APA Table.docx", note=NULL, landscape=FALSE, save=TRUE) {

  est = apaTable(data, level1.header, level1.colspan, level2.header, number, title, filename, note, landscape, save)
  est$call = match.call()
  class(est) = "apa.table"
  est

}

##
#' Define a print method
#'
#' @param  x A \code{apa.table} object
#' @export
##
print.apa.table = function(x, ...) {
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
apaTable = function(data, level1.header, level1.colspan, level2.header, number, title, filename, note, landscape, save) {

  # Initialize function
  options(warn = 0)

  # Define variables
  level2 = FALSE
  apa.signif = NULL
  apa.italics = c("B", "d", "df", "F", "M", "n", "N", "p", "r", "R^2", "SD", "SE", "t", "z" )
  apa.tableName = ifelse(is.numeric(number), paste("Table", number, sep = "", collapse = ""), "Table XX")

  # Check if a valid data frame is supplied
  if ((!is.data.frame(data)) || (is.data.frame(data) && nrow(data) == 0)) {
    error = "Invalid data is supplied."
    warning(error)
    return(list(succes = error))
  }

  # Check if valid headers are supplied
  if(!is.character(level1.header)) {
    error = "No valid headers are specified."
    warning(error)
    return(list(succes = error))
  }

  # Check if the save argument is a valid type
  if(!is.logical(save)) {
    error = "The save argument is not of logical type."
    warning(error)
    return(list(succes = error))
  } else {

    if (save == TRUE) {

      # Check if a valid filename is supplied
      if((!is.character(filename)) || (!grepl(".docx", filename))) {
        error = "The supplied filename is not valid. Please specify a valid 'docx' file."
        warning(error)
        return(list(succes = error))
      } else {
        apa.filename = filename
      }

      # Check if the landscape argument is a valid type
      if(!is.logical(landscape)) {
        error = "The landscape argument is not of logical type."
        warning(error)
        return(list(succes = error))
      }

    }
  }

  # Check the size of the dataset
  if (ncol(data) > 20 | nrow(data) > 100) {
    error = "The supplied data is too big to generate an APA formatted table."
    warning(error)
    return(list(succes = error))
  } else {

    # Convert factors to characters
    i = sapply(data, is.factor)
    data[i] = lapply(data[i], as.character)

    # Prepare table headers

    # Check if level 2 headers are supplied
    if(is.character(level2.header)) {
      if((is.null(level1.colspan)) || (sum(level1.colspan) != length(level2.header))) {
        error = "The level 1 colspan doesn't match the number of level 2 headers."
        warning(error)
        return(list(succes = error))
      } else {

        level2 = TRUE

        # Internal function to save where to insert empty columns
        sequence = function(x) {
          n = length(x)
          y = x[-1L] != x[-n] + 1
          i = c(which(y|is.na(y)),n)
          list(lengths = diff(c(0L,i)), values = x[utils::head(c(0L,i)+1L,-1L)])
        }

        # Insert empty columns between level1 headers
        test1 = sequence(grep("\\w|\\S", level1.header, perl = TRUE))

        if (test1$values[1] > 1) {
          start = paste(level1.header[1:test1$values[[1]]-1], sep = "", collapse = ";")
        }

        tmp.header1 = tmp.colspan = tmp.header2 = NULL

        for (i in 1:length(test1$values)) {
          n = test1$values[i] + (test1$lengths[i] - 1)
          tmp.header1[i] = paste(level1.header[test1$values[i]:n], sep = "", collapse = "; ;")
        }
        tmp.header1 = paste(tmp.header1, sep = "", collapse = ";;")
        tmp.header1 = c(start, tmp.header1)
        tmp.header1 = unlist(strsplit(paste(tmp.header1, sep = "", collapse = ";"), split=";"))

        # Insert empty columns for level 2 headers
        test2 = first = last = grep(" ", tmp.header1, perl = TRUE)

        for (i in 1:length(test2)) {
          if (i == 1) {
            tmp.colspan = c(level1.colspan[1:test2[i]-1], 1, level1.colspan[test2[i]:length(level1.colspan)])
            tmp.header2 = c(level2.header[1:sum(tmp.colspan[1:test2[i]-1])], " ", level2.header[sum(tmp.colspan[1:test2[i]]):length(level2.header)])
          } else {
            tmp.colspan = c(tmp.colspan[1:first[i]], 1, level1.colspan[last[i]:length(level1.colspan)])
            tmp.header2 = c(tmp.header2[1:sum(tmp.colspan[1:first[i]])], " ", level2.header[sum(tmp.colspan[1:last[i]]):length(level2.header)])
          }
          first = first - 1
          last = last - i
        }

        # Save new generated headers
        level1.colspan = tmp.colspan
        level1.header = tmp.header1
        level2.header = tmp.header2

        if (sum(level1.colspan) != length(level2.header)) {
          error = "The generated level 1 colspan doesn't match the number of level 2 headers."
          warning(error)
          return(list(succes = error))
        }

        apa.header = header = level2.header

      }
    } else {
      apa.header = header = level1.header
    }

    signif = which(header == "*", arr.ind = TRUE) - 1
    script = which(header == "_", arr.ind = TRUE) - 1
    bridge = which(header == " ", arr.ind = TRUE)

    # Check significance columns
    if(length(signif) > 0) {

      # Convert "+" symbol to unicode dagger symbol
      data[which(data == "+", arr.ind = TRUE)] = "\u2020"

      # Create a footnote indicating significant values
      if (save == TRUE) {
        apa.signif = apaStyle::apa.signif(data)$signif
      }
    }

    # Create a user defined footnote
    if(!is.null(note) && save == TRUE) {
      apa.note = ReporteRs::pot("Note. ", ReporteRs::textProperties(font.family = "Times", font.size = 12, font.style = "italic")) +
        ReporteRs::pot(note, ReporteRs::textProperties(font.family = "Times", font.size = 12))
    } else {
      apa.note = ""
    }

    # Include empty columns where column spaces are requested
    if (length(bridge) > 0) {
      apa.void = rep("", nrow(data))
      for(i in 1:length(bridge)) {
        if(bridge[i] == length(header)) {
          data = data.frame(data, apa.void)
        } else if (bridge[i] == 1) {
          data = data.frame(apa.void, data)
        } else {
          data = data.frame(data[1:bridge[i]-1], apa.void, data[(bridge[i]):ncol(data)])
        }
      }
    }

    # Check if the length of the dataframe matches with the length of the header
    if (ncol(data) != length(header)) {
      error = "The supplied data doesn't match the specified table header."
      warning(error)
      return(list(succes = error))
    } else {

      # Text default for the APA Table
      options('ReporteRs-fontsize' = 10, 'ReporteRs-default-font' = 'Times')

      # Define header properties
      headerCellProps = ReporteRs::cellProperties(padding = 7, border.bottom.style = "solid", border.top.style = "solid", border.left.style = "none", border.right.style = "none")
      headerParProps = ReporteRs::parProperties()

      # Create APA table
      apa.table = ReporteRs::FlexTable(data, header.columns = FALSE, header.text.props = ReporteRs::textNormal(), header.par.props = ReporteRs::parCenter(), header.cell.props = headerCellProps)
      apa.table[, 2:length(data)] = ReporteRs::parCenter()

      colspan = c()
      merged = 0
      skip = FALSE

      for(j in 1:length(data)) {

        if((length(signif) > 0 & signif[1] == j) | (length(script) > 0 & script[1] == j)) {

          apa.table[, j] = ReporteRs::parRight()
          apa.table[, j] = ReporteRs::chprop(ReporteRs::cellProperties(padding.right = 0, padding.left = 7, padding.top = 7, padding.bottom = 7, border.style = "none"))
          apa.table[nrow(data), j] = ReporteRs::chprop(ReporteRs::cellProperties(padding.right = 0, padding.left = 7, padding.top = 7, padding.bottom = 7, border.bottom.style = "solid", border.top.style = "none", border.left.style = "none", border.right.style = "none"))

          index = j + 1
          remove = index - merged

          colspan = c(colspan, 2)
          apa.header = apa.header[-remove]

          apa.table[, index] = ReporteRs::parLeft()
          apa.table[, index] = ReporteRs::chprop(ReporteRs::cellProperties( padding.left = 0, padding.right = 7, padding.top = 7, padding.bottom = 7, border.style = "none"))
          apa.table[nrow(data), index] = ReporteRs::chprop(ReporteRs::cellProperties( padding.left = 0, padding.right = 7, padding.top = 7, padding.bottom = 7, border.bottom.style = "solid", border.top.style = "none", border.left.style = "none", border.right.style = "none"))

          merged = merged + 1
          skip = TRUE

          if (length(signif) > 0 & signif[1] == j) {
            signif = signif[-1]
          } else {
            script = script[-1]
            apa.table[, index] = ReporteRs::textProperties(vertical.align = 'subscript')
          }

        } else {

          if (TRUE == skip) {
            skip = FALSE
          } else {
            colspan = c(colspan, 1)
            apa.table[, j] = ReporteRs::chprop(ReporteRs::cellProperties( padding = 7, border.style = "none"))
            apa.table[nrow(data), j] = ReporteRs::chprop(ReporteRs::cellProperties(padding = 7, border.bottom.style = "solid", border.top.style = "none", border.left.style = "none", border.right.style = "none"))
          }

        }
      }

      if (length(data) != sum(colspan)) {
        error = "The sum of colspan is different from the number of columns of the dataset."
        warning(error)
        return(list(succes = error))
      }

      if (level2 == TRUE) {

        if (length(level2.header) != sum(level1.colspan)) {
          error = "The sum of colspan is different from the number of columns of the dataset."
          warning(error)
          return(list(succes = error))
        }

        apa.table = ReporteRs::addHeaderRow(apa.table, value = level1.header, colspan = level1.colspan)
        apa.table = ReporteRs::addHeaderRow(apa.table, value = apa.header, colspan = colspan)

        apa.table[1, 1:length(level2.header), to = 'header', side = 'bottom'] = ReporteRs::borderNone()
        apa.table[2, 1:length(level2.header), to = 'header', side = 'top'] = ReporteRs::borderNone()

        borders = sapply(grep("\\w", level1.header, perl = TRUE), function(x) sum(level1.colspan[1:x-1])+1)
        catch = sapply(borders, function(x) apa.table[1, borders, to = 'header', side = 'bottom'] = ReporteRs::borderProperties(width = 1, style = 'solid'))

      } else {
        apa.table = ReporteRs::addHeaderRow(apa.table, value = apa.header, colspan = colspan)
      }

      apa.table[, 1, to = 'header'] = ReporteRs::parLeft()

      # Put APA reserved abbreviations in italic
      for(i in 1:length(header)) {
        if(is.element(header[i], apa.italics)) {
          row = ifelse(level2 == TRUE, 2, 1)
          apa.table[row, i, to = 'header'] = ReporteRs::textProperties(font.style = 'italic')
        }
      }

      if (save == TRUE) {

        # Generate MS Word document
        apa.doc = ReporteRs::docx(title = title)

        if(landscape == TRUE) {
          apa.doc = ReporteRs::addSection(apa.doc, landscape = landscape)
        }

        # Add content to word document
        apa.doc = ReporteRs::addParagraph(apa.doc, ReporteRs::pot(apa.tableName, ReporteRs::textProperties(font.family="Times", font.size=12)), stylename = "Normal")
        apa.doc = ReporteRs::addParagraph(apa.doc, ReporteRs::pot(title, ReporteRs::textProperties(font.family="Times", font.size=12, font.style="italic")), stylename = "Normal")
        apa.doc = ReporteRs::addFlexTable(apa.doc, apa.table)

        # Add table footers
        if(!is.null(apa.signif)) {
          apa.doc = ReporteRs::addParagraph(apa.doc, apa.signif, stylename = "Normal")
        }

        if(!is.null(apa.note)) {
          apa.doc = ReporteRs::addParagraph(apa.doc, apa.note, stylename = "Normal")
        }

        if(landscape == TRUE) {
          apa.doc = ReporteRs::addSection(apa.doc)
        }

        if (file.exists(apa.filename)) {
          if(!file.create(apa.filename, overwrite = TRUE, showWarnings = FALSE)[1]) {
            error = "The specified filename already exists and is used by another application. Make sure you close this application first."
            warning(error)
            return(list(succes = error))
          }
        }

        ReporteRs::writeDoc(apa.doc, apa.filename)

      }

      return(list(succes = TRUE, save = save, table = apa.table))

    }

  }

}
