
# Functions to export a mldr object to several formats

#' Export an mldr object or set of mldr objects to different file formats
#' @description Writes one or more files in the specified formats with the content of the \code{mldr} or \code{mldr.folds} given as parameter
#' @param mld The \code{mldr/mldr.folds} object to be exported
#' @param format A vector of strings stating the desired file formats. It can contain the values \code{'MULAN'}, \code{'MEKA'},
#' \code{'KEEL'}, \code{'CSV'} and \code{'LIBSVM'}
#' @param sparse Boolean value indicating if sparse representation has to be used for ARFF-based file formats
#' @param basename Base name for the files. \code{'unnamed_mldr'} is used by default
#' @examples
#'\dontrun{
#' library(mldr.datasets)
#' write.mldr(emotions, format = c('CSV', 'KEEL'))
#' }
#' @export
write.mldr <- function(mld, format = c("MULAN", "MEKA"), sparse = FALSE, basename = ifelse(!is.null(mld$name) && nchar(mld$name) > 0,
                                                                                           regmatches(mld$name, regexpr("(\\w)+", mld$name)),
                                                                                           "unnamed_mldr")) {
  format <- toupper(format)
  available.formats <- c("MULAN", "MEKA", "KEEL", "CSV", "LIBSVM")

  # Parameter checks
  if (!all(format %in% available.formats)) {
    stop("Invalid format found. Allowed formats: ", paste(available.formats, collapse = ", "))
  }
  if (!"mldr" %in% class(mld)) {
    if ("mldr.folds" %in% class(mld))
      invisible(lapply(1:length(mld), function(i) {
        write.mldr(mld[[i]]$train, format, sparse, basename = paste0(basename, "-", i, 'x', length(mld), '-tra') )
        write.mldr(mld[[i]]$test, format, sparse, basename = paste0(basename, "-", i, 'x', length(mld), '-test') )
      }))
    else
      stop("Object must be of class mldr or mldr.folds")
  } else {

    inform <- function(name) cat(paste0("Wrote file ", name, "\n"))

    if ("MULAN" %in% format) {
      # Open and write ARFF file
      if (!"MEKA" %in% format) {
        name <- paste0(basename, ".arff")
        arffConnection <- file(name)
        writeLines(export.mulan(mld, sparse), arffConnection)
        close(arffConnection)
        inform(name)
      }

      # Open and write XML file
      name <- paste0(basename, ".xml")
      xmlConnection <- file(name)
      writeLines(export.xml(mld), xmlConnection)
      close(xmlConnection)
      inform(name)
    }

    if ("MEKA" %in% format) {
      # Open and write ARFF file
      name <- paste0(basename, ".arff")
      arffConnection <- file(name)
      writeLines(export.meka(mld, sparse), arffConnection)
      close(arffConnection)
      inform(name)
    }

    if ("KEEL" %in% format) {
      # Open and write DAT file
      name <- paste0(basename, ".dat")
      datConnection <- file(name)
      writeLines(export.keel(mld, sparse), datConnection)
      close(datConnection)
      inform(name)
    }

    if ("CSV" %in% format) {
      # Open and write CSV file
      name <- paste0(basename, ".csv")
      csvConnection <- file(name)
      writeLines(export.csv(mld, sparse), csvConnection)
      close(csvConnection)
      inform(name)

      name <- paste0(basename, "_labels.csv")
      labelConnection <- file(name)
      writeLines(export.csv.labels(mld), labelConnection)
      close(labelConnection)
      inform(name)
    }

    if ("LIBSVM" %in% format) {
      # Open and write SVM file
      name <- paste0(basename, ".svm")
      svmConnection <- file(name)
      writeLines(export.libsvm(mld), svmConnection)
      close(svmConnection)
      inform(name)
    }
  }
}

export.mulan <- function(mld, sparse) {
  paste(
    export.mulan.header(mld),
    export.arff.attributes(mld),
    export.arff.data(mld, sparse),
    sep = "\n"
  )
}

export.meka <- function(mld, sparse) {
  paste(
    export.meka.header(mld),
    export.arff.attributes(mld),
    export.arff.data(mld, sparse),
    sep = "\n"
  )
}

export.keel <- function(mld, sparse) {
  paste(
    export.keel.header(mld),
    export.arff.attributes(mld),
    export.arff.inputs(mld),
    export.arff.outputs(mld),
    export.arff.data(mld, sparse),
    sep = "\n"
  )
}

export.mulan.header <- function(mld) {
  paste0("@relation ", ifelse(!is.null(mld$name) && nchar(mld$name) > 0,
                             mld$name,
                             "unnamed_mldr"))
}

export.meka.header <- function(mld) {
  paste0("@relation '", ifelse(!is.null(mld$name) && nchar(mld$name) > 0,
                             mld$name,
                             "unnamed_mldr"),
        ": -C ", mld$measures$num.labels, "'")
}

export.keel.header <- export.mulan.header

export.arff.attributes <- function(mld) {
  paste0("@attribute ", names(mld$attributes), " ", mld$attributes, collapse = "\n")
}

export.arff.inputs <- function(mld) {
  paste0(
    "@inputs ", do.call(
      paste,
      c(
        as.list(names(mld$attributes[mld$attributesIndexes])),
        sep = ", "
      )
    )
  )
}

export.arff.outputs <- function(mld) {
  paste0(
    "@outputs ", do.call(
      paste,
      c(
        as.list(rownames(mld$labels)),
        sep = ", "
      )
    )
  )
}

export.arff.data <- function(mld, sparse, header = "@data\n") {
  data <- mld$dataset[, 1:mld$measures$num.attributes]
  data[is.na(data)] <- '?'
  paste0(
    header,
    ifelse(sparse, export.sparse.arff.data(data), export.dense.arff.data(data))
  )
}


export.dense.arff.data <- function(data) {
  paste(
    do.call("paste", c(unname(data), list(sep = ','))),
    collapse = "\n"
  )
}

export.sparse.arff.data <- function(data) {
  paste(
    apply(
      # 'as.matrix' implicit conversion of a data.frame will insert spaces to adjust
      # width of values (when the inferred data type is 'character'). To prevent
      # this, a workaround needs to be done by manually formatting the data.frame.
      # Source: http://stackoverflow.com/a/15618761
      sapply(data, format, trim = TRUE, justify = "none"),
      1, function(instance)
        paste0("{",
               paste(which(instance != 0) - 1, # features start counting at 0
                     instance[instance != 0],
                     sep = " ", collapse = ","
                     ),
              "}"
        )
    ),
    collapse = "\n"
  )
}

export.csv <- function(mld, sparse) export.arff.data(mld, sparse, header = "")

export.csv.labels <- function(mld) {
  paste(rownames(mld$labels), mld$labels$index, sep = ", ", collapse = "\n")
}

export.xml <- function(mld) {
  xmlheader <- '<?xml version="1.0" encoding="utf-8"?>'
  labelstag <- '<labels xmlns="http://mulan.sourceforge.net/labels">'
  labeltags <- paste(c('<label name="'), rownames(mld$labels), c('"></label>'), collapse = "\n", sep = "")
  labelsend <- '</labels>'

  paste(xmlheader, labelstag, labeltags, labelsend, sep = "\n")
}

export.libsvm <- function(mld) {
  paste(
    apply(mld$dataset, 1, function(instance) {
        inputs <- instance[mld$attributesIndexes]
        outputs <- instance[mld$labels$index]
        paste(
          paste(which(outputs == 1) - 1, collapse = ","), # libSVM counts labels starting from zero
          paste(which(instance != 0), instance[instance != 0], sep = ":", collapse = " "),
          sep = " "
        )
      }
    ),
    collapse = "\n"
  )
}
