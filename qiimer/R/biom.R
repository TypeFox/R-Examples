#' Extract raw data from a BIOM object.
#' 
#' This function extracts the raw data from a BIOM object, using the correct 
#' row and column names in the result.
#' 
#' The BIOM object can be any list-like representation of the JSON 
#' source code in a BIOM-format file produced by QIIME.  There are several 
#' options for creating BIOM objects from QIIME output files.  The official
#' library for BIOM files, \code{biom}, can create compatible objects via the 
#' \code{read_biom} function.  Alternately, the \code{fromJSON} function from
#' either \code{RJSONIO} or \code{rjson} may be used.
#' 
#' @param b A BIOM object.
#' @return For sparse biom objects, returns a 3-column data frame of row names,
#'   column names, and the data value.  The first column is named using the 
#'   first word in the BIOM object's \code{type} attribute (e.g. "OTU" for OTU
#'   tables).  The second and third columns are named "SampleID" and "value",
#'   respectively.  For dense biom objects, returns a matrix.
#' @export
#' @examples
#' data(relmbeta_biom)
#' head(biom_raw_data(relmbeta_biom))
biom_raw_data <- function (b) {
  # Retrieve the first word from the biom type
  # One of: OTU, Taxon, Gene, Function, Ortholog, Pathway, Metabolite
  row_label <- strsplit(b$type, " ")[[1]][1]
  rnames <- sapply(b$rows, `[[`, "id")
  cnames <- sapply(b$columns, `[[`, "id")
  if (b$matrix_type == "sparse") {
    df <- as.data.frame(do.call(rbind, b$data))
    colnames(df) <- c(row_label, "SampleID", "value")
    df[[row_label]] <- factor(rnames[df[[row_label]] + 1])
    df$SampleID <- factor(cnames[df$SampleID + 1])
    df
  } else if (b$matrix_type == "dense") {
    dnames <- list()
    dnames[[row_label]] <- rnames
    dnames[["SampleID"]] <- cnames
    matrix(unlist(b$data), nrow=b$shape[1], byrow=T, dimnames=dnames)
  } else {
    message("unsupported matrix type")
    NULL
  }
}

#' Extract taxonomy info from a biom object.
#'
#' For BIOM objects representing OTU tables, this function will extract the 
#' taxonomy for each OTU.  The taxonomy info is typically stored in an 
#' attribute of the row metadata named "taxonomy".
#'
#' This function may be used more generally to extract metadata from rows in a 
#' BIOM object.  The \code{attr} argument may be adjusted to match the desired
#' attribute name.  If the metadata has a nested structure, a character vector 
#' may be supplied to the \code{attr} argument.
#'
#' The BIOM object can be any list-like representation of the JSON 
#' source code in a BIOM-format file produced by QIIME.  There are several 
#' options for creating BIOM objects from QIIME output files.  The official
#' library for BIOM files, \code{biom}, can create compatible objects via the 
#' \code{read_biom} function.  Alternately, the \code{fromJSON} function from
#' either \code{RJSONIO} or \code{rjson} may be used.
#' 
#' @param b A BIOM object.
#' @param attr The metadata attribute under which the taxonomy information 
#'   can be found for each row item in the biom file.
#' @return A list of character vectors, one per row.
#' @export
#' @examples
#' data(relmbeta_biom)
#' head(biom_taxonomy(relmbeta_biom))
biom_taxonomy <- function (b, attr="taxonomy") {
  bt <- lapply(b$rows, `[[`, c("metadata", attr))
  names(bt) <- sapply(b$rows, `[[`, "id")
  bt
}
