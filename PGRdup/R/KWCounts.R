### This file is part of 'PGRdup' package for R.

### Copyright (C) 2014, ICAR-NBPGR.
#
# PGRdup is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# PGRdup is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Generate keyword counts
#' 
#' \code{KWCounts} generates keyword counts from PGR passport database 
#' fields(columns).
#' 
#' This function computes the keyword counts from PGR passport database 
#' fields(columns) specified in the \code{fields} argument. The first field is 
#' considered as the primary key or identifier and is not used for counting the 
#' keywords. Any strings given in the \code{excep} argument are ignored for 
#' generating the counts.
#' 
#' The keyword counts can give a rough indication of the completeness of the 
#' data in the database fields being used for identification of probable 
#' duplicates.
#' 
#' @param x A data frame.
#' @param fields A character vector with the names of fields(columns) of the 
#'   data frame from which KWIC index is to be generated. The first field is 
#'   considered as the primary key or identifier (see \strong{Details}).
#' @param excep A vector of the keywords not to be considered for the counts
#'   (see \strong{Details}).
#' @examples
#' # Load PGR passport database
#' GN <- GN1000
#'
#' # Specify database fields to use as a vector
#' GNfields <- c("NationalID", "CollNo", "DonorID", "OtherID1", "OtherID2")
#' 
#' # Specify the exceptions as a vector
#' exep <- c("A", "B", "BIG", "BOLD", "BUNCH", "C", "COMPANY", "CULTURE", 
#'          "DARK", "E", "EARLY", "EC", "ERECT", "EXOTIC", "FLESH", "GROUNDNUT", 
#'          "GUTHUKAI", "IMPROVED", "K", "KUTHUKADAL", "KUTHUKAI", "LARGE", 
#'          "LIGHT", "LOCAL", "OF", "OVERO", "P", "PEANUT", "PURPLE", "R", 
#'          "RED", "RUNNER", "S1", "SAM", "SMALL", "SPANISH", "TAN", "TYPE", 
#'          "U", "VALENCIA", "VIRGINIA", "WHITE")
#'          
#' # Compute the keyword counts
#' GNKWCouts <- KWCounts(GN, GNfields, exep)
#' 
#' # Plot the keyword counts
#' bp <- barplot(table(GNKWCouts$COUNT),
#'               xlab = "Word count", ylab = "Frequency", col = "#CD5555")
#' text(bp, 0, table(GNKWCouts$COUNT),cex=1,pos=3)
#' @seealso \code{\link[stringi]{stri_count}}
#' @return A data frame with the keyword counts for each record.
#' @import data.table
#' @importFrom stringi stri_count
#' @export KWCounts
KWCounts <- function(x, fields, excep) {
  if (is.data.frame(x) == FALSE) {
    # Check if x is a data.frame and stop if not
    stop("x is not a data.frame")
  }
  if (is.vector(fields) == FALSE) {
    # Check if fields is a vector or not
    stop("fields is not a vector")
  }
  if (length(fields) == 1) {
    # Check if more than one field is given as input and stop if not
    stop("Only one field given as input")
  }
  if (is.element(FALSE, fields %in% colnames(x)) == TRUE) {
    # Check if fields are present in x and stop if not
    stop("One or more fields are missing in x")
  }
  # Check excep argument
  if (!is.null(excep) && is.vector(excep, mode = "character") == FALSE) {
    stop("'excep' is not a character vector")
  }
  if (!is.null(excep)) {
    excep <- toupper(excep)
  } else {
    excep <- ""
  }
  #setDT(x)
  x <- as.data.table(x[fields])
  # Convert the fields in x to character
  for (col in fields) set(x, j = col, value = as.character(x[[col]]))
  # Convert NAs to empty strings
  for (j in fields) set(x , which(is.na(x[[j]])), j, "")
  setDF(x)
  if (is.element("", x[fields[1]]) | is.element(TRUE,
                                                duplicated(x[fields[1]]))) {
    # Check primary key/ID is unique and not NULL
    stop("Primary key/ID field should be unique and not NULL\n Use PGRdup::ValidatePrimKey() to identify and rectify the aberrant records first")
  }
  #setDT(x)
  x <- as.data.table(x)
  # Remove exceptions
  x[, fields[-1] := lapply(.SD, function(x) gsub(paste0(excep,
                                                        collapse = "|"), "",
                                                 x)), .SDcols = fields[-1]]
  # Get the word counts
  x[, COUNT := stri_count(do.call(paste, c(.SD, sep = " ")),
                          regex = "\\S+"), .SDcols = fields[-1]]
  # Prepare output
  x[, fields[-1] := NULL]
  setnames(x, names(x), c("PRIM_ID", "COUNT"))
  setDF(x)
  return(x)
}
