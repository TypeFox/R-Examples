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



#' Validate if a data frame column confirms to primary key/ID constraints
#' 
#' \code{ValidatePrimKey} checks if a column in a data frame confirms to the 
#' primary key/ID constraints of absence of duplicates and NULL values. Aberrant
#' records if encountered are returned in the output list.
#' 
#' The function checks whether a field(column) in a data frame of PGR passport 
#' database confirms to the primary key/ID constraints of absence of duplicates 
#' and NULL values. If records with nonconforming values in the column are 
#' encountered, they are returned in the output list for rectification.
#' 
#' If multiple fields(columns) are given as a character vector in 
#' \code{prim.key} field, only the first element will be considered as the 
#' primary key/ID field(column).
#' 
#' Cleaning of the data in the input field(column) using the 
#' \code{\link[PGRdup]{DataClean}} function with appropriate arguments is 
#' suggested before running this function.
#' 
#' It is recommended to run this function and rectify aberrant records in a PGR 
#' passport database before creating a KWIC index using the 
#' \code{\link[PGRdup]{KWIC}} function.
#' 
#' @param x A data frame.
#' @param prim.key A character vector indicating the name of the data frame 
#'   column to be validated for primary key/ID constraints (see 
#'   \strong{Details}).
#' @return A list with containing the following components: \tabular{ll}{ 
#'   \code{message1} \tab Indicates whether duplicated values were encountered 
#'   in \code{prim.key} field(column) of data frame \code{x} or not. \cr 
#'   \code{Duplicates} \tab A data frame of the records with duplicated prim.key
#'   values if they were encountered. \cr \code{message2} \tab Indicates whether
#'   NULL values were encountered in \code{prim.key} field(column) of data frame
#'   \code{x} or not. \cr \code{NullRecords} \tab A data frame of the records 
#'   with NULL prim.key values if they were encountered. \cr }
#' @seealso \code{\link[PGRdup]{DataClean}}, \code{\link[PGRdup]{KWIC}}
#' @examples
#' GN <- GN1000
#' ValidatePrimKey(x=GN, prim.key="NationalID")
#' \dontrun{
#' # Show error in case of duplicates and NULL values 
#' # in the primary key/ID field "NationalID"
#' GN[1001:1005,] <- GN[1:5,]
#' GN[1001,3] <- ""
#' ValidatePrimKey(x=GN, prim.key="NationalID")}
#' @export
ValidatePrimKey <- function(x, prim.key) {
  if (is.data.frame(x) == FALSE) {
    # Check if x is a data frame and stop if not
    stop("x is not a data frame")
  }
  if (is.vector(prim.key) == FALSE) {
    # Check if prim.key is a vector or not
    stop("prim.key is not a vector")
  }
  if (length(prim.key) > 1) {
    # Check if only one field is given as input and use first element if not
    prim.key <- prim.key[1]
    warning("prim.key length >1; Only the first element is used")
  }
  if (is.element("FALSE", prim.key %in% colnames(x)) == TRUE) {
    # Check if prim.key field is present in x and stop if not
    stop("prim.key field missing in x")
  }
  Result <- list(message1 = NULL, Duplicates = NULL, message2 = NULL,
                 NullRecords = NULL)
  # Convert NAs to empty strings
  x[prim.key][is.na(x[,prim.key])] <- ""
  if (is.element("TRUE", duplicated(x[prim.key]))) {
    # Check if duplicated records are there in prim.key
    Result$message1 <- "ERROR: Duplicated records found in prim.key field"
    message(Result$message1)
    x$primdup <- duplicated(x[prim.key]) | duplicated(x[prim.key],
                                                      fromLast = TRUE)
    Result$Duplicates <- subset(x, primdup == TRUE)
    Result$Duplicates[, "primdup"] <- NULL
    Result$Duplicates <- Result$Duplicates[order(Result$Duplicates[prim.key]),]
  } else {
    Result$message1 <- "OK: No duplicated records found in prim.key field"
    message(Result$message1)
  }
  if (is.element("", x[,prim.key]) == TRUE) {
    # Check if empty characters are present in prim.key field
    Result$message2 <- "ERROR: NULL records found in prim.key field"
    message(Result$message2)
    Result$NullRecords <- subset(x, get(prim.key) == "")
  } else {
    Result$message2 <- "OK: No NULL records found in prim.key field"
    message(Result$message2)
  }
  return(Result)
}
