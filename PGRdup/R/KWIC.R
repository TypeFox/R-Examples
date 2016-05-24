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



#' Create a KWIC Index
#' 
#' \code{KWIC} creates a Keyword in Context index from PGR passport database 
#' fields.
#' 
#' The function generates a Keyword in Context index from a data frame of a PGR 
#' passport database based on the fields(columns) stated in the arguments, using
#' \code{\link[data.table]{data.table}} package.
#' 
#' The first element of vector \code{fields} is considered as the primary key or
#' identifier which uniquely identifies all rows in the data frame.
#' 
#' Cleaning of the data the input fields(columns) using the 
#' \code{\link[PGRdup]{DataClean}} function with appropriate arguments is 
#' suggested before running this function.
#' 
#' @param x A data frame from which KWIC index is to be generated.
#' @param fields A character vector with the names of fields(columns) of the 
#'   data frame from which KWIC index is to be generated. The first field is 
#'   considered as the primary key or identifier (see \strong{Details}).
#' @param min.freq Frequency of keywords are not computed if below 
#'   \code{min.freq}. Default is 10.
#' @return A list of class \code{KWIC} containing the following components: 
#'   \tabular{ll}{ \code{KWIC} \tab The KWIC index in the form of a data frame. 
#'   \cr \code{KeywordFreq} \tab A data frame of the keywords detected with 
#'   frequency greater than \code{min.freq}. \cr \code{Fields} \tab A character 
#'   vector with the names of the PGR database fields from which the keywords 
#'   were extracted. \cr }
#' @seealso \code{\link[data.table]{data.table}},
#'   \code{\link[PGRdup]{print.KWIC}}
#' @references Knupffer, H. (1988). The European Barley Database of the ECP/GR: 
#'   an introduction. Die Kulturpflanze, 36(1), 135-162.\cr \cr Knupffer, H., 
#'   Frese, L., & Jongen, M. W. M. (1997). Using central crop databases: 
#'   searching for duplicates and gaps. In: E. Lipman, M. W. M. Jongen, T. J. L.
#'   van Hintum, T. Gass, and L. Maggioni (Eds.), Central Crop Databases: Tools 
#'   for Plant Genetic Resources Management (pp. 67-77). Rome, Italy: 
#'   International Plant Genetic Resources Institute/ Wageningen, The 
#'   Netherlands: Centre for Genetic Resources.
#' @examples
#' # Load PGR passport database
#' GN <- GN1000
#' 
#' # Specify as a vector the database fields to be used
#' GNfields <- c("NationalID", "CollNo", "DonorID", "OtherID1", "OtherID2")
#' 
#' # Clean the data
#' GN[GNfields] <- lapply(GN[GNfields], function(x) DataClean(x))
#' 
#' \dontrun{
#' 
#' # Generate KWIC index
#' GNKWIC <- KWIC(GN, GNfields)
#' GNKWIC
#' 
#' # Retrieve the KWIC index from the KWIC object
#' KWIC <- GNKWIC[[1]]
#' 
#' # Retrieve the keyword frequencies from the KWIC object
#' KeywordFreq <- GNKWIC[[2]]
#' 
#' # Show error in case of duplicates and NULL values 
#' # in the primary key/ID field "NationalID"
#' GN[1001:1005,] <- GN[1:5,]
#' GN[1001,3] <- ""
#' GNKWIC <- KWIC(GN, GNfields)
#' }
#' @import data.table
#' @importFrom stringi stri_split_fixed
#' @export KWIC
#' @export print.KWIC
#' @rdname KWIC
KWIC <- function(x, fields, min.freq = 10) {
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
  #setDT(x)
  x <- as.data.table(x)
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
  # Create context fields
  x[, KWIC := do.call(paste, c(.SD, sep = " = ")), .SDcols = fields]
  x[, COMBINED := do.call(paste, .SD), .SDcols = fields]
  # Create KWIC index using data.table
  K <-  as.list(rep(NA, length(fields)))
  for (i in 1:(length(fields))) {
    K[[i]] <-  x[, list(KEYWORD = unlist(strsplit(get(fields[i]), " ")),
                        FIELD = fields[i]),
                 by = list(PRIM_ID = get(fields[1]), KWIC)]
    K[[i]] <- K[[i]][!is.na(K[[i]]$KEYWORD),]
  }
  KWIC <- rbindlist(K)
  rm(K, x)
  #KWIC$KEYWORD[is.na(KWIC$KEYWORD)] <- ""
  set(KWIC, which(is.na(KWIC[["KEYWORD"]])), "KEYWORD", "")
  KWIC <- setkey(KWIC, KEYWORD)
  # Remove all '\' from KWIC
  KWIC[, KWIC := gsub(pattern = "([\\])", replacement = "", x = KWIC)]
  KWIC[, KEYWORD := gsub(pattern = "([\\])", replacement = "", x = KEYWORD)]
  # Remove records with blank keywords
  KWIC <- subset(KWIC, KEYWORD != "")
  # Remove duplicate records
  KWIC <- setkey(KWIC, NULL)
  KWIC <- unique(KWIC)
  # Add padding space in KWIC
  KWIC[, KWIC := paste(" ", KWIC, " ")]
  # Escape all Regex special characters
  KWIC[, KEYWORD := gsub(pattern = "([.|()\\^{}+$*?]|\\[|\\])",
                         replacement = "\\\\\\1", x = KEYWORD)]
  # Highlight keywords in KWIC
  KWIC[, KWIC := mapply(gsub, pattern = paste0(" ", KEYWORD, " "),
                        replacement = paste0(" <<", KEYWORD, ">> "), KWIC)]
  KWIC[, KWIC := gsub("^\\s+|\\s+$", "", KWIC)]
  # Unescape all Regex special characters
  KWIC[, KEYWORD := gsub(pattern = "\\\\(.)", replacement = "\\1", x = KEYWORD)]
  # Split KWIC
  KWIC[, c("KWIC_L", "KW1") := do.call(rbind.data.frame,
                                 stri_split_fixed(KWIC, "<<", 2))][]
  KWIC[, c("KWIC_KW", "KWIC_R") := do.call(rbind.data.frame,
                                stri_split_fixed(KW1, ">>", 2))][]
  cols <- c("KWIC_L", "KWIC_KW", "KWIC_R")
  KWIC[, (cols) := lapply(.SD, as.character), .SDcols = cols]
  KWIC[, KW1 := NULL]
  # Clean output data.frame
  KWIC <- setkey(KWIC, FIELD)
  KWIC <- setkey(KWIC, PRIM_ID)
  setcolorder(KWIC, c("PRIM_ID", "FIELD", "KEYWORD", "KWIC", "KWIC_L",
                      "KWIC_KW", "KWIC_R"))
  KWICIndex <- list(KWIC = NULL, KeywordFreq = NULL, Fields = fields)
  #KWICIndex[[1]] <- as.data.frame(KWIC)
  KWICIndex[[1]] <- setDF(KWIC)
  # Get keyword freq
  kwf <- as.data.frame(table(KWIC$KEYWORD))
  kwf <- subset(kwf, Freq > min.freq)
  kwf <- kwf[order(-kwf$Freq), ]
  rownames(kwf) <- NULL
  setnames(kwf, old = "Var1", new = "Keyword")
  KWICIndex[[2]] <- kwf
  # Set Class
  class(KWICIndex) <- "KWIC"
  return(KWICIndex)
}
